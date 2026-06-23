"""
bp2BedSV 0.1
============

Copyright (C) 2026-current Veronique Geoffroy (veronique.geoffroy@inserm.fr)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; If not, see <http://www.gnu.org/licenses/>.
"""
import subprocess
import sys
import gzip
from pathlib import Path

def open_variant_stream(path):
    """
    Open a variant file as a text stream.

    Supports:
      - .vcf
      - .vcf.gz
      - .bcf

    Returns
    -------
    fin : file-like object
        Readable text stream containing VCF lines.
    proc : subprocess.Popen or None
        bcftools process if input is BCF, otherwise None.
    """
    path = Path(path)
    suffixes = path.suffixes

    if suffixes and suffixes[-1] == ".bcf":
        proc = subprocess.Popen(
            ["bcftools", "view", "-Ov", str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return proc.stdout, proc

    elif len(suffixes) >= 2 and suffixes[-2:] == [".vcf", ".gz"]:
        return gzip.open(path, "rt"), None

    elif suffixes and suffixes[-1] == ".vcf":
        return open(path, "r"), None

    else:
        raise ValueError(
            f"Unsupported input format for {path}. Expected .vcf, .vcf.gz or .bcf"
        )

def normalize_shorthand_notation_in_alt(svfile_in, svfile_out, chunk_size=50000):
    """
    Normalize ALT fields for variant-extractor compatibility while preserving non-key metadata tags.
    Specifically removes tags with key=value attributes (e.g. SVSIZE=59) but keeps standalone tags (e.g. AGGREGATED).
    Example:
        <DUP:SVSIZE=59:AGGREGATED> >> <DUP:AGGREGATED>
    Supports:
        - .vcf
        - .vcf.gz
        - .bcf
    """

    def fix_alt(alt):
        if not (alt.startswith("<") and alt.endswith(">")):
            return alt

        content = alt[1:-1]
        parts = content.split(":")
        cleaned = [p for p in parts if "=" not in p]

        return "<" + ":".join(cleaned) + ">"


    buffer = []
    fin, proc = open_variant_stream(svfile_in)

    try:
        with fin, open(svfile_out, "w") as fout:
            for line in fin:
                if line.startswith("#"):
                    buffer.append(line)
                    continue

                fields = line.rstrip("\n").split("\t")

                # ligne non standard / vide : on la recopie telle quelle
                if len(fields) < 5:
                    buffer.append(line)
                    continue

                fields[4] = fix_alt(fields[4])
                buffer.append("\t".join(fields) + "\n")

                if len(buffer) >= chunk_size:
                    fout.writelines(buffer)
                    buffer = []

            if buffer:
                fout.writelines(buffer)

        # si on est passé par bcftools, vérifier qu'il s'est terminé correctement
        if proc is not None:
            stderr = proc.stderr.read()
            retcode = proc.wait()
            if retcode != 0:
                raise RuntimeError(
                    f"bcftools failed on {svfile_in} (exit code {retcode}):\n{stderr}"
                )

    finally:
        if proc is not None and proc.poll() is None:
            proc.terminate()




def write_bed(extractor, out_path, chunk_size=5000):

    buffer = []

    def flush(out):
        if buffer:
            out.writelines(buffer)
            buffer.clear()

    with open(out_path, "w") as out:

        for sv in extractor:

            chrom1 = sv.contig
            pos1 = sv.pos
            svid = sv.id or "."

            # -------------------------
            # BND / TRANSLOCATION FIRST
            # -------------------------
            is_bnd = False

            if hasattr(sv, "mate_contig") and hasattr(sv, "mate_pos"):
                chrom2 = sv.mate_contig
                pos2 = sv.mate_pos
                is_bnd = True

            elif hasattr(sv, "alt") and sv.alt and ":" in str(sv.alt):
                import re
                m = re.search(r'([^:\[\]]+):(\d+)', str(sv.alt))
                if m:
                    chrom2 = m.group(1)
                    pos2 = int(m.group(2))
                    is_bnd = True

            if is_bnd:
                buffer.append(f"{chrom1}\t{pos1-1}\t{pos1}\t{svid}\n")
                buffer.append(f"{chrom2}\t{pos2-1}\t{pos2}\t{svid}\n")

                if len(buffer) >= chunk_size:
                    flush(out)

                continue

            # -------------------------
            # SV NON-BND (DEL/DUP/INV)
            # -------------------------
            buffer.append(f"{chrom1}\t{pos1-1}\t{pos1}\t{svid}\n")

            if getattr(sv, "end", None) and sv.end != pos1:
                buffer.append(f"{chrom1}\t{sv.end-1}\t{sv.end}\t{svid}\n")

            if len(buffer) >= chunk_size:
                flush(out)

        # flush final
        flush(out)



