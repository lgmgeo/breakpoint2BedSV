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
import pysam
from pathlib import Path
import tempfile


def open_variant_stream(path):
    """
    Open a variant file as a text stream.

    Supports:
      - .vcf
      - .vcf.gz (gzip or bgzip)
      - .bcf

    If a .vcf.gz is not BGZF-compressed, it is transparently converted
    to a temporary BGZF file.

    Returns
    -------
    fin : file-like object
        Readable text stream containing VCF lines.
    proc : subprocess.Popen or None
        bcftools process if input is BCF, otherwise None.
    """
    path = Path(path)
    suffixes = path.suffixes

    # BCF 
    if suffixes and suffixes[-1] == ".bcf":
        proc = subprocess.Popen(
            ["bcftools", "view", "-Ov", str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return proc.stdout, proc

    # VCF.GZ
    elif len(suffixes) >= 2 and suffixes[-2:] == [".vcf", ".gz"]:
        # Vérifie si le fichier est déjà BGZF
        try:
            pysam.VariantFile(str(path)).close()
            bgzf_path = path
        except (OSError, ValueError):
            # Conversion automatique en BGZF
            tmp = tempfile.NamedTemporaryFile(
                suffix=".vcf.gz",
                delete=False,
            )
            tmp.close()

            pysam.tabix_compress(str(path), tmp.name, force=True)
            pysam.tabix_index(tmp.name, preset="vcf", force=True)

            bgzf_path = Path(tmp.name)

        proc = subprocess.Popen(
            ["bcftools", "view", "-Ov", str(bgzf_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return proc.stdout, proc        

    # Plain VCF
    elif suffixes and suffixes[-1] == ".vcf":
        return open(path, "r"), None

    else:
        raise ValueError(f"[ERROR] Unsupported input format for {path}. Expected .vcf, .vcf.gz or .bcf")

def normalize_shorthand_notation_in_alt(svfile_in, svfile_out, chunk_size=50000):
    """
    Normalize <ALT> fields for variant-extractor compatibility:
    Removes tags with key=value attributes (e.g. SVSIZE=59) but keeps standalone tags (e.g. AGGREGATED).
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
    n_mcnv = 0
    n_cpx = 0
    n_ctx = 0
    n_bnd = 0

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

                # FILTER: MULTIALLELIC (as in gnomAD SV v4: same ID and ALT = <CN0> or <CN1> or <CN2> or <CN3> or...)
                if "MULTIALLELIC" in fields[6]:
                    n_mcnv += 1
                    continue

                # TYPE: BND/CPX/CTX (as in gnomAD SV v4)
                if fields[4] == "<BND>":
                    n_bnd += 1
                    continue
                if fields[4] == "<CPX>":
                    n_cpx += 1
                    continue  
                if fields[4] == "<CTX>":
                    n_ctx += 1
                    continue

                # ALT: normalisation de la notation abrégée
                fields[4] = fix_alt(fields[4])

                buffer.append("\t".join(fields) + "\n")

                if len(buffer) >= chunk_size:
                    fout.writelines(buffer)
                    buffer = []

            if buffer:
                fout.writelines(buffer)
            if n_mcnv > 0:
                print(f"[{time.strftime('%H:%M:%S')}] - {n_mcnv} MULTIALLELIC lines were removed from the input file")
            if n_cpx > 0:   
                print(f"[{time.strftime('%H:%M:%S')}] - {n_cpx} CPX lines were removed from the input file")
            if n_bnd > 0:
                print(f"[{time.strftime('%H:%M:%S')}] - {n_bnd} BND lines were removed from the input file")
            if n_ctx > 0:
                print(f"[{time.strftime('%H:%M:%S')}] - {n_ctx} CTX lines were removed from the input file")

        # si on est passé par bcftools, vérifier qu'il s'est terminé correctement
        if proc is not None:
            stderr = proc.stderr.read()
            retcode = proc.wait()
            if retcode != 0:
                raise RuntimeError(f"[ERROR] bcftools failed on {svfile_in} (exit code {retcode}):\n{stderr}")

    finally:
        if proc is not None and proc.poll() is None:
            proc.terminate()




def write_bed(extractor, out_path, chunk_size=5000):

    buffer = []

    def flush(out):
        if buffer:
            out.writelines(buffer)
            buffer.clear()

    # Suppress repeated htslib/pysam warnings (e.g. contig/header issues)
    # These warnings can be emitted multiple times because the VCF is parsed previously by pysam (def has_only_valid_variants) and here by VariantExtractor 
    pysam.set_verbosity(0)

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



from collections import defaultdict
from pathlib import Path
import time


def merge_and_sort_bed(input_bed_path: str, output_bed_path: str) -> None:
    """
    Merge identical genomic coordinates (chr, start, end) and sort the BED file.
    For duplicated intervals, the 4th column is merged using comma separation.

    Example
    -------
    Input:
        chr1  100  200  SV1
        chr1  100  200  SV2
        chr1  300  400  SV3

    Output:
        chr1  100  200  SV1,SV2
        chr1  300  400  SV3
    """

    print(f"[{time.strftime('%H:%M:%S')}] Merging and sorting BED file")

    bed_dict = defaultdict(set)

    # Parse input BED file
    with open(input_bed_path, "r") as f:
        for line in f:
            if not line.strip():
                continue

            cols = line.rstrip().split("\t")

            # Key = genomic coordinates
            key = (cols[0], int(cols[1]), int(cols[2]))

            # Value = annotation (col 4 if exists)
            value = cols[3] if len(cols) > 3 else "."

            bed_dict[key].add(value)

    # Build merged BED lines
    merged_lines = [
        f"{chrom}\t{start}\t{end}\t{','.join(sorted(values))}\n"
        for (chrom, start, end), values in bed_dict.items()
    ]

    # Sort final BED
    merged_lines.sort(key=lambda x: (x.split("\t")[0], int(x.split("\t")[1])))

    # Write output
    with open(output_bed_path, "w") as out:
        out.writelines(merged_lines)

