#!/usr/bin/env python3

import argparse
import re
import tempfile
from variant_extractor import VariantExtractor



def normalize_shorthand_notation_in_alt(vcf_in, vcf_out, chunk_size=50000):
    """
    Normalize ALT fields for variant-extractor compatibility while preserving non-key metadata tags.

    Specifically removes key=value attributes (e.g. SVSIZE=59) but keeps standalone tags (e.g. AGGREGATED).

    Example:
        <DUP:SVSIZE=59:AGGREGATED> >> <DUP:AGGREGATED>
    """
    def fix_alt(alt):
        if not (alt.startswith("<") and alt.endswith(">")):
            return alt

        content = alt[1:-1]
        parts = content.split(":")

        cleaned = []
        for p in parts:
            if "=" in p:
                continue
            cleaned.append(p)

        return "<" + ":".join(cleaned) + ">"

    buffer = []

    with open(vcf_in) as fin, open(vcf_out, "w") as fout:

        for line in fin:

            if line.startswith("#"):
                buffer.append(line)
                continue

            fields = line.rstrip().split("\t")
            fields[4] = fix_alt(fields[4])

            buffer.append("\t".join(fields) + "\n")

            # flush par batch
            if len(buffer) >= chunk_size:
                fout.writelines(buffer)
                buffer = []

        # flush final
        if buffer:
            fout.writelines(buffer)



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



def main():

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Consvert Ssv svCF to BED breakpoints using variant-extractor"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input svCF file"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output BED file with all the Ssv breakpoints"
    )

    args = parser.parse_args()

    # Normalise ALT (for shorthand notation)
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
    tmp_path = tmp.name
    tmp.close()
    normalize_shorthand_notation_in_alt(args.input, tmp_path)
    
    # Load and parse the input VCF using VariantExtractor.
    extractor = VariantExtractor(tmp_path, ensure_pairs=False)


    # Write parsed SV breakpoints to the output BED file. 
    tmp_bed = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    tmp_bed_path = tmp_bed.name
    tmp_bed.close()

    write_bed(extractor, tmp_bed_path)

    # SORT BED FINAL
    with open(tmp_bed_path) as f:
        lines = f.readlines()

    lines.sort(key=lambda x: (x.split("\t")[0], int(x.split("\t")[1])))

    with open(args.output, "w") as out:
        out.writelines(lines)




if __name__ == "__main__":
    main()
