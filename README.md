
<div align="center">
  <h1 style="font-weight: bold; margin-bottom: 0.2em;">breakpoint2BedSV</h1>
  <h3 style="margin-top: 0;">Convert SV breakpoints from VCF/BCF to BED</h3>
</div>

- [Why extracting start/end SV breakpoints from a VCF is not trivial](#why-extracting-startend-sv-breakpoints-from-a-vcf-is-not-trivial)
- [Requirements](#requirements)
- [Quick Installation](#quick-installation)
- [Command line usage / Options](#command-line-usage--options)
- [Outputs](#outputs)
- [How to cite?](#how-to-cite)
- [Tests](#tests)
- [Example use: Assessment of SV presence/absence in a cohort relative to gnomAD v4 SV](#example-use-assessment-of-sv-presenceabsence-in-a-cohort-relative-to-gnomad-v4-sv)

## Why extracting start/end SV breakpoints from a VCF is not trivial
In an SV VCF, the first breakpoint is usually straightforward to retrieve from the `CHROM` and `POS` columns.
However, the second breakpoint is not encoded in a single standardized way and may appear in different fields depending on the SV type or the caller.

| SV representation                                          | First breakpoint | Second breakpoint                        |
| ---------------------------------------------------------- | ---------------- | ---------------------------------------- |
| Symbolic allele (`<DEL>`, `<DUP>`, `<INV>`, `<CNV>`)       | `CHROM:POS`      | usually `INFO/END`                       |
| Breakend notation (e.g. `]chr13:53040041]ATATATATACACACA)` | `CHROM:POS`      | embedded in the `ALT` field              |
| Sequence notation (e.g. `REF=A` and `ALT=ATGATTCGTTCTG...`)| `CHROM:POS`      | embedded in the `REF` field              |
| Sequence notation (e.g. `REF=TGGAATTAGCCTG...` and `ALT=T`)| `CHROM:POS`      | embedded in the `REF` field              |
| Caller-specific representations                            | `CHROM:POS`      | may use alternative tags such as `SVEND` |

As a consequence, extracting both breakpoints from an SV VCF requires handling multiple representations.
`breakpoint2BedSV` addresses this issue by converting heterogeneous SV representations into a unified BED-like breakpoint format.


## Requirements
```
python >=3.8
#poetry #(https://python-poetry.org/docs/#installation)
pysam==0.22.1
variant_extractor==5.1.0
```

## Quick Installation
```
conda create -n breakpoint2BedSV python=3.8 pysam==0.22.1
pip install variant-extractor

conda activate breakpoint2BedSV

# Install with poetry
(not yet available)

# Install from PyPI
pip3 install breakpoint2BedSV (not yet available)

# Install from GitHub
git clone git@github.com:lgmgeo/breakpoint2BedSV.git

# Upgrade 
pip3 install breakpoint2BedSV --upgrade (not yet available)
```

## Command line usage / Options
```bash
usage: breakpoint2BedSV.py [-h] [-V] -i <File> [-d <Dir>] -o <File> [-T <Dir>] [-v]

Convert SV breakpoints from VCF/BCF to BED

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit

Input files:
  -i <File>, --input-file <File>
                        the SV VCF/BCF input file
                        VCF/VCF.gz/BCF files are supported
                        multi-allelic lines are not allowed
                        required

Output options:
  -d <Dir>, --output-dir <Dir>
                        the output directory
                        default: current directory
  -o <File>, --output-file <File>
                        Output BED file with all the SV breakpoints
                        required

Behavior:
  -T <Dir>, --tmp-dir <Dir>
                        Directory where temporary files will be created.
                        If not provided, the system default temporary directory is used.
  -v, --verbose         enable verbose output

```

## Outputs
Running the tool will generate a BED output file with SV start/end coordinates and the associated VCF ID.

## How to cite?
Please cite the following doi if you are using this tool in your research:<br>
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.14922213-blue.svg)](https://doi.org/)

## Tests
A set of tests is included to ensure the correct functioning of breakpoint2BedSV.
The test scripts are located in the `tests/breakpoint2BedSV` directory.

Running them after installation or before submitting changes is recommended to verify that everything works as expected.

## Example use: Assessment of SV presence/absence in a cohort relative to gnomAD v4 SV

**Aim**

Annotate SVs in a VCF with a `gnomAD_excl` flag when at least one breakpoint overlaps a gnomAD v4 SV exclusion region.

**Workflow**

```text
SV VCF
  │
  ├── breakpoint2BedSV
  │      → convert all SVs into breakpoint-level BED intervals
  │
  ├── bedtools intersect
  │      → overlap SV breakpoints with gnomAD v4 SV exclusion regions
  │
  ├── collect overlapping SV IDs
  │
  └── annotate VCF
         → add INFO flag: gnomAD_excl
```

**GRCh38 gnomAD exclusion resources**

SV calling is less reliable in some genomic regions due to:
- low mappability / depth bias
- peri-centromeric or peri-telomeric repeats
- known problematic regions in population datasets such as gnomAD

Two GRCh38 gnomAD exclusion regions:
* `depth_blacklist.sorted.bed.gz`
* `PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed.gz`
```bash
curl -O https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/depth_blacklist.sorted.bed.gz
curl -O  https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed.gz
```

**Output**

SVs with at least one breakpoint overlapping one of these exclusion regions are flagged in the VCF with:

```vcf
##INFO=<ID=gnomAD_excl,Number=0,Type=Flag,Description="At least one SV breakpoint overlaps a gnomAD exclusion region">
```

**Implementation**

 1. Convert SV VCF to breakpoint BED format

```bash
breakpoint2BedSV \
  --vcf input.vcf \
  --output sv.breakpoints.bed
```

This step standardizes all SV representations (DEL/DUP/INV/BND/SVEND) into a unified breakpoint BED format.

---

 2. Identify SVs overlapping gnomAD v4 exclusion regions

```bash
bedtools intersect \
  -a sv.breakpoints.bed \
  -b depth_blacklist.sorted.bed.gz \
  -wa | cut -f4 | sort -u > excluded_ids.txt

bedtools intersect \
  -a sv.breakpoints.bed \
  -b PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed.gz \
  -wa | cut -f4 | sort -u >> excluded_ids.txt

sort -u excluded_ids.txt > excluded_ids.final.txt
```

---

 3. Annotate original VCF with `gnomAD_excl` flag

```bash
awk -F'\t' '
BEGIN {
    OFS="\t"
    while ((getline line < "excluded_ids.final.txt") > 0)
        excl[line] = 1
}
{
    if ($0 ~ /^#/) {
        print
        next
    }

    id = $3

    if (id in excl) {
        if ($8 == "." || $8 == "") {
            $8 = "gnomAD_excl"
        } else {
            $8 = $8 ";gnomAD_excl"
        }
    }

    print
}
' input.vcf > input.gnomAD_excl.vcf
```
