"""
liftoverSV 0.1
==============

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

import gzip
import functools
import builtins
import re
from itertools import islice
import os
import subprocess
import time


print_flush = functools.partial(builtins.print, flush=True)


def is_an_empty_vcf_file(vcf_file: str) -> bool:
    """
    Check if a VCF file (.vcf or .vcf.gz) is empty (ignoring header lines),
    or if it isn't a VCF file based on its extension.

    Returns:
        True if the file is empty or not a VCF, False otherwise.
    """
    # Quick check on the file extension
    if not re.search(r"\.vcf(\.gz)?$", vcf_file, re.IGNORECASE):
        return True

    # Determine whether to use gzip or normal open
    open_func = gzip.open if vcf_file.endswith(".gz") else open

    try:
        with open_func(vcf_file, "rt") as f:  # 'rt' = text mode
            for line in f:
                line = line.strip()
                # Ignore empty lines and header lines
                if line and not line.startswith("#"):
                    return False  # Found at least one data line
    except FileNotFoundError:
        # File doesn't exist → treat as empty
        return True
    except Exception as e:
        # Other read errors → treat as non-VCF/empty
        print(f"[WARNING] Could not read file {vcf_file}: {e}")
        return True

    # No data lines found → file is empty
    return True



def check_vcf_variant_line_format(vcf_path):
    """
    Check:
     - whether the #CHROM line exist
     - whether variant lines in a VCF file have the same number of fields
        as defined by the #CHROM header line.
     - wether empty variant line exists 

    Returns:
      "OK" if the VCF is well-formatted; otherwise, returns the error message.
    """
    
    chrom_header = None
    expected_field_number = None

    with open_any_text_file(vcf_path) as f:
        for line_number, line in enumerate(f, start=1):

            # Skip metadata
            if line.startswith("##"):
                continue

            # Extract #CHROM line
            if line.startswith("#CHROM"):
                chrom_header = line.rstrip("\n")
                expected_field_number = len(chrom_header.split("\t"))
                continue

            # Skip if #CHROM not found yet
            if chrom_header is None:
                return "Header line starting with #CHROM not found in the VCF!"

            # Stop if we reach an empty line
            if not line.strip():
                return "Empty line present in the VCF"

            # Process variant line
            fields = line.rstrip("\n").split("\t")
            num_fields = len(fields)

            if num_fields != expected_field_number:
                return f"Incorrect number of fields: {expected_field_number} expected, {num_fields} on line {line_number}"
                
    return "OK"




# Natural sorting 
#################
# chromosomes = ["chr1", "chr10", "chr2", "chrX"]
# sorted_chromosomes = sorted(chromosomes, key=natural_sort_key)
# print(sorted_chromosomes)
# => Output: ['chr1', 'chr2', 'chr10', 'chrX']
def natural_sort_key(s):
    """Return a key for natural sorting (like Tcl -dictionary)."""
    # split into list of ints and non-ints
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]



def open_any_text_file(path):
    """
    Open a text file transparently, whether it is plain text or gzip-compressed.

    This function automatically detects gzip compression by checking the magic
    bytes (0x1F, 0x8B), so it works regardless of the file extension. It is intended
    for reading text-based files (VCF, BED, TSV, CSV, config files, etc.) that may
    be distributed either uncompressed or compressed with gzip.

    Args:
        path (str): Path to the file to open.

    Returns:
        file object: A text-mode file handle. 
        If the file is gzip-compressed, a gzip.open() handle is returned; otherwise, a standard open().

    Raises:
        OSError: If the file cannot be opened or read.
    """
    # Detect gzip by reading the first two bytes of the file
    with open(path, "rb") as test:
        start = test.read(2)

    # Magic bytes for gzip compression = 1F 8B
    if start == b"\x1f\x8b":
        return gzip.open(path, "rt")   # gzip-compressed text file
    else:
        return open(path, "rt")        # plain-text file




def is_multi_allelic(g_liftoverSV):
    """
    Check if the VCF input file contains multi-allelic lines.
    """

    print(f"[{time.strftime('%H:%M:%S')}] Ensuring that the input VCF contains only biallelic variants")

    input_file = g_liftoverSV["input_file"]

    if input_file.endswith(".gz"):
        cmd = f"zcat {input_file} | grep -v ^# | cut -f 4-5 | grep -c ,"
    else:
        cmd = f"grep -v ^# {input_file} | cut -f 4-5 | grep -c ,"

    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        output = result.stdout.strip()
        n_multiallelic_line = output[0] if output else "" 
        if n_multiallelic_line != "0":
            print("####################################################################################")
            print("Please split the multi-allelic lines of the VCF input file before to run liftoverSV")
            print("Exit without error.")
            print("####################################################################################")
            sys.exit(0)
    except subprocess.CalledProcessError:
        pass

