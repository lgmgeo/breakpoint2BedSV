"""
bp2BedSV 0.1
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
import re
import subprocess
import time
import pysam
import sys


def has_only_valid_variants(sv_file: str) -> bool:
    """
    Check if a VCF/VCF.gz/BCF file:
    - has the good extension
    - contains at least 1 SV
    - exists
    - is valid

    Returns:
        True if file is empty or invalid, False otherwise.
    """

    # Quick check on extension
    if not re.search(r"\.vcf(\.gz)?$|\.bcf$", sv_file, re.IGNORECASE):
        print(f"[WARNING] Not the correct extension: {sv_file}")
        sys.exit(2)

    try:
        with pysam.VariantFile(sv_file) as vf:
            # Try to get first record
            for _ in vf:
                return True  # Found at least one variant

    except FileNotFoundError:
        # File doesn't exist 
        print(f"[WARNING] File doesn't exist: {sv_file}: {e}")
        sys.exit(2)
    except ValueError:
        # pysam raises ValueError for invalid format
        sys.exit(2)
    except Exception as e:
        print(f"[WARNING] Could not read file: {sv_file}: {e}")
        sys.exit(2)

    # No records found → file is empty
    sys.exit(2)




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




def is_multi_allelic(g_bp2BedSV):
    """
    Check if the SV input file contains multi-allelic lines.
    """

    print(f"[{time.strftime('%H:%M:%S')}] Ensuring that the SV input file contains only biallelic variants")

    input_file = g_bp2BedSV["input_file"]
    print(input_file)
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
            print("Please split the multi-allelic lines of the SV input file before to run bp2BedSV")
            print("Exit without error.")
            print("####################################################################################")
            sys.exit(0)
    except subprocess.CalledProcessError:
        pass

