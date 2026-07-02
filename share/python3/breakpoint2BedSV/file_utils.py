"""
breakpoint2BedSV 1.0
====================

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
import shutil
import tempfile
from pathlib import Path
import re
import subprocess
import time
import pysam
import sys


def ensure_bgzf(path):
    """
    Ensure that a .vcf.gz file is BGZF-compressed.

    If the input file is already BGZF-compressed and readable by
    ``pysam.VariantFile``, its path is returned unchanged.

    Otherwise, the file is decompressed and recompressed into a temporary
    BGZF-compressed VCF, which is indexed with tabix.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the input VCF/VCF.gz file.

    Returns
    -------
    tuple[str, str | None]
        - Path to a BGZF-compatible VCF.
        - Path to the temporary BGZF file if one was created, otherwise None.

    Notes
    -----
    When a temporary BGZF file is created, the caller is responsible for
    removing both the ``.vcf.gz`` file and its associated ``.tbi`` index.
    """
    path = Path(path)

    # If the file is not a .vcf.gz, return it as is
    if not str(path).endswith(".vcf.gz"):
        return str(path), None

    try:
        # Check if the .vcf.gz is already BGZF-compressed
        pysam.VariantFile(str(path))
        return str(path), None
    
    except Exception:
        # If not, decompress and recompress to BGZF in a temporary file

        # Create a temporary file for the decompressed VCF
        tmp_vcf = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False)
        tmp_vcf.close()

        try:
            # Decompress the gzip file to the temporary VCF file
            with gzip.open(path, "rb") as fin, open(tmp_vcf.name, "wb") as fout:
                shutil.copyfileobj(fin, fout)

            # Create a temporary file for the BGZF-compressed VCF
            tmp_bgzf = tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False)
            tmp_bgzf.close()

            # Recompress the temporary VCF file to BGZF
            pysam.tabix_compress(tmp_vcf.name, tmp_bgzf.name, force=True)
            pysam.tabix_index(tmp_bgzf.name, preset="vcf", force=True)

            return tmp_bgzf.name, tmp_bgzf.name
        finally:
            # Clean up the temporary decompressed VCF file
            Path(tmp_vcf.name).unlink(missing_ok=True)



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

    print(f"[{time.strftime('%H:%M:%S')}] Ensuring that the SV input file contains only valid variants")

    # Quick check on extension
    if not re.search(r"\.vcf(\.gz)?$|\.bcf$", sv_file, re.IGNORECASE):
        raise ValueError(f"[ERROR] Not the correct extension: {sv_file}")

    try:

        with pysam.VariantFile(sv_file) as vf:
            # Try to get first record
            for _ in vf:
                return True  # Found at least one variant
    except FileNotFoundError:
        # File doesn't exist 
        raise ValueError(f"[ERROR] File doesn't exist: {sv_file}")
    except (ValueError, OSError) as e:
        # pysam raises ValueError for invalid format
        raise ValueError(f"[ERROR] Invalid VCF/BCF file {sv_file}: {e}")
    except Exception as e:
        raise ValueError(f"[ERROR] Could not read file {sv_file}: {e}")

    # No records found → file is empty
    raise ValueError(f"[ERROR] No SV found in file: {sv_file}")




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
