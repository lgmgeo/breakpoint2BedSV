#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
breakpoint2bedsv
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

import os
import platform
import tempfile
from pathlib import Path
from variant_extractor import VariantExtractor
from breakpoint2bedsv import __version__

import logging
logger = logging.getLogger(__name__)


def run_pipeline(args):

    # Import the different modules
    # (to keep here after the definition of the correct relative path to sys.path)
    ##############################################################################
    from breakpoint2bedsv.workflow import normalize_and_filter_vcf, write_bed, merge_and_sort_bed
    from breakpoint2bedsv.file_utils import ensure_bgzf, has_only_valid_variants



    # Display
    #########
    logger.info("breakpoint2bedsv %s", __version__)
    logger.info("Copyright (C) 2026-current GEOFFROY Veronique")
    logger.info("Please feel free to create a Github issue for any suggestions or bug reports")
    logger.info("https://github.com/lgmgeo/breakpoint2bedsv/issues")
    logger.info("Python version: %s", platform.python_version())
    


    # Arguments display
    ###################
    logger.info("Listing arguments")
    logger.info("           ***************************************************")
    logger.info("           breakpoint2bedsv has been run with these arguments:")
    logger.info("           ***************************************************")

    for key, value in sorted(vars(args).items()):
        if key == "version":
            continue
        if value in ("", None):
            continue

        key = key.replace("_", "-")
        logger.info("           --%s %s", key, value)

    logger.info("           ***************************************************")


    # Ensure that the input SV file is pysam-compatible (VCF/BCF)
    #############################################################
    tmp_bgzf = None
    args.input_file, tmp_bgzf = ensure_bgzf(args.input_file)


    # Check the input_file
    ######################
    try:
        has_only_valid_variants(args.input_file)

    except ValueError:
        raise

    except Exception:
        logger.exception("Unexpected error  while checking the input file")
        raise


    # Normalise ALT (for shorthand notation)
    # e.g. <DUP:SVSIZE=59:AGGREGATED> >> <DUP:AGGREGATED>
    # => creation of the "tmp_normalize_path" normalized VCF
    #####################################################
    logger.info("Normalizing ALT field (for shorthand notation interpretation)")
    tmp_normalized = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
    tmp_normalized_path = tmp_normalized.name
    tmp_normalized.close()
    normalize_and_filter_vcf(args.input_file, tmp_normalized_path)


    # Load and parse the normalized input VCF using VariantExtractor
    ################################################################
    logger.info("Loading the normalized data")
    extractor = VariantExtractor(tmp_normalized_path, ensure_pairs=False)
    os.remove(tmp_normalized_path)


    # Write parsed SV breakpoints to the output BED file
    ####################################################
    logger.info("Parsing SV breakpoints")
    tmp_bed = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    tmp_bed_path = tmp_bed.name
    tmp_bed.close()
    write_bed(extractor, tmp_bed_path)


    # REMOVE REDUNDANCY + SORT BED FINAL
    ####################################
    logger.info("Writing the sorted output BED file")
    merge_and_sort_bed(tmp_bed_path, args.output_file)
    os.remove(tmp_bed_path)
    

    # Cleanup temporary files
    #########################
    if tmp_bgzf is not None:
        Path(tmp_bgzf).unlink(missing_ok=True)
        Path(tmp_bgzf + ".tbi").unlink(missing_ok=True)


    # Finished
    ##########
    logger.info("breakpoint2bedsv completed successfully")




