#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
breakpoint2BedSV 0.1
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

import os
import sys
import re
import time
import platform
import tempfile
from variant_extractor import VariantExtractor



def get_script_directory():
    return(os.path.dirname(os.path.realpath(sys.argv[0])))


def main(argv):

    # Initialisation of the g_bp2BedSV dictionary
    #####################################################
    g_bp2BedSV = {}


    # FHS directories setup:
    # - g_bp2BedSV["install_dir"]
    # - g_bp2BedSV["doc_dir"]
    # - g_bp2BedSV["python_dir"]
    # - g_bp2BedSV["bash_dir"]
    ######################################
    g_bp2BedSV["install_dir"] = get_script_directory()
    # $breakpoint2BedSV/share/python3/breakpoint2BedSV/bin/breakpoint2BedSV.py
    if g_bp2BedSV["install_dir"].endswith("share/python3/breakpoint2BedSV/bin"):
        g_bp2BedSV["install_dir"] = g_bp2BedSV["install_dir"][:-34]
    # $breakpoint2BedSV/bin/breakpoint2BedSV.py
    if g_bp2BedSV["install_dir"].endswith("bin"):
        g_bp2BedSV["install_dir"] = g_bp2BedSV["install_dir"][:-4]

    g_bp2BedSV["doc_dir"]  = os.path.join(g_bp2BedSV["install_dir"], "share", "doc", "breakpoint2BedSV")
    g_bp2BedSV["python_dir"] = os.path.join(g_bp2BedSV["install_dir"], "share", "python3", "breakpoint2BedSV")
    g_bp2BedSV["bash_dir"] = os.path.join(g_bp2BedSV["install_dir"], "share", "bash", "breakpoint2BedSV")


    # Add the correct relative path to sys.path so that Python can locate breakpoint2BedSV_*.py
    # The path is computed relative to the location of this script (__file__), ensuring
    # it works regardless of the current working directory from which the script is executed.
    #########################################################################################
    sys.path.insert(0, os.path.join(g_bp2BedSV['python_dir']))


    # Import the different modules
    # (to keep here after the definition of the correct relative path to sys.path)
    ##############################################################################
    from config import configure_bp2BedSV
    from workflow import normalize_shorthand_notation_in_alt, write_bed


    # Search for the breakpoint2BedSV VERSION
    #########################################
    if "version" not in g_bp2BedSV:
        runFile = os.path.join(g_bp2BedSV["install_dir"], "bin", "breakpoint2BedSV.py")
        if os.path.exists(runFile):
            try:
                # Open the file in text mode
                with open(runFile, "rt") as f:
                    for line in f:
                        m = re.match(r"^breakpoint2BedSV ([0-9]+\.[0-9]+)", line)
                        if m:
                            g_bp2BedSV["version"] = m.group(1)
                            break
            except Exception as e:
                print(f"[WARNING] Could not read {runFile}: {e}")

    if "version" not in g_bp2BedSV:
        g_bp2BedSV["version"] = "X.X"


    # Downloading configuration
    ###########################
    configure_bp2BedSV(argv, g_bp2BedSV)


    # Display
	#########
    print(f"\nbreakpoint2BedSV {g_bp2BedSV['version']}")
    print("Copyright (C) 2026-current GEOFFROY Veronique")
    print("Please feel free to create a Github issue for any suggestions or bug reports")
    print("https://github.com/lgmgeo/breakpoint2BedSV/issues\n")
    print("\nPython version:", platform.python_version(), "\n")
    print("Application name used:")
    print(g_bp2BedSV["install_dir"], "\n")


    # Arguments display
    ###################
    print(f"\n[{time.strftime('%H:%M:%S')}] Listing arguments")
    print("           ***************************************************")
    print("           breakpoint2BedSV has been run with these arguments:")
    print("           ***************************************************")

    for key in sorted(g_bp2BedSV.keys()):
        if key in ["bash_dir", "doc_dir", "install_dir", "python_dir", "tcl_dir", "version"]:
            continue
        val = g_bp2BedSV[key]
        if val == "":
            continue
        key = key.replace("_", "-")
        print(f"           --{key} {val}")

    print("           ***************************************************")



    # Normalise ALT (for shorthand notation)
    # e.g. <DUP:SVSIZE=59:AGGREGATED> >> <DUP:AGGREGATED>
    # => creation of the "tmp_normalize_path" normalized VCF
    #####################################################
    print(f"[{time.strftime('%H:%M:%S')}] Normalizing ALT field (for shorthand notation interpretation)")
    tmp_normalized = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
    tmp_normalized_path = tmp_normalized.name
    tmp_normalized.close()
    normalize_shorthand_notation_in_alt(g_bp2BedSV["input_file"], tmp_normalized_path)


    # Load and parse the normalized input VCF using VariantExtractor
    ################################################################
    print(f"[{time.strftime('%H:%M:%S')}] Loading the normalized data")
    extractor = VariantExtractor(tmp_normalized_path, ensure_pairs=False)
    os.remove(tmp_normalized_path)


    # Write parsed SV breakpoints to the output BED file
    ####################################################
    print(f"[{time.strftime('%H:%M:%S')}] Parsing SV breakpoints")
    tmp_bed = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    tmp_bed_path = tmp_bed.name
    tmp_bed.close()
    write_bed(extractor, tmp_bed_path)


    # SORT BED FINAL
    ################
    print(f"[{time.strftime('%H:%M:%S')}] Writing the sorted output BED file")
    with open(tmp_bed_path) as f:
        lines = f.readlines()

    os.remove(tmp_bed_path)

    lines.sort(key=lambda x: (x.split("\t")[0], int(x.split("\t")[1])))

    with open(g_bp2BedSV["output_file"], "w") as out:
        out.writelines(lines)




	# Finished
	##########
    print(f"[{time.strftime('%H:%M:%S')}] breakpoint2BedSV completed successfully")


if __name__ == "__main__":
    main(sys.argv[1:])



