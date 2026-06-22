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

import os
import sys
import re
import subprocess
import shutil
import argparse
import tempfile
from file_utils import is_an_empty_vcf_file, check_vcf_variant_line_format, print_flush as print
from functools import partial

def valid_vcf_input_file(vcf_input_file):
    """
    Check if the VCF input file:
    - exists
    - has .vcf or .vcf.gz extension
    - contains at least 1 SV line (not empty)
    - contains the #CHROM line
    - has the same number of fields in the variant lines as defined by the #CHROM header line.
    - empty variant line exists 
    """
    # Check if the file exists
    if not os.path.isfile(vcf_input_file):
        # Normal error -> argparse will handle and exit with code 2
        raise argparse.ArgumentTypeError(f"File '{vcf_input_file}' does not exist")

    # Check if the file has a valid VCF extension
    if not (vcf_input_file.endswith(".vcf") or vcf_input_file.endswith(".vcf.gz")):
        raise argparse.ArgumentTypeError(f"File '{vcf_input_file}' must have .vcf or .vcf.gz extension")

    # Special case: empty VCF file -> print message and exit with code 0
    if is_an_empty_vcf_file(vcf_input_file):
        print(f"File '{vcf_input_file}' is empty, no SV to lift")
        sys.exit(0)

    # Check if the input VCF file contains multi-allelic lines
    message = check_vcf_variant_line_format(vcf_input_file)
    if message != "OK":
        print(message)
        sys.exit(1)

    # If all checks pass, return the vcf_input_file
    return vcf_input_file
    



def valid_tool_path(tool_path, tool_name):
    """
    Validate that the given tool is installed and executable.
    - If a full path is provided, it must exist
    - If only a command name is provided, it must be found in $PATH
    - Validate that a CLI tool exists and prints 'usage' or 'help' when run without arguments
    """
        
    # Resolve "full path" / "command name" in PATH
    resolved_path = shutil.which(tool_path)
    if resolved_path is None:
        print(f"\nError: {tool_name} not found in PATH ('{tool_path}').")
        sys.exit(2)

    # Try running the tool
    try:
        # Run the tool without arguments, capture stdout and stderr
        result = subprocess.run(
            [resolved_path],
            stdout=subprocess.PIPE,  # capture stdout
            stderr=subprocess.STDOUT,  # redirect stderr to stdout
            text=True,                # return string instead of bytes
            timeout=5                  # optional: avoid hanging
        )
 
        # Check if 'usage' or 'help' appears in output
        output = result.stdout.lower()
        if "usage" not in output and "help" not in output:
            print(f"\nError: {tool_name} does not seem valid ('{tool_path}').")
            sys.exit(2)

    except Exception:
        print(f"\nError: Cannot execute {tool_name} ('{tool_path}'). {str(e)}")
        sys.exit(2)

    return resolved_path

 

    
def configure_bp2BedSV(argv, g_bp2BedSV):
    """
    Configure bp2BedSV options from argv.
    """

    # Creation of the parser
    ########################
    parser = argparse.ArgumentParser(
        description="Convert SV breakpoints from VCF to BED using variant-extractor",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Definition of the arguments
    #############################

    # ───────────────────────────────────────────
    # 0) HELP & VERSION
    # (argparse automatically adds -h/--help)
    # ───────────────────────────────────────────
    parser.add_argument(
        "-V", "--version",
        action="version",
        version=f"bp2BedSV {g_bp2BedSV['version']}",
        help="show program's version number and exit"
    )    

    # ───────────────────────────────────────────
    # 1) INPUT FILES
    # ───────────────────────────────────────────
    group_input = parser.add_argument_group("Input files")

    group_input.add_argument(
        "-i", "--input-file", dest="input_file",
        metavar="<File>",
        required=True,
        type=valid_vcf_input_file,
        help="""the SV VCF input file
gzipped VCF file is supported
multi-allelic lines are not allowed
required"""
    )

    # ───────────────────────────────────────────
    # 2) OUTPUT OPTIONS
    # ───────────────────────────────────────────
    group_output = parser.add_argument_group("Output options")

    group_output.add_argument(
        "-d", "--output-dir", dest="output_dir",
        type=str, 
        metavar="<Dir>",
        help="""the output directory
default: current directory"""
    )

    group_output.add_argument(
        "-o", "--output-file", dest="output_file",
        required=True,
        metavar="<File>",
        help="""Output BED file with all the SV breakpoints
required"""
    )

    # ───────────────────────────────────────────
    # 3) BEHAVIORAL PARAMETERS
    # ───────────────────────────────────────────
    group_behavior = parser.add_argument_group("Behavior")



    group_behavior.add_argument(
        "-T", "--tmp-dir", dest="tmp_dir",
        type=str,
        metavar="<Dir>",
        help="""Directory where temporary files will be created.
If not provided, the system default temporary directory is used."""
    )

    group_behavior.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="enable verbose output"
    )

    # Parsing of the arguments
    ##########################
    # WARNING:
    # Argparse converts option names to valid Python identifiers:
    # - long names become lowercase
    # - hyphens '-' are replaced with underscores '_'
    # Access the value via args.option_name, e.g., input-file >> args.input_file
    args = parser.parse_args()

    # Completion of the g_bp2BedSV dictionary
    ###########################################
    g_bp2BedSV.update(vars(args))

    # Check tmp_dir
    ###############
    # Determine tmp_dir
    if g_bp2BedSV["tmp_dir"] is None:
        g_bp2BedSV["tmp_dir"] = tempfile.gettempdir()  # default system tmp
    else:
        # Ensure directory exists
        if not os.path.isdir(g_bp2BedSV["tmp_dir"]):
            raise ValueError(f"Temporary directory does not exist: {g_bp2BedSV['tmp_dir']}")
    
    # Determine output_dir if not given in argument
    ###############################################
    if g_bp2BedSV["output_dir"] == None:
        if "/" in g_bp2BedSV["output_file"]:
            output_dir = os.path.dirname(g_bp2BedSV["output_file"])
            if not os.path.exists(output_dir):
                output_dir = "."    
        else:
            output_dir = "."
    # Store output_dir in global dictionary
    g_bp2BedSV["output_dir"] = output_dir

    # Determine output_file
    #######################
    if not g_bp2BedSV["output_file"].lower().endswith(".bed"):
        g_bp2BedSV["output_file"] += ".bed"


