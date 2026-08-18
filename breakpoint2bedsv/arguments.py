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

from breakpoint2bedsv import __version__
import os
import argparse
import tempfile

import logging
logger = logging.getLogger(__name__)


def valid_tool_path(tool_path, tool_name):
    """
    Validate that the given tool is installed and executable.
    - If a full path is provided, it must exist
    - If only a command name is provided, it must be found in $PATH
    - Validate that a CLI tool exists and prints 'usage' or 'help' when run without arguments
    """
        
    # Resolve "full path" / "command name" in PATH
    resolved_path = shutil.which(tool_path)
    logger.debug("Checking executable: %s", resolved_path)

    if resolved_path is None:
        raise ValueError(f"{tool_name} not found in PATH ('{tool_path}').")

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
            raise ValueError(f"{tool_name} does not seem valid ('{tool_path}').")

    except Exception as e:
        raise ValueError(f"Cannot execute {tool_name} ('{tool_path}'). {str(e)}")

    logger.debug("Resolved %s to: %s", tool_name, resolved_path)

    return resolved_path

 

    
def parse_args(argv=None):
    """
    Configure bp2BedSV options from argv.
    """

    # Creation of the parser
    ########################
    parser = argparse.ArgumentParser(
        description="Convert SV breakpoints from VCF/BCF to BED",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Definition of the arguments
    #############################

    # ───────────────────────────────────────────
    # 0) HELP, VERSION & LOGGING
    # (argparse automatically adds -h/--help)
    # ───────────────────────────────────────────
    parser.add_argument(
        "-V", "--version",
        action="version",
        version=f"breakpoint2bedsv {__version__}",
        help="show program's version number and exit"
    )

    parser.add_argument(
        "--log-file",
        metavar="<File>",
        help="write log messages to the specified file"
    )

    # ───────────────────────────────────────────
    # 1) INPUT FILES
    # ───────────────────────────────────────────
    group_input = parser.add_argument_group("Input files")

    group_input.add_argument(
        "-i", "--input-file", dest="input_file",
        metavar="<File>",
        required=True,
        help="""the SV VCF/BCF input file
VCF/VCF.gz/BCF files are supported
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
        help="""output BED file containing non redundant SV breakpoints
(VCF/BCF IDs are merged as a comma-separated list when
multiple variants share the same coordinates)
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
        help="""directory where temporary files will be created
if not provided, the system default temporary directory is used"""
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
    args = parser.parse_args(argv)

        
    # Check tmp_dir
    ###############
    # Determine tmp_dir
    if args.tmp_dir is None:
        args.tmp_dir = tempfile.gettempdir()  # default system tmp
    else:
        # Ensure directory exists
        if not os.path.isdir(args.tmp_dir):
            raise ValueError(
                f"Temporary directory does not exist: {args.tmp_dir}"
            )

    
    # Determine output_dir if not given in argument
    ###############################################
    if args.output_dir is None:
        if "/" in args.output_file:
            output_dir = os.path.dirname(args.output_file)
            if not os.path.exists(output_dir):
                output_dir = "."
        else:
            output_dir = "."
    else:
        output_dir = args.output_dir

    # Store output_dir in args
    args.output_dir = output_dir


    # Determine output_file
    #######################
    if not args.output_file.lower().endswith(".bed"):
        args.output_file += ".bed"


    return args