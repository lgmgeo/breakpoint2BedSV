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

import sys
import logging

from breakpoint2bedsv.arguments import parse_args
from breakpoint2bedsv.logging_utils import setup_logging
from breakpoint2bedsv.core import run_pipeline

# Module-level logger.
# The logging configuration is initialized in main() after parsing
# the command-line arguments.
logger = logging.getLogger(__name__)

def main(argv=None):
    """
    Main entry point for the breakpoint2bedsv command-line interface.

    Parameters
    ----------
    argv : list[str] or None
        Command-line arguments.
        If None, arguments are read from sys.argv.

    Returns
    -------
    int
        Exit status:
        - 0: successful execution
        - 1: unexpected error
        - 2: input or argument validation error
    """

    # Use the arguments provided by the caller.
    # When called from the command line, use sys.argv instead.
    if argv is None:
        argv = sys.argv[1:]

    
    try:

        # Parse and validate command-line arguments.
        # This also sets default values for optional arguments.
        args = parse_args(argv)

        # Configure the logging system according to the command-line
        # options (--verbose and --log-file).
        setup_logging(
            verbose=args.verbose,
            log_file=args.log_file
        )

        # Execute the main breakpoint2bedsv processing pipeline.
        # Errors are allowed to propagate and are handled below.
        run_pipeline(args)

    except ValueError as e:
        # Input or configuration error.
        # Return a dedicated exit code so that the caller can
        # distinguish validation errors from unexpected failures.
        logger.error("%s", e)
        return 2

    except Exception:
        # Catch any unexpected error.
        # logger.exception() automatically includes the traceback,
        # which is useful for debugging.
        logger.exception("Unexpected error")
        return 1

    # The pipeline completed successfully.
    return 0


if __name__ == "__main__":
    sys.exit(main())