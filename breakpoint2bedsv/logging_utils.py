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

import logging
import sys


def setup_logging(verbose=False, log_file=None):
    """
    Configure the application logging system.

    Parameters
    ----------
    verbose : bool
        If True, display detailed log messages on stderr.
        Otherwise, only informational messages are displayed.

    log_file : str or None
        Optional path to a log file.
        If provided, log messages are also written to this file.
    """

    # Set the logging level according to the verbose option.
    # DEBUG messages are displayed when verbose mode is enabled.
    level = logging.DEBUG if verbose else logging.INFO

    # Create the main logging formatter.
    # The timestamp is useful for tracking the progress of long-running
    # analyses and identifying where a problem occurred.
    formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%H:%M:%S"
    )

    # Create the console handler.
    # Messages are written to stderr so that standard output can remain
    # available for program results if needed.
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)

    # Configure the root logger.
    # force=True removes handlers that may have been installed previously
    # (for example by a test framework or another library).
    logging.basicConfig(
        level=level,
        handlers=[console_handler],
        force=True
    )

    # If a log file was specified, add a file handler.
    # The file receives the same messages as the console.
    if log_file is not None:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)

        logging.getLogger().addHandler(file_handler)
