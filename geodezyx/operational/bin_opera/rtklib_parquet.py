#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2026/03/04

@author: psakic

CLI interface for converting RTKLIB output files to Parquet format.
"""

import sys
import argparse
import logging
from pathlib import Path
from geodezyx import utils, files_rw
from geodezyx.operational.soft_frontend.rtklib_frontend import rtklib_parquet

log = logging.getLogger("geodezyx")

# Extract defaults from rtklib_parquet function at module level for synchronization
RTKLIB_PARQUET_DEFAULTS = utils.fct_def_args(rtklib_parquet)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert RTKLIB output files to Parquet format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  rtklib_parquet -r /path/to/results
  rtklib_parquet -r /path/to/results --pattern "*.out"
  rtklib_parquet -r /path/to/results --force
  rtklib_parquet -r /path/to/results --sample 15min
  rtklib_parquet -r /path/to/results -p "*.out" -f -s 1H
""",
    )

    # Required arguments
    parser.add_argument(
        "-r",
        "--resdir",
        default=None,
        required=True,
        help="Results directory (alternative to positional argument)",
    )

    parser.add_argument(
        "-p",
        "--pattern",
        default="*out",
        help="File pattern to search for (default: *out)",
    )

    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force conversion even if parquet file already exists",
    )

    parser.add_argument(
        "-s",
        "--sample",
        default=RTKLIB_PARQUET_DEFAULTS.get("sample"),
        help=(
            "Resampling interval for position data (default: None, no resampling). "
            "Examples: '1min', '15min', '1H' (1 hour), '1D' (1 day)"
        ),
    )

    return parser.parse_args()


def rtklib_prq_main():
    """Main entry point for the CLI."""
    args = parse_args()
    # Configure logging
    # Determine results directory
    resdir = args.resdir

    if not resdir:
        log.error("Error: Results directory must be provided")
        return 1

    resdir_path = Path(resdir)
    if not resdir_path.exists():
        log.error(f"Error: Directory not found: {resdir}")
        return 1

    if not resdir_path.is_dir():
        log.error(f"Error: Path is not a directory: {resdir}")
        return 1

    log.info(f"Processing directory: {resdir}")
    log.info(f"Pattern: {args.pattern}")
    log.info(f"Force conversion: {args.force}")
    if args.sample:
        log.info(f"Resampling interval: {args.sample}")

    try:
        f_prq_lis = rtklib_parquet(resdir, pattern=args.pattern, force=args.force, sample=args.sample)

        if f_prq_lis:
            log.info(f"Successfully converted {len(f_prq_lis)} file(s) to Parquet format")
            return 0
        else:
            log.warning(f"No files matched pattern '{args.pattern}' or all parquet files already exist")
            return 0

    except Exception as e:
        log.error(f"Error during conversion: {e}", exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(rtklib_prq_main())





