#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2026/03/04

@author: psakic

CLI interface for converting RTKLIB output files to Parquet format.
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from geodezyx import utils, files_rw

log = logging.getLogger("geodezyx")

def rtklib_parquet(resdir, pattern="*out", force=False):
    """
    Convert RTKLIB output files to Parquet format.

    Parameters
    ----------
    resdir : str
        Results directory containing RTKLIB output files.
    pattern : str, optional
        File pattern to search for (default: "*out").
    force : bool, optional
        Force conversion even if parquet file already exists (default: False).

    Returns
    -------
    list
        List of created/updated parquet files.
    """
    l_out = utils.find_recursive(resdir, pattern)
    f_prq_lis = []
    for f in l_out:
        f_prq = f.replace(".out", ".parquet")
        if not os.path.isfile(f_prq) or force:
            df_out2prq = files_rw.read_rtklib(f, return_df=True)
            df_out2prq.to_parquet(f_prq, engine="auto")
            f_prq_lis.append(f_prq)
            log.info(f"Created parquet file: {f_prq}")

    return f_prq_lis

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert RTKLIB output files to Parquet format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  rtklib_parquet /path/to/results
  rtklib_parquet -r /path/to/results --pattern "*.out"
  rtklib_parquet -r /path/to/results --force
  rtklib_parquet -r /path/to/results -p "*.out" -f
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

    try:
        f_prq_lis = rtklib_parquet(resdir, pattern=args.pattern, force=args.force)

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





