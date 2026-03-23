#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2026/03/23

@author: psakic

CLI interface for merging individual RTKLIB parquet files into a single
consolidated parquet file.
"""

import sys
import argparse
import logging
from pathlib import Path
from geodezyx.operational.soft_frontend.rtklib_frontend import rtklib_merge_parquet

log = logging.getLogger("geodezyx")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge individual RTKLIB parquet files into a single consolidated parquet file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  rtklib_merge_parquet -d /path/to/results
  rtklib_merge_parquet -d /path/to/results --exp_prefix myexp
  rtklib_merge_parquet -d /path/to/results -e myexp --fast_parquet_merge -o /path/to/results/2024/001/run1.out /path/to/results/2024/002/run2.out
""",
    )

    parser.add_argument(
        "-d",
        "--parquet_dir",
        required=True,
        help="Root directory where parquet files are located and where the merged file is saved.",
    )

    parser.add_argument(
        "-e",
        "--exp_prefix",
        default="",
        help="Prefix used to name the merged output file (<exp_prefix>_all.parquet). Default: ''",
    )

    parser.add_argument(
        "-o",
        "--rtklib_out_files",
        nargs="+",
        default=None,
        metavar="FILE",
        help=(
            "List of .out file paths produced by a previous RTKLIB run. "
            "Only used when --fast_parquet_merge is set."
        ),
    )

    parser.add_argument(
        "-fpm",
        "--fast_parquet_merge",
        action="store_true",
        help=(
            "If set, only merge the parquet files corresponding to --rtklib_out_files "
            "and append them to an already-existing _all.parquet file. "
            "If not set, scan the whole parquet_dir recursively for parquet files."
        ),
    )

    return parser.parse_args()


def rtklib_merge_prq_main():
    """Main entry point for the CLI."""
    args = parse_args()

    parquet_dir = args.parquet_dir
    parquet_dir_path = Path(parquet_dir)

    if not parquet_dir_path.exists():
        log.error(f"Error: Directory not found: {parquet_dir}")
        return 1

    if not parquet_dir_path.is_dir():
        log.error(f"Error: Path is not a directory: {parquet_dir}")
        return 1

    if args.fast_parquet_merge and not args.rtklib_out_files:
        log.error(
            "Error: --fast_parquet_merge requires --rtklib_out_files to be provided."
        )
        return 1

    log.info(f"Parquet directory:    {parquet_dir}")
    log.info(f"Experiment prefix:    '{args.exp_prefix}'")
    log.info(f"Fast parquet merge:   {args.fast_parquet_merge}")
    if args.rtklib_out_files:
        log.info(f"RTKLIB out files:     {len(args.rtklib_out_files)} file(s)")

    try:
        all_prq_path = rtklib_merge_parquet(
            parquet_dir=parquet_dir,
            exp_prefix=args.exp_prefix,
            rtklib_out_files=args.rtklib_out_files,
            fast_parquet_merge=args.fast_parquet_merge,
        )
        log.info(f"Merged parquet saved to: {all_prq_path}")
        return 0

    except Exception as e:
        log.error(f"Error during merge: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(rtklib_merge_prq_main())

