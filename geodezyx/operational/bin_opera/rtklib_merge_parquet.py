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
  # scan a whole directory
  rtklib_merge_parquet -d /path/to/results

  # scan a directory with an experiment prefix
  rtklib_merge_parquet -d /path/to/results -e myexp

  # fast merge: append specific files to an existing _all.parquet
  rtklib_merge_parquet -d /path/to/results -e myexp --fast_merge \\
      -o /path/to/results/2024/001/run1.out /path/to/results/2024/002/run2.out

  # explicit list of parquet files (no directory scan)
  rtklib_merge_parquet -d /path/to/a.parquet /path/to/b.parquet -e myexp
""",
    )

    parser.add_argument(
        "-d",
        "--parquet_inp",
        nargs="+",
        required=True,
        metavar="PATH",
        help=(
            "Either a single directory path (all *.parquet files inside are collected "
            "recursively) or an explicit list of parquet file paths to merge."
        ),
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
            "Only used when --fast_merge is set and parquet_inp is a directory."
        ),
    )

    parser.add_argument(
        "-fm",
        "--fast_merge",
        action="store_true",
        help=(
            "If set, only merge the parquet files corresponding to --rtklib_out_files "
            "(or those in the explicit list) and append them to an already-existing "
            "_all.parquet file. "
            "If not set, scan the whole directory recursively for parquet files."
        ),
    )

    return parser.parse_args()


def rtklib_merge_prq_main():
    """Main entry point for the CLI."""
    args = parse_args()

    # Resolve parquet_inp: single directory or explicit list of files
    inp = args.parquet_inp
    if len(inp) == 1 and Path(inp[0]).is_dir():
        parquet_inp = inp[0]
        inp_label = f"directory: {parquet_inp}"
    else:
        parquet_inp = inp
        for f in parquet_inp:
            if not Path(f).exists():
                log.error(f"File not found: {f}")
                return 1
        inp_label = f"{len(parquet_inp)} explicit parquet file(s)"

    if args.fast_merge and args.rtklib_out_files is None and isinstance(parquet_inp, list):
        log.warning(
            "--fast_merge with an explicit file list: "
            "--rtklib_out_files is ignored; the provided files are merged directly."
        )

    if args.fast_merge and args.rtklib_out_files is None and isinstance(parquet_inp, str):
        log.error(
            "--fast_merge with a directory requires --rtklib_out_files to be provided."
        )
        return 1

    log.info(f"Parquet input:        {inp_label}")
    log.info(f"Experiment prefix:    '{args.exp_prefix}'")
    log.info(f"Fast merge:           {args.fast_merge}")
    if args.rtklib_out_files:
        log.info(f"RTKLIB out files:     {len(args.rtklib_out_files)} file(s)")

    try:
        all_prq_path = rtklib_merge_parquet(
            parquet_inp,
            rtklib_out_files=args.rtklib_out_files,
            fast_merge=args.fast_merge,
            exp_prefix=args.exp_prefix,
        )
        log.info(f"Merged parquet saved to: {all_prq_path}")
        return 0

    except Exception as e:
        log.error(f"Error during merge: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(rtklib_merge_prq_main())

