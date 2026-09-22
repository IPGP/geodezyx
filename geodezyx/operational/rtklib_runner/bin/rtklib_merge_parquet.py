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
from geodezyx import utils, conv
from geodezyx.operational.rtklib_runner.rtklib_parquet import rtklib_merge_prq

log = logging.getLogger("geodezyx")

# Extract defaults from rtklib_merge_prq function at module level for synchronization
RTKLIB_MERGE_DEFAULTS = utils.fct_def_args(rtklib_merge_prq)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge individual RTKLIB parquet files into a single consolidated parquet file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # scan a whole directory
  rtklib_merge_prq -i /path/to/results

  # scan a directory with an experiment prefix
  rtklib_merge_prq -i /path/to/results -x myexp

  # specify output directory
  rtklib_merge_prq -i /path/to/results -o /path/to/output

  # fast merge: append specific files to an existing _all.parquet
   rtklib_merge_prq -i /path/to/results -x myexp --fast_merge \\
       -rof /path/to/results/2024/001/run1.out /path/to/results/2024/002/run2.out

   # explicit list of parquet files (no directory scan)
   rtklib_merge_prq -i /path/to/a.parquet /path/to/b.parquet -x myexp
 
   # with resampling (15min)
   rtklib_merge_prq -i /path/to/results -x myexp -p 15min
 
   # with resampling (1 hour)
   rtklib_merge_prq -i /path/to/results -x myexp -p 1H

   # with date filtering (start and end)
   rtklib_merge_prq -i /path/to/results -x myexp -s 2024-01-01 -e 2024-01-31

   # with date filtering (start + days)
   rtklib_merge_prq -i /path/to/results -x myexp -s 2024-01-01 -d 7

   # with date filtering (end + days)
   rtklib_merge_prq -i /path/to/results -x myexp -e 2024-01-31 -d 7

   # with output directory and resampling
   rtklib_merge_prq -i /path/to/results -x myexp -o /output/dir -p 15min
""",
    )

    parser.add_argument(
        "-i",
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
        "-x",
        "--exp_prefix",
        default="",
        help="Prefix used to name the merged output file (<exp_prefix>_all.parquet). Default: ''",
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

    parser.add_argument(
        "-rof",
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
        "-p",
        "--sample",
        default=RTKLIB_MERGE_DEFAULTS.get("sample"),
        help=(
            "Resampling interval for position data (default: None, no resampling). "
            "If provided, resamples each table to the specified interval before merging. "
            "Examples: '1min', '15min', '1H' (1 hour), '1D' (1 day)"
        ),
    )

    parser.add_argument(
        "-s",
        "--start_date",
        dest="start_date",
        type=conv.date_pattern2dt,
        default=None,
        help=(
            "Start date for filtering parquet files (flexible format: YYYY-MM-DD, YYYY-DDD, etc.). "
            "Assumes directory structure: <path>/<year>/<doy>/*.parquet"
        ),
    )

    parser.add_argument(
        "-e",
        "--end_date",
        type=conv.date_pattern2dt,
        default=None,
        help=(
            "End date for filtering parquet files (flexible format: YYYY-MM-DD, YYYY-DDD, etc.). "
            "Assumes directory structure: <path>/<year>/<doy>/*.parquet"
        ),
    )

    parser.add_argument(
        "-d",
        "--days",
        type=int,
        default=None,
        help=(
            "Number of days to process. Only used if start_date XOR end_date is provided. "
            "If start_date is given, processes N days starting from start_date. "
            "If end_date is given, processes N days ending at end_date."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        default=None,
        metavar="PATH",
        help=(
            "Output directory where the merged parquet file will be saved. "
            "If not provided, the output is saved to the input directory "
            "(or the directory of the first file if an explicit list is provided)."
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

    # Validate date parameters
    if (args.start_date is not None or args.end_date is not None) and args.days is not None:
        if args.start_date is not None and args.end_date is not None:
            log.error(
                "Cannot use both --start_date and --end_date together with --days. "
                "Use either (--start_date and optionally --days) or (--end_date and optionally --days)."
            )
            return 1

    log.info(f"Parquet input:        {inp_label}")
    log.info(f"Experiment prefix:    '{args.exp_prefix}'")
    log.info(f"Fast merge:           {args.fast_merge}")
    if args.rtklib_out_files:
        log.info(f"RTKLIB out files:     {len(args.rtklib_out_files)} file(s)")
    if args.output_dir:
        log.info(f"Output directory:     {args.output_dir}")
    if args.sample:
        log.info(f"Resampling interval:  {args.sample}")
    if args.start_date or args.end_date:
        log.info(f"Date range:           {args.start_date or 'N/A'} to {args.end_date or 'N/A'}")
    if args.days:
        log.info(f"Days:                 {args.days}")

    try:
        all_prq_path = rtklib_merge_prq(
            parquet_inp,
            exp_prefix=args.exp_prefix,
            fast_merge=args.fast_merge,
            rtklib_out_files=args.rtklib_out_files,
            sample=args.sample,
            start_date=args.start_date,
            end_date=args.end_date,
            days=args.days,
            output_dir=args.output_dir,
        )
        log.info(f"Merged parquet saved to: {all_prq_path}")
        return 0

    except Exception as e:
        log.error(f"Error during merge: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(rtklib_merge_prq_main())

