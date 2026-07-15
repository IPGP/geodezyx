#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2026/06/03

@author: psakic

CLI interface for converting merged RTKLIB parquet files to CSV format,
processed by rover/base pairs.
"""

import sys
import argparse
import logging
from pathlib import Path
from geodezyx.operational.soft_frontend.rtklib_frontend import parquet2csv
from geodezyx import utils

log = logging.getLogger("geodezyx")

# Extract defaults from parquet2csv function at module level for synchronization
PARQUET2CSV_DEFAULTS = utils.fct_def_args(parquet2csv)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert merged RTKLIB parquet file to CSV format, processed by rover/base pairs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert parquet to CSV with default 15min resampling
  parquet2csv -i /path/to/merged_all.parquet -o /path/to/output_csv

  # Convert with 1min resampling
  parquet2csv -i /path/to/merged_all.parquet -o /path/to/output_csv -s 1min

  # Convert with 1 hour resampling
  parquet2csv -i /path/to/merged_all.parquet -o /path/to/output_csv -s 1H

  # Custom resampling (daily)
  parquet2csv -i /path/to/merged_all.parquet -o /path/to/output_csv -s 1D
""",
    )

    parser.add_argument(
        "-i",
        "--prq_inp",
        required=True,
        metavar="FILE",
        help="Path to the merged parquet file containing rover/base pair GNSS solutions",
    )

    parser.add_argument(
        "-o",
        "--out_dir",
        required=True,
        metavar="DIR",
        help="Output directory where CSV files will be saved",
    )

    parser.add_argument(
        "-s",
        "--sample",
        default=PARQUET2CSV_DEFAULTS.get("sample", "15min"),
        help=(
            "Resampling interval for position data "
            "(default: %(default)s). "
            "Examples: '1min', '15min', '1H' (1 hour), '1D' (1 day)"
        ),
    )

    return parser.parse_args()


def parquet2csv_main():
    """Main entry point for the parquet2csv CLI."""
    args = parse_args()

    # Validate input file
    prq_inp_path = Path(args.prq_inp)
    if not prq_inp_path.exists():
        log.error(f"Input parquet file not found: {args.prq_inp}")
        return 1

    if not prq_inp_path.is_file():
        log.error(f"Input path is not a file: {args.prq_inp}")
        return 1

    # Create output directory if it doesn't exist
    try:
        out_dir_path = Path(args.out_dir)
        out_dir_path.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        log.error(f"Failed to create output directory '{args.out_dir}': {e}")
        return 1

    log.info("parquet2csv conversion parameters:")
    log.info(f"  Input parquet file: {args.prq_inp}")
    log.info(f"  Output directory:   {args.out_dir}")
    log.info(f"  Resampling interval: {args.sample}")

    try:
        parquet2csv(
            prq_inp=args.prq_inp,
            out_dir=args.out_dir,
            sample=args.sample,
        )
        log.info("parquet2csv conversion completed successfully")
        return 0

    except Exception as e:
        log.error(f"Error during parquet2csv conversion: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(parquet2csv_main())

