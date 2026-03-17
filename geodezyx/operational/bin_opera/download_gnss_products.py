#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI entry point for GNSS products download.

Mirrors the ``download_gnss_rinex`` CLI but targets orbital / clock products.

@author: psakic
"""

import argparse
import sys

from geodezyx import conv
from geodezyx.operational.download_gnss.download_prods import download_gnss_products


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download GNSS orbit / clock products from various IGS data centers. "
                    "Supports multiple analysis centers (wum, cod, jpl, …) and product types "
                    "(sp3, clk, erp, bia, …) with optional parallel downloading."
    )

    # ------------------------------------------------------------------ dates
    parser.add_argument(
        "-s", "--startdate",
        required=True,
        help="Start date of the download period (YYYY-MM-DD or YYYY-DDD format)."
    )
    parser.add_argument(
        "-e", "--enddate",
        default=None,
        help="End date of the download period (YYYY-MM-DD or YYYY-DDD format). "
             "If omitted, only startdate is downloaded."
    )

    # --------------------------------------------------------------- required
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="Root directory where products will be stored."
    )
    parser.add_argument(
        "-c", "--datacenter",
        required=True,
        help="IGS data / archive center to use. Supported values: "
             "'ign' (IGN St Mandé – default), "
             "'cddis' (NASA CDDIS, SFTP), "
             "'ign_ensg' (IGN ENSG secondary server), "
             "'bkg' (BKG data center), "
             "'wuhan' (WHU/IGG data center)."
    )

    # ---------------------------------------------------------- product config
    parser.add_argument(
        "-ac", "--ac_names",
        nargs="+",
        default=["grg", "cod"],
        metavar="AC",
        help="One or more analysis-center names (e.g., grg cod jpl igs). "
             "For the new long-name convention you may embed the latency tag "
             "directly, e.g. IGS0OPSRAP. Default: grg cod."
    )
    parser.add_argument(
        "-pt", "--prod_types",
        nargs="+",
        default=["sp3", "clk"],
        metavar="PROD",
        help="One or more product types to download (e.g., sp3 clk erp bia snx). "
             "Default: sp3 clk."
    )
    parser.add_argument(
        "-rp", "--remove_patterns",
        nargs="+",
        default=["ULA"],
        metavar="PAT",
        help="Filename patterns to exclude from downloads. Default: ULA."
    )

    # ---------------------------------------------------------- archive layout
    parser.add_argument(
        "-a", "--archtype",
        default="week",
        help="Local archive sub-directory structure. Examples: "
             "'week' (GPS week/DOW), 'year/doy', 'none'. "
             "See effective_save_dir_orbit for full details. Default: 'week'."
    )

    # ---------------------------------------------------------- naming / repro
    parser.add_argument(
        "--no_new_name_conv",
        action="store_true",
        help="Disable the new long-filename naming convention search "
             "(only use old short names)."
    )
    parser.add_argument(
        "--mgex",
        action="store_true",
        help="Target MGEX products instead of standard IGS products."
    )
    parser.add_argument(
        "-r", "--repro",
        type=int,
        default=0,
        help="Reprocessing campaign number (0 = operational). Default: 0."
    )
    parser.add_argument(
        "--dow_manu",
        default="auto",
        help="Day-of-week override for weekly products. "
             "'auto' (default) = computed from date; "
             "'none' = no DOW in regex (weekly search); "
             "integer 0-7 = specific DOW."
    )

    # --------------------------------------------------------------- download
    parser.add_argument(
        "-pd", "--parallel_download",
        type=int,
        default=1,
        help="Number of parallel download threads. "
             "Recommended: 1-8 (CDDIS is forced to 1). Default: 1."
    )
    parser.add_argument(
        "-f", "--force",
        action="store_true",
        help="Force re-download even if files already exist locally."
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Quiet mode: list available files without downloading them."
    )

    # ------------------------------------------------- crawl table save/load
    parser.add_argument(
        "--save_crawl",
        default=None,
        metavar="CSV_PATH",
        help="Save the FTP crawl table to this CSV file path."
    )
    parser.add_argument(
        "--load_crawl",
        default=None,
        metavar="CSV_PATH",
        help="Load a previously saved FTP crawl table instead of crawling live."
    )
    parser.add_argument(
        "--skip_crawl",
        action="store_true",
        help="Skip the FTP crawl step entirely (requires --load_crawl)."
    )
    parser.add_argument(
        "--save_all_files",
        default=None,
        metavar="CSV_PATH",
        help="Save ALL remote filenames found on the FTP server to this CSV path."
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # ------------------------------------------------------------------ dates
    try:
        startdate, enddate = conv.minmax_pattern_dt(args.startdate, args.enddate)
    except Exception as e:
        print(f"Error parsing dates: {e}", file=sys.stderr)
        sys.exit(1)

    # ----------------------------------------------------------------- dow_manu
    dow_manu_val = args.dow_manu
    if dow_manu_val == "auto":
        dow_manu_val = False
    elif dow_manu_val == "none":
        dow_manu_val = None
    else:
        try:
            dow_manu_val = int(dow_manu_val)
        except ValueError:
            print(
                f"Error: --dow_manu must be 'auto', 'none', or an integer, got '{args.dow_manu}'",
                file=sys.stderr,
            )
            sys.exit(1)

    download_gnss_products(
        archive_dir=args.output_dir,
        startdate=startdate,
        enddate=enddate,
        ac_names=tuple(args.ac_names),
        prod_types=tuple(args.prod_types),
        remove_patterns=tuple(args.remove_patterns),
        archtype=args.archtype,
        new_name_conv=not args.no_new_name_conv,
        parallel_download=args.parallel_download,
        archive_center=args.datacenter,
        mgex=args.mgex,
        repro=args.repro,
        dow_manu=dow_manu_val,
        quiet_mode=args.quiet,
        force=args.force,
        path_ftp_crawled_files_save=args.save_crawl,
        path_ftp_crawled_files_load=args.load_crawl,
        skip_crawl=args.skip_crawl,
        path_all_ftp_files_save=args.save_all_files,
    )


if __name__ == "__main__":
    main()

