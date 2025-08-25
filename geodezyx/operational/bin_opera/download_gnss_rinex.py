#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06/06/2025 11:23:06

@author: psakic
"""

import argparse
import datetime as dt
import sys
from geodezyx import conv

from geodezyx.operational.download_rinex import download_gnss_rinex


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download GNSS RINEX files from various data centers. "
                   "Supports multiple archive centers (IGS, EUREF, etc.) and handles "
                   "both RINEX2 and RINEX3 formats with parallel downloading capabilities."
    )

    parser.add_argument(
        "-l", "--stations_list",
        nargs="+",
        required=True,
        help="List of 4 or 9-character GNSS station names (e.g., ZIMM TLSE00FRA BRUS)"
    )
    parser.add_argument(
        "-c", "--datacenter",
        required=True,
        help="Data center/archive to download from. Supported values: "
             "'igs_cddis' or 'igs' (CDDIS data center), "
             "'igs_sopac' (SOPAC/UCSD/SIO data center), "
             "'igs_ign' (IGN's main server at St Mandé), "
             "'igs_ign_ensg' (IGN's secondary server at ENSG), "
             "'igs_bkg' (BKG's IGS data center), "
             "'rgp' (IGN's RGP main server at St Mandé), "
             "'rgp_ensg' (IGN's RGP, secondary server at ENSG, Marne-la-Vallée), "
             "'sonel' (SONEL data center), "
             "'euref' (EPN data center hosted at ROB), "
             "'nav' or 'brdc' (navigation files from ROB server), "
             "'nav_rt' or 'brdc_rt' (real-time navigation files from BKG server)"
    )
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="Root directory where RINEX files will be stored in organized subdirectories"
    )
    parser.add_argument(
        "-s", "--startdate",
        required=True,
        help="Start date for download period (YYYY-MM-DD or YYYY-DDD format)"
    )
    parser.add_argument(
        "-e", "--enddate",
        required=True,
        help="End date for download period (YYYY-MM-DD or YYYY-DDD format)"
    )
    parser.add_argument(
        "-a", "--archtype",
        default="stat",
        help="Archive directory structure type. Options: 'stat' (station-based), "
             "'year' (year-based), 'daily' (daily structure). "
             "Exemple: 'stat/year/doy'. "
             "Default: 'year/doy'"
    )
    parser.add_argument(
        "-u", "--user",
        default="anonymous",
        help="FTP username for authentication. Use 'anonymous' for public access"
    )
    parser.add_argument(
        "-p", "--passwd",
        default="anonymous@isp.com",
        help="FTP password for authentication. Use email for anonymous access"
    )
    parser.add_argument(
        "-nr2", "--no_rnx2",
        action="store_true",
        help="Skip downloading RINEX2 format files (short filename convention)"
    )
    parser.add_argument(
        "-nr3", "--no_rnx3",
        action="store_true",
        help="Skip downloading RINEX3 format files (long filename convention)"
    )
    parser.add_argument(
        "-f", "--force",
        action="store_true",
        help="Force re-download files even if they already exist locally"
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Quiet mode: list available files on server without downloading them"
    )
    parser.add_argument(
        "-pd", "--parallel_download",
        type=int,
        default=1,
        help="Number of parallel download threads to use. Higher values speed up "
             "downloads but increase server load. Recommended: 1-8. Default: 1"
    )

    return parser.parse_args()

def main():
    args = parse_args()
    try:
        srtdate, enddate = conv.minmax_pattern_dt(args.startdate, args.enddate)
    except Exception as e:
        print(f"Error parsing dates: {e}")
        sys.exit(1)

    statdico = {args.datacenter: args.stations_list}

    download_gnss_rinex(
        statdico=statdico,
        output_dir=args.output_dir,
        startdate=srtdate,
        enddate=enddate,
        archtype=args.archtype,
        user=args.user,
        passwd=args.passwd,
        no_rnx2=args.no_rnx2,
        no_rnx3=args.no_rnx3,
        force=args.force,
        quiet_mode=args.quiet,
        parallel_download=args.parallel_download,
    )


if __name__ == "__main__":
    main()
