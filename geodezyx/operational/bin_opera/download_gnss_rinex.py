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
        description="CLI frontend for download_gnss_rinex"
    )
    parser.add_argument("-l", "--stations_list", nargs="+", required=True, help="List of station names")
    parser.add_argument("-c", "--datacenter", required=True, help="Archive center (e.g. igs_cddis, euref, etc.)")
    parser.add_argument("-o","--output_dir", required=True, help="Root directory to store RINEX files")
    parser.add_argument("-s","--startdate", required=True, help="Start date (YYYY-MM-DD)")
    parser.add_argument("-e","--enddate", required=True, help="End date (YYYY-MM-DD)")
    parser.add_argument("-a","--archtype", default="stat", help="Archive directory structure type")
    parser.add_argument("-u","--user", default="anonymous", help="FTP username")
    parser.add_argument("-p","--passwd", default="anonymous@isp.com", help="FTP password")
    parser.add_argument("-nr2","--no_rnx2", action="store_true", help="Download RINEX2 files")
    parser.add_argument("-nr3","--no_rnx3", action="store_true", help="Download RINEX3 files")
    parser.add_argument("-f","--force", action="store_true", help="Force download even if file exists")
    parser.add_argument("-q","--quiet", action="store_true", help="List available RINEXs without downloading")
    return parser.parse_args()

def main():
    args = parse_args()
    try:
        startdate = conv.date_pattern_2_dt(args.startdate)
        enddate = conv.date_pattern_2_dt(args.enddate)
        startdate = min((startdate, enddate))
        enddate = min((startdate, enddate))
    except Exception as e:
        print(f"Error parsing dates: {e}")
        sys.exit(1)

    statdico = {args.datacenter: args.stations_list}

    download_gnss_rinex(
        statdico=statdico,
        output_dir=args.output_dir,
        startdate=startdate,
        enddate=enddate,
        archtype=args.archtype,
        user=args.user,
        passwd=args.passwd,
        no_rnx2=args.no_rnx2,
        no_rnx3=args.no_rnx3,
        force=args.force,
        quiet_mode=args.quiet,
    )

if __name__ == "__main__":
    main()
