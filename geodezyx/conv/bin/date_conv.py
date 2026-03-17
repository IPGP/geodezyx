#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
geodezyx.bin.doy
================
CLI tool: convert a date string (any format accepted by date_pattern2dt)
and print all equivalent date representations.

Entry point registered as ``date_conv`` in pyproject.toml.
"""

import argparse

from geodezyx.conv.conv_time import (
    date_pattern2dt,
    dt2gpstime,
    dt2jjul_cnes,
    dt2mjd,
    dt2posix,
    dt2year_decimal,
)


def main():
    """
    Command-line interface for date_pattern2dt.

    Convert a date string in various formats to a Python datetime and print
    all equivalent representations.

    Accepted input formats
    ----------------------
    YYYY-MM-DD    : ISO date              e.g.  2026-03-17
    YYYY-DDD      : Year + DayOfYear      e.g.  2026-076
    WWWW-D        : GPS Week + DoW        e.g.  2408-2
    YYYYMMDD      : compact ISO           e.g.  20260317
    YYMMDD        : 2-digit year ISO      e.g.  260317
    JJJJJ         : CNES Julian day       e.g.  29665
    <free text>   : parsed by dateparser  e.g.  "17 March 2026"
    """
    parser = argparse.ArgumentParser(
        description=(
            "Convert a date string to a Python datetime and display all "
            "equivalent date representations."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Accepted input formats:
  YYYY-MM-DD   ISO date               e.g.  2026-03-17
  YYYY-DDD     Year + Day-of-Year     e.g.  2026-076
  WWWW-D       GPS Week + Day-of-Week e.g.  2408-2
  YYYYMMDD     Compact ISO            e.g.  20260317
  YYMMDD       2-digit year ISO       e.g.  260317
  JJJJJ        CNES Julian day        e.g.  29665
  <free text>  Parsed by dateparser   e.g.  "17 March 2026"
""",
    )
    parser.add_argument(
        "date_str",
        type=str,
        help="Input date string to convert (see accepted formats above).",
    )
    args = parser.parse_args()

    # --- parse the input ---
    try:
        date = date_pattern2dt(args.date_str)
    except Exception as exc:
        parser.error(f"Could not parse '{args.date_str}': {exc}")
        return  # unreachable, but satisfies linters

    # --- gather all representations ---
    gpsweek, gpsdow = dt2gpstime(date, secinweek=False, outputtype=int)
    _, gpssow       = dt2gpstime(date, secinweek=True,  outputtype=int)
    doy             = int(date.strftime("%j"))
    mjd             = dt2mjd(date)
    jjcnes          = int(dt2jjul_cnes(date, onlydays=True))
    posix_ts        = dt2posix(date)
    year_dec        = dt2year_decimal(date)

    # --- formatted output ---
    sep = "-" * 55
    print(sep)
    print(f"  Input string  : {args.date_str}")
    print(sep)
    print(f"  {'YYYY-MM-DD':<28}: {date.strftime('%Y-%m-%d')}")
    print(f"  {'YYYYMMDD':<28}: {date.strftime('%Y%m%d')}")
    print(f"  {'YYMMDD':<28}: {date.strftime('%y%m%d')}")
    print(f"  {'YYYY-DDD (Year-DayOfYear)':<28}: {date.strftime('%Y')}-{doy:03d}")
    print(f"  {'WWWW-D (GPS Week-DayOfWeek)':<28}: {gpsweek:04d}-{gpsdow}")
    print(f"  {'GPS Week':<28}: {gpsweek}")
    print(f"  {'GPS Day of Week':<28}: {gpsdow}")
    print(f"  {'GPS Second of Week':<28}: {gpssow}")
    print(f"  {'MJD':<28}: {mjd:.6f}")
    print(f"  {'CNES Julian Day (JJJJJ)':<28}: {jjcnes}")
    print(f"  {'POSIX timestamp':<28}: {posix_ts:.0f}")
    print(f"  {'Decimal Year':<28}: {year_dec:.8f}")
    print(f"  {'ISO 8601 datetime':<28}: {date.strftime('%Y-%m-%dT%H:%M:%S')}")
    print(sep)


if __name__ == "__main__":
    main()

