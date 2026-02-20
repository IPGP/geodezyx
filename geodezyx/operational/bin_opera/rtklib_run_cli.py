#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2026/02/19

@author: psakic

CLI interface for running RTKLIB processing with the GeodeZYX toolbox.
"""

import argparse
import sys
import yaml
from pathlib import Path
from geodezyx.conv import minmax_pattern_dt
from geodezyx.operational.soft_frontend import rtklib_frontend

import logging
log = logging.getLogger('geodezyx')

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run RTKLIB (rnx2rtkp) processing for multiple rover/base RINEX pairs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  rtklib_run --rnx_dir /data/rinex --cfgfile_generik config.conf --sites_rovers TLSE ZIMM --site_base GRAS0 --out_dir /results
  rtklib_run -y config.yaml
  rtklib_run -y config.yaml --sites_rovers TLSE ZIMM  # Override YAML settings
""",
    )

    # YAML config file
    parser.add_argument(
        "-y",
        "--config_yaml",
        default=None,
        help="YAML configuration file with all parameters. CLI args override YAML settings.",
    )

    # Required arguments (unless provided in YAML)
    parser.add_argument(
        "-r", "--rnx_dir", help="Directory containing RINEX files to process"
    )
    parser.add_argument(
        "-c",
        "--cfgfile_generik",
        help="Path to generic RTKLIB configuration file (.conf)",
    )
    parser.add_argument(
        "-ro",
        "--sites_rovers",
        nargs="+",
        help="List of rover station names (4-char codes, e.g., TLSE ZIMM BRUS)",
    )
    parser.add_argument(
        "-b",
        "--site_base",
        help="Base station name (4-char or 9-char code, e.g., GRAS or GRAS00FRA)",
    )
    parser.add_argument("-o", "--out_dir", help="Output directory for results")

    # Optional arguments    parser.add_argument("-td","--tmp_dir", default=None, help="Temporary directory")
    parser.add_argument(
        "-pd", "--prod_dir", default=None, help="GNSS products directory"
    )
    parser.add_argument(
        "-ip",
        "--igs_prods",
        default="GRG0OPSULT",
        help="IGS product ID (default: GRG0OPSULT)",
    )
    parser.add_argument(
        "-exp", "--exp_prefix", default="", help="Output filename prefix"
    )
    parser.add_argument(
        "-s", "--date_srt", default=None, help="Start date (YYYY-MM-DD or YYYY-DDD)"
    )
    parser.add_argument(
        "-e", "--date_end", default=None, help="End date (YYYY-MM-DD or YYYY-DDD)"
    )
    parser.add_argument(
        "-m",
        "--posmode",
        default=None,
        choices=[
            "single",
            "dgps",
            "kinematic",
            "static",
            "static-start",
            "movingbase",
            "fixed",
            "ppp-kine",
            "ppp-static",
            "ppp-fixed",
        ],
        help="Position mode",
    )
    parser.add_argument(
        "-fmt",
        "--solformat",
        default=None,
        choices=["llh", "xyz", "enu", "nmea"],
        help="Solution format",
    )
    parser.add_argument(
        "-eph",
        "--sateph",
        default=None,
        choices=["brdc", "precise"],
        help="Satellite ephemeris",
    )
    parser.add_argument("-f", "--force", action="store_true", help="Force reprocessing")
    parser.add_argument(
        "-cln", "--clean_tmp", action="store_true", help="Clean temp directory"
    )
    parser.add_argument(
        "-p", "--procs", type=int, default=8, help="Parallel workers"
    )
    parser.add_argument(
        "-x", "--exe_path", default="rnx2rtkp", help="Path to rnx2rtkp executable"
    )

    return parser.parse_args()


def main():
    """Main entry point for the RTKLIB CLI."""
    args = parse_args()

    # Load YAML config if provided
    config = {}
    if args.config_yaml:
        yaml_path = Path(args.config_yaml)
        if not yaml_path.exists():
            print(
                f"Error: YAML config file not found: {args.config_yaml}",
                file=sys.stderr,
            )
            sys.exit(1)
        with open(yaml_path, "r") as f:
            config = yaml.safe_load(f) or {}

    # Build kwargs: merge YAML and CLI args (CLI overrides YAML)
    kwargs = {}
    for arg_name in vars(args):
        if arg_name == "config_yaml":  # Skip the YAML config file argument itself
            continue

        cli_val = getattr(args, arg_name)
        yaml_val = config.get(arg_name)

        # CLI takes precedence if explicitly set
        if cli_val is not None:
            kwargs[arg_name] = cli_val
        elif yaml_val is not None:
            kwargs[arg_name] = yaml_val

    # Parse dates if provided
    if "date_srt" in kwargs and "date_end" in kwargs:
        try:
            kwargs["date_srt"], kwargs["date_end"] = minmax_pattern_dt(
                kwargs["date_srt"], kwargs["date_end"]
            )
        except Exception as e:
            print(f"Error parsing dates: {e}", file=sys.stderr)
            sys.exit(1)

    log.info(f"RTKLIB RUN CLI parameters:")
    for key in kwargs:
        log.info(f"* {key}: {kwargs[key]}")

    # Call the function
    try:
        results = rtklib_frontend.rtklib_run(**kwargs)
        successful = sum(1 for r in results if r)
        return 0 if successful > 0 else 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
