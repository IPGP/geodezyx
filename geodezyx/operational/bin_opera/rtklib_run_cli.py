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
from geodezyx.conv.conv_time_orig import minmax_pattern_dt
from geodezyx.operational.soft_frontend import rtklib_frontend



def parse_args():
    parser = argparse.ArgumentParser(
        description="Run RTKLIB (rnx2rtkp) processing for multiple rover/base RINEX pairs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  rtklib_run --rnx_dir /data/rinex --cfgfile_generik config.conf --sites_rovers TLSE ZIMM --site_base GRAS0 --out_dir /results
  rtklib_run -y config.yaml
  rtklib_run -y config.yaml --sites_rovers TLSE ZIMM  # Override YAML settings
"""
    )

    # YAML config file
    parser.add_argument("-y", "--config_yaml", default=None,
                       help="YAML configuration file with all parameters. CLI args override YAML settings.")

    # Required arguments (unless provided in YAML)
    parser.add_argument("--rnx_dir", help="Directory containing RINEX files to process")
    parser.add_argument("--cfgfile_generik", help="Path to generic RTKLIB configuration file (.conf)")
    parser.add_argument("--sites_rovers", nargs="+", help="List of rover station names (4-char codes, e.g., TLSE ZIMM BRUS)")
    parser.add_argument("--site_base", help="Base station name (4-char or 9-char code, e.g., GRAS or GRAS00FRA)")
    parser.add_argument("--out_dir", help="Output directory for results")

    # Optional arguments
    parser.add_argument("--tmp_dir", default=None, help="Temporary directory")
    parser.add_argument("--prod_dir", default=None, help="GNSS products directory")
    parser.add_argument("--igs_prods", default="GRG0OPSULT", help="IGS product ID (default: GRG0OPSULT)")
    parser.add_argument("--exp_prefix", default="", help="Output filename prefix")
    parser.add_argument("--date_srt", default=None, help="Start date (YYYY-MM-DD or YYYY-DDD)")
    parser.add_argument("--date_end", default=None, help="End date (YYYY-MM-DD or YYYY-DDD)")
    parser.add_argument("--posmode", default=None,
                       choices=["single", "dgps", "kinematic", "static", "moving-base", "fixed", "ppp-kinematic", "ppp-static"],
                       help="Position mode")
    parser.add_argument("--solformat", default=None, choices=["llh", "xyz", "enu", "nmea"], help="Solution format")
    parser.add_argument("--sateph", default=None,
                       choices=["brdc", "precise", "brdc+sbas", "brdc+ssrapc", "brdc+ssrcom"],
                       help="Satellite ephemeris")
    parser.add_argument("--force", action="store_true", help="Force reprocessing")
    parser.add_argument("--clean_tmp", action="store_true", help="Clean temp directory")
    parser.add_argument("--procs", type=int, default=4, help="Parallel workers (default: 4)")
    parser.add_argument("--exe_path",
                       default="/home/sakic/SOFTWARE/RTKLIB_explorer/RTKLIB/app/consapp/rnx2rtkp/gcc/rnx2rtkp",
                       help="Path to rnx2rtkp executable")

    return parser.parse_args()


def main():
    """Main entry point for the RTKLIB CLI."""
    args = parse_args()

    # Load YAML config if provided
    config = {}
    if args.config_yaml:
        yaml_path = Path(args.config_yaml)
        if not yaml_path.exists():
            print(f"Error: YAML config file not found: {args.config_yaml}", file=sys.stderr)
            sys.exit(1)
        with open(yaml_path, 'r') as f:
            config = yaml.safe_load(f) or {}

    # Build kwargs: merge YAML and CLI args (CLI overrides YAML)
    kwargs = {}
    for arg_name in vars(args):
        if arg_name == 'config_yaml':  # Skip the YAML config file argument itself
            continue

        cli_val = getattr(args, arg_name)
        yaml_val = config.get(arg_name)

        # CLI takes precedence if explicitly set
        if cli_val is not None:
            kwargs[arg_name] = cli_val
        elif yaml_val is not None:
            kwargs[arg_name] = yaml_val

    # Parse dates if provided
    if 'date_srt' in kwargs and 'date_end' in kwargs:
        try:
            kwargs['date_srt'], kwargs['date_end'] = minmax_pattern_dt(kwargs['date_srt'], kwargs['date_end'])
        except Exception as e:
            print(f"Error parsing dates: {e}", file=sys.stderr)
            sys.exit(1)

    # Call the function
    try:
        results = rtklib_frontend.rtklib_run(**kwargs)
        successful = sum(1 for r in results if r)
        print(f"Processed {len(results)} pairs: {successful} successful, {len(results) - successful} failed")
        return 0 if successful > 0 else 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())









