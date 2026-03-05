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
from geodezyx import utils

import logging
log = logging.getLogger("geodezyx")

# Extract defaults from rtklib_run function at module level for synchronization
RTKLIB_RUN_DEFAULTS = utils.fct_def_args(rtklib_frontend.rtklib_run)


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
        "-d", "--rnx_dir", help="Directory containing RINEX files to process"
    )
    parser.add_argument(
        "-c",
        "--cfgfile_generik",
        help="Path to generic RTKLIB configuration file (.conf)",
    )
    parser.add_argument(
        "-r",
        "--sites_rovers",
        nargs="+",
        help="List of rover station names (4-char codes, e.g., TLSE ZIMM BRUS)",
    )
    parser.add_argument(
        "-b",
        "--sites_bases",
        nargs="+",
        help="Base station name(s) (4-char or 9-char code, e.g., GRAS or GRAS00FRA)",
    )
    parser.add_argument("-o", "--out_dir", help="Output directory for results")

    # Optional arguments
    parser.add_argument(
        "-td",
        "--tmp_dir",
        default=None,
        help="Temporary directory (default: None, will use out_dir/TMP)",
    )
    parser.add_argument(
        "-pd",
        "--prod_dir",
        default=None,
        help="GNSS products directory (default: None, will use tmp_dir)",
    )
    parser.add_argument(
        "-ip",
        "--igs_prods",
        default=RTKLIB_RUN_DEFAULTS.get("igs_prods"),
        help=f"IGS product ID (default: {RTKLIB_RUN_DEFAULTS.get('igs_prods')})",
    )
    parser.add_argument(
        "-exp",
        "--exp_prefix",
        default="",
        help="Output filename prefix (default: empty string)",
    )
    parser.add_argument(
        "-s",
        "--date_srt",
        default=None,
        help="Start date (YYYY-MM-DD or YYYY-DDD, default: None = all files)",
    )
    parser.add_argument(
        "-e",
        "--date_end",
        default=None,
        help="End date (YYYY-MM-DD or YYYY-DDD, default: None = all files)",
    )
    parser.add_argument(
        "-m",
        "--posmode",
        default=RTKLIB_RUN_DEFAULTS.get("posmode"),
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
        help=f"Position mode (default: {RTKLIB_RUN_DEFAULTS.get('posmode')})",
    )
    parser.add_argument(
        "-sol",
        "--solformat",
        default=RTKLIB_RUN_DEFAULTS.get("solformat"),
        choices=["llh", "xyz", "enu", "nmea"],
        help=f"Solution format (default: {RTKLIB_RUN_DEFAULTS.get('solformat')})",
    )
    parser.add_argument(
        "-eph",
        "--sateph",
        default=RTKLIB_RUN_DEFAULTS.get("sateph"),
        choices=["brdc", "precise"],
        help=f"Satellite ephemeris (default: {RTKLIB_RUN_DEFAULTS.get('sateph')})",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help=f"Force reprocessing (default: {RTKLIB_RUN_DEFAULTS.get('force')})",
    )
    parser.add_argument(
        "-ktmp",
        "--keep_tmp",
        action="store_true",
        help=f"Clean temp directory (default: {RTKLIB_RUN_DEFAULTS.get('clean_tmp')})",
    )
    parser.add_argument(
        "-p",
        "--procs",
        type=int,
        default=RTKLIB_RUN_DEFAULTS.get("procs"),
        help=f"Parallel workers (default: {RTKLIB_RUN_DEFAULTS.get('procs')})",
    )
    parser.add_argument(
        "-x",
        "--exe_path",
        default=RTKLIB_RUN_DEFAULTS.get("exe_path"),
        help=f"Path to rnx2rtkp executable (default: {RTKLIB_RUN_DEFAULTS.get('exe_path')})",
    )

    # Convert Namespace to dictionary
    args_namespace = parser.parse_args()
    return vars(args_namespace)


def rtklib_run_main():
    """Main entry point for the RTKLIB CLI."""
    kwargs_cli = parse_args()

    # Load YAML config if provided
    kwargs_cfg = {}
    if kwargs_cli.get("config_yaml"):
        yaml_path = Path(kwargs_cli["config_yaml"])
        if not yaml_path.exists():
            log.error(f"Error: YAML config file not found: {kwargs_cli['config_yaml']}")
            sys.exit(1)
        with open(yaml_path, "r") as f:
            kwargs_cfg = yaml.safe_load(f) or {}

    # Start with function defaults
    kwargs_out = dict(RTKLIB_RUN_DEFAULTS)

    # Override with YAML config
    if kwargs_cfg:
        kwargs_out.update(kwargs_cfg)

    # Override with CLI args (only if explicitly provided - not None)
    for arg_name, arg_value in kwargs_cli.items():
        if arg_name == "config_yaml":  # Skip the YAML config file argument itself
            continue
        if arg_value is not None:
            kwargs_out[arg_name] = arg_value

    # Parse dates if provided
    if kwargs_out.get("date_srt") and kwargs_out.get("date_end"):
        try:
            kwargs_out["date_srt"], kwargs_out["date_end"] = minmax_pattern_dt(
                kwargs_out["date_srt"], kwargs_out["date_end"]
            )
        except Exception as e:
            print(f"Error parsing dates: {e}", file=sys.stderr)
            sys.exit(1)

    log.info(f"RTKLIB RUN CLI parameters:")
    for key, value in sorted(kwargs_out.items()):
        log.info(f"  {key}: {value}")

    # Call the function
    try:
        results = rtklib_frontend.rtklib_run(**kwargs_out)
        successful = sum(1 for r in results if r)
        return 0 if successful > 0 else 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(rtklib_run_main())
