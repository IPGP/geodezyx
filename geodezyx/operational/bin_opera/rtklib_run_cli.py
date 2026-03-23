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
# then, no need of default values in the argparse arguments,
# as they will be overridden by YAML or CLI if provided,
# otherwise will use the function defaults


def parse_site_pair(sites_rovers_raw, sites_bases_raw=None):
    """
    Parse sites_rovers that may contain comma-separated 'ROVER,BASE' pairs.

    If any item in sites_rovers_raw contains a comma, **all** items are interpreted
    as explicit rover/base pairs (2-tuples), and sites_bases is forced to None so
    that make_pairs uses the explicit-pairs branch.

    Parameters
    ----------
    sites_rovers_raw : list of str
        Raw rover list, e.g. ``['TLSE', 'ZIMM']`` or ``['TLSE,GRAS', 'ZIMM,GRAS']``.
    sites_bases_raw : list of str or None
        Raw base list; ignored (set to None) when pair syntax is detected.

    Returns
    -------
    sites_rovers : list of str or list of (str, str)
    sites_bases  : list of str or None
    """
    if not sites_rovers_raw:
        return sites_rovers_raw, sites_bases_raw

    has_pairs = any("," in str(s) for s in sites_rovers_raw)

    if has_pairs:
        pairs = []
        for s in sites_rovers_raw:
            parts = str(s).split(",")
            if len(parts) != 2:
                raise ValueError(
                    f"Invalid rover/base pair: '{s}'. "
                    f"Expected 'ROVER,BASE' format (exactly one comma)."
                )
            pairs.append((parts[0].strip(), parts[1].strip()))
        return pairs, None

    return sites_rovers_raw, sites_bases_raw


def parse_args():
    """Parse CLI arguments and return only explicitly provided values."""
    parser = argparse.ArgumentParser(
        description="Run RTKLIB (rnx2rtkp) processing for multiple rover/base RINEX pairs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Standard mode: all rovers paired with all bases (cartesian product)
  rtklib_run --rnx_dir /data/rinex --cfgfile_generik config.conf --sites_rovers TLSE ZIMM --sites_bases GRAS --out_dir /results

  # Explicit pair mode: comma-separated ROVER,BASE pairs (sites_bases is ignored)
  rtklib_run --rnx_dir /data/rinex --cfgfile_generik config.conf --sites_rovers TLSE,GRAS ZIMM,BRUS --out_dir /results

  # YAML config (all CLI args override YAML settings)
  rtklib_run -y config.yaml
  rtklib_run -y config.yaml --sites_rovers TLSE,GRAS ZIMM,BRUS
""",
    )

    # YAML config file
    parser.add_argument(
        "-y",
        "--config_yaml",
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
        help=(
            "List of rover station names (4-char or 9-char codes, e.g., TLSE ZIMM BRUS). "
            "Alternatively, provide explicit ROVER,BASE pairs separated by a comma "
            "(e.g., TLSE,GRAS ZIMM,BRUS); in that case --sites_bases is ignored. "
            "Both standard names and ROVER,BASE syntax are also accepted in the YAML config."
        ),
    )
    parser.add_argument(
        "-b",
        "--sites_bases",
        nargs="+",
        help=(
            "Base station name(s) (4-char or 9-char code, e.g., GRAS or GRAS00FRA). "
            "Not needed when --sites_rovers uses ROVER,BASE pair syntax."
        ),
    )
    parser.add_argument("-o", "--out_dir", help="Output directory for results")

    # Optional arguments (no defaults - only explicit values)
    parser.add_argument(
        "-td",
        "--tmp_dir",
        help="Temporary directory (default: None, will use out_dir/TMP)",
    )
    parser.add_argument(
        "-pd",
        "--prod_dir",
        help="GNSS products directory (default: None, will use tmp_dir)",
    )
    parser.add_argument(
        "-ip",
        "--igs_prods",
        help=f"IGS product ID (default: {RTKLIB_RUN_DEFAULTS.get('igs_prods')})",
    )
    parser.add_argument(
        "-exp",
        "--exp_prefix",
        help="Output filename prefix (default: empty string)",
    )
    parser.add_argument(
        "-s",
        "--date_srt",
        help="Start date (YYYY-MM-DD or YYYY-DDD, default: None = all files)",
    )
    parser.add_argument(
        "-e",
        "--date_end",
        help="End date (YYYY-MM-DD or YYYY-DDD, default: None = all files)",
    )
    parser.add_argument(
        "-m",
        "--posmode",
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
        choices=["llh", "xyz", "enu", "nmea"],
        help=f"Solution format (default: {RTKLIB_RUN_DEFAULTS.get('solformat')})",
    )
    parser.add_argument(
        "-eph",
        "--sateph",
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
        help=f"Parallel workers (default: {RTKLIB_RUN_DEFAULTS.get('procs')})",
    )
    parser.add_argument(
        "-x",
        "--exe_path",
        help=f"Path to rnx2rtkp executable (default: {RTKLIB_RUN_DEFAULTS.get('exe_path')})",
    )

    args_namespace = parser.parse_args()
    args_dict = vars(args_namespace)

    # Filter: keep only arguments explicitly provided (not None, not False for action="store_true")
    return {k: v for k, v in args_dict.items() if v is not None and v is not False}


def rtklib_run_main():
    """Main entry point for the RTKLIB CLI."""
    kwargs_cli = parse_args()

    # Load YAML config if provided
    kwargs_cfg = {}
    if "config_yaml" in kwargs_cli:
        yaml_path = Path(kwargs_cli["config_yaml"])
        if not yaml_path.exists():
            log.error(f"Error: YAML config file not found: {kwargs_cli['config_yaml']}")
            sys.exit(1)
        with open(yaml_path, "r") as f:
            kwargs_cfg = yaml.safe_load(f) or {}

    # Merge: function defaults → YAML → CLI (each layer overrides the previous)
    kwargs_out = dict(RTKLIB_RUN_DEFAULTS)
    kwargs_out.update(kwargs_cfg)
    kwargs_out.update({k: v for k, v in kwargs_cli.items() if k != "config_yaml"})

    # Parse rover/base pairs if provided as 'ROVER,BASE' comma-separated strings
    # (works for both CLI and YAML sources)
    if "sites_rovers" in kwargs_out:
        kwargs_out["sites_rovers"], kwargs_out["sites_bases"] = parse_site_pair(
            kwargs_out.get("sites_rovers"),
            kwargs_out.get("sites_bases"),
        )

    # Parse dates if provided
    if kwargs_out.get("date_srt") and kwargs_out.get("date_end"):
        try:
            kwargs_out["date_srt"], kwargs_out["date_end"] = minmax_pattern_dt(
                kwargs_out["date_srt"], kwargs_out["date_end"]
            )
        except Exception as e:
            print(f"Error parsing dates: {e}", file=sys.stderr)
            sys.exit(1)

    log.info(f"RTKLIB RUN parameters:")
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
