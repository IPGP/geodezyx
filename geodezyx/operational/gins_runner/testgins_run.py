#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 05/01/2026 19:12:11

@author: psakic
"""

from geodezyx.operational.gins_runner import get_rnx_eost, singugins_run
from geodezyx import utils


def run_testgins(dir_res_testgins=None, dir_rnxs_testgins=None):
    """
    Executes the test GINS workflow by downloading RINEX files and running the GINS process.

    Parameters:
    dir_res_testgins (str, optional): Path to the directory where the results will be stored.
        If not provided, a timestamped folder under '/root/030_RESULTS/' will be created.
    dir_rnxs_testgins (str, optional): Path to the directory where input RINEX files will be
        downloaded or read. Defaults to '/root/020_BDRNX/rnxs_testGINS_from_EOST' if not provided.

    Returns:
    None
    """
    # If no directory for RINEX files is provided, use the default path
    if not dir_rnxs_testgins:
        dir_rnxs_testgins = "/root/020_BDRNX/rnxs_testGINS_from_EOST"

    # Download or read RINEX files for the specified years
    get_rnx_eost(dir_rnxs_testgins, year_start=2022, year_end=2023)

    # If no results directory is provided, create a timestamped folder
    if not dir_res_testgins:
        timstp = utils.get_timestamp()
        dir_res_testgins = "/root/030_RESULTS/testGINS_" + timstp

    # Run the GINS process with the specified directories
    singugins_run(
        results_folder=dir_res_testgins,
        bdrnx_folder=dir_rnxs_testgins,
    )

    # The function does not return any value
    return None


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description=(__doc__ or "").strip())

    parser.add_argument(
        "-o",
        "--dir_res_testgins",
        type=str,
        default=None,
        help=(
            "Results directory. If omitted, a timestamped folder under "
            "'/root/030_RESULTS/' is created. (optional)"
        ),
    )
    parser.add_argument(
        "-r",
        "--dir_rnxs_testgins",
        type=str,
        default=None,
        help=(
            "Directory where input RINEX files will be downloaded/read. "
            "If omitted, defaults to '/root/020_BDRNX/rnxs_testGINS_from_EOST'. (optional)"
        ),
    )

    args = parser.parse_args(argv)
    run_testgins(
        dir_res_testgins=args.dir_res_testgins,
        dir_rnxs_testgins=args.dir_rnxs_testgins,
    )


if __name__ == "__main__":
    main()
