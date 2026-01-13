#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 05/01/2026 19:12:11

@author: psakic
"""

from geodezyx.operational.gins_runner import get_rnx_eost
from geodezyx.operational.gins_runner.singugins_run import singugins_run
from geodezyx import utils

def run_testgins(results_folder=None, rnxs_folder=None, no_download_rnxs=False):
    """
    Executes the test GINS workflow by downloading RINEX files and running the GINS process.

    Parameters
    ----------
    results_folder : str, optional
        Path to the directory where the results will be stored.
        If None, a timestamped folder under '/root/030_RESULTS/' will be created.
        Default is None.
    rnxs_folder : str, optional
        Path to the directory where input RINEX files will be downloaded or read.
        If None, defaults to '/root/020_BDRNX/rnxs_testGINS_from_EOST'.
        Default is None.
    no_download_rnxs : bool, optional
        If True, the script will not download RINEX files from EOST server
        and will use existing files in the specified rnxs_folder.
        Default is False.

    Returns
    -------
    None

    See Also
    --------
    get_rnx_eost : Downloads RINEX files from the EOS Strasbourg server.
    singugins_run : Runs the GINS processing workflow.

    Examples
    --------
    >>> run_testgins()

    >>> run_testgins(results_folder="/custom/results", rnxs_folder="/custom/rnxs")

    >>> run_testgins(no_download_rnxs=True, rnxs_folder="/existing/rnxs")
    """

    # If no directory for RINEX files is provided, use the default path
    if not rnxs_folder:
        rnxs_folder = "/root/020_BDRNX/rnxs_testGINS_from_EOST"

    # Download or read RINEX files for the specified years
    if not no_download_rnxs:
        get_rnx_eost(rnxs_folder, year_start=None, year_end=None)

    # If no results directory is provided, create a timestamped folder
    if not results_folder:
        timstp = utils.get_timestamp()
        results_folder = "/root/030_RESULTS/testGINS_" + timstp

    # Run the GINS process with the specified directories
    singugins_run(
        results_folder=results_folder,
        bdrnx_folder=rnxs_folder,
    )

    # The function does not return any value
    return None


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description=(__doc__ or "").strip())

    parser.add_argument(
        "-o",
        "--results_folder",
        type=str,
        default=None,
        help=(
            "Results directory. If omitted, a timestamped folder under "
            "'/root/030_RESULTS/' is created. (optional)"
        ),
    )
    parser.add_argument(
        "-r",
        "--rnxs_folder",
        type=str,
        default=None,
        help=(
            "Directory where input RINEX files will be downloaded/read. "
            "If omitted, defaults to '/root/020_BDRNX/rnxs_testGINS_from_EOST'. (optional)"
        ),
    )

    parser.add_argument(
        "-n",
        "--no_download_rnxs",
        action="store_true",
        help=(
            "If set, the script will not download RINEX files from EOST server"
            "and will use existing files in the specified rnxs_folder."
        ),
    )

    args = parser.parse_args(argv)
    run_testgins(
        results_folder=args.results_folder,
        rnxs_folder=args.rnxs_folder,
        no_download_rnxs=args.no_download_rnxs,
    )

if __name__ == "__main__":
    main()
