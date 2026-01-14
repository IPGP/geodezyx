#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 05/01/2026 19:12:11

@author: psakic
"""

from geodezyx.operational.gins_runner import get_rnx_eost
from geodezyx.operational.gins_runner.singugins_run import singugins_run
from geodezyx import utils

def run_testgins(results_folder=None, rnxs_folder=None, download_rnxs=True, gins_run=True):
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
    download_rnxs : bool, optional
        If True, the script will not download RINEX files from EOST server
        and will use existing files in the specified rnxs_folder.
        Default is True.
    gins_run : bool, optional
        If True, the GINS processing workflow will be executed.
        Default is True.

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
    if download_rnxs:
        get_rnx_eost(rnxs_folder, year_start=None, year_end=None)

    # If no results directory is provided, create a timestamped folder
    if not results_folder:
        timstp = utils.get_timestamp()
        results_folder = "/root/030_RESULTS/testGINS_" + timstp

    # Run the GINS process with the specified directories
    if gins_run:
        singugins_run(
            results_folder=results_folder,
            bdrnx_folder=rnxs_folder,
            quick_mode=False
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
        "-nd",
        "--no_download_rnxs",
        action="store_true",
        help=(
            "If set, the script will not download RINEX files from EOST server"
            "and will use existing files in the specified rnxs_folder."
        ),
    )

    parser.add_argument(
        "-ng",
        "--no_gins_run",
        action="store_true",
        help=(
            "If set, the GINS processing workflow will not be executed."
        ),
    )

    args = parser.parse_args(argv)
    run_testgins(
        results_folder=args.results_folder,
        rnxs_folder=args.rnxs_folder,
        download_rnxs=not args.no_download_rnxs,
        gins_run=not args.no_gins_run,
    )

if __name__ == "__main__":
    main()
