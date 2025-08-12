#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 21/03/2025 21:57:04

@author: psakic
"""

import geodezyx.operational as opera
import geodezyx.conv as conv
import datetime as dt
import argparse
import os


def singugins_run(
    results_folder,
    specific_sites=[],
    start_epoch=dt.datetime(2000, 5, 3),
    end_epoch=dt.datetime(2099, 12, 31),
    bdrnx_folder="/root/020_BDRNX/",
    nprocs=8,
    no_updatebd=False,
    no_concat_orb_clk=False,
    verbose=False,
    force=False,
    spotgins_run_kwargs={},
    no_rnx2=False,
    no_rnx3=False,
):
    """
    Run the SPOTGINS process within a SINGUGINS container.

    Parameters
    ----------
    results_folder : str
        Path to the folder where the results will be stored.
        We strongly recommend using a subfolder within
        the default SINGUGINS results folder
        `/root/030_RESULTS/<your_subfolder>`.
    specific_sites : list, optional
        List of specific RINEXs site IDs to include in the process.
        Default is an empty list.
    start_epoch : datetime, optional
        Start date and time for the selected RINEXs. Default is May 3, 2000.
    end_epoch : datetime, optional
        End date and time for the selected RINEXs. Default is December 31, 2099.
    bdrnx_folder : str, optional
        Path to the folder containing the RINEX files. Default is "/root/020_BDRNX/".
    nprocs : int, optional
        Number of processes to use. Default is 8.
    no_updatebd : bool, optional
        Flag to indicate whether to update the database. Default is False.
    no_concat_orb_clk : bool, optional
        Flag to indicate whether to concatenate the orbit and
        clock files prior to the main processing
        Default is False.
    verbose : bool, optional
        Verbose mode. If True, additional information will be printed.
    force : bool, optional
        Flag to force the process to run even if the results already exist.
    spotgins_run_kwargs : dict, optional
        Additional keyword arguments to pass to the `spotgins_run` function.
        Default is an empty dictionary.
    no_rnx2 : bool, optional
        If True, RINEX2 files will not be processed.
        Default is False (RINEX2 files will be processed).
    no_rnx3 : bool, optional
        If True, RINEX3 files will not be processed.
        Default is False (RINEX3 files will be processed).

    Returns
    -------
    None
    """

    import numpy as np
    srt_epoch_ok = min([start_epoch, end_epoch])
    end_epoch_ok = max([start_epoch, end_epoch])

    rnxs_lis = opera.rinex_finder(
        bdrnx_folder,
        specific_sites=specific_sites,
        start_epoch=srt_epoch_ok,
        end_epoch=end_epoch_ok,
        short_name=not no_rnx2,
        long_name=not no_rnx3,
    )

    opera.spotgins_run(
        rnxs_path_inp=rnxs_lis,
        results_folder_inp=results_folder,
        nprocs=nprocs,
        no_updatebd=no_updatebd,
        no_concat_orb_clk=no_concat_orb_clk,
        verbose=verbose,
        force=force,
        **spotgins_run_kwargs,
    )

    return

def main():
    parser = argparse.ArgumentParser(
        description="Run the SPOTGINS process within a SINGUGINS container."
    )
    parser.add_argument(
        "-o",
        "--results_folder",
        type=str,
        help="Path to the folder where the results will be stored.",
        required=True,
    )
    parser.add_argument(
        "-l",
        "--specific_sites",
        type=str,
        nargs="*",
        default=[],
        help="List of specific RINEXs site IDs to include in the process.",
    )
    parser.add_argument(
        "-s",
        "--start_epoch",
        type=lambda s: conv.date_pattern_2_dt(s),
        default=dt.datetime(2000, 5, 3),
        help="Start date and time for the selected RINEXs.",
    )
    parser.add_argument(
        "-e",
        "--end_epoch",
        type=lambda s: conv.date_pattern_2_dt(s),
        default=dt.datetime(2099, 12, 31),
        help="End date and time for the selected RINEXs.",
    )
    parser.add_argument(
        "-r",
        "--bdrnx_folder",
        type=str,
        default="/root/020_BDRNX/",
        help="Path to the folder containing the RINEX files.",
    )
    parser.add_argument(
        "-k",
        "--spotgins_run_kwargs",
        type=str,
        default="{}",
        help="Additional keyword arguments to pass to the `spotgins_run` function (as a JSON string).",
    )

    parser.add_argument(
        "-p",
        "--nprocs",
        type=int,
        default=8,
        help="Number of processes to use.",
    )
    parser.add_argument(
        "-nu",
        "--no_updatebd",
        action="store_true",
        help="Flag to indicate whether to update the database.",
    )

    parser.add_argument(
        "-nc",
        "--no_concat_orb_clk",
        action="store_true",
        help="Flag to indicate whether to concatenate the"
             "orbit and clock files prior to the main processing",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbose mode.",
    )

    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Flag to force the process to run even if the results already exist.",
    )

    parser.add_argument(
        "-nr2",
        "--no_rnx2",
        action="store_true",
        help="If True, RINEX2 files will not be processed.",
    )

    parser.add_argument(
        "-nr3",
        "--no_rnx3",
        action="store_true",
        help="If True, RINEX3 files will not be processed.",
    )

    args = parser.parse_args()

    # Convert spotgins_run_kwargs from JSON string to dictionary
    import json

    spotgins_run_kwargs = json.loads(args.spotgins_run_kwargs)

    singugins_run(
        results_folder=os.path.abspath(args.results_folder),
        specific_sites=args.specific_sites,
        start_epoch=args.start_epoch,
        end_epoch=args.end_epoch,
        bdrnx_folder=args.bdrnx_folder,
        spotgins_run_kwargs=spotgins_run_kwargs,
        no_updatebd=args.no_updatebd,
        no_concat_orb_clk=args.no_concat_orb_clk,
        nprocs=args.nprocs,
        verbose=args.verbose,
        force=args.force,
        no_rnx2=args.no_rnx2,
        no_rnx3=args.no_rnx3,
    )


if __name__ == "__main__":
    main()
