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


def singugins_run(
    results_folder,
    specific_sites=[],
    start_epoch=dt.datetime(2000, 5, 3),
    end_epoch=dt.datetime(2099, 12, 31),
    bdrnx_folder="/root/020_BDRNX/",
    spotgins_run_kwargs={},
):
    """
    Run the SPOTGINS process within a SINGUGINS container.

    Parameters
    ----------
    results_folder : str
      Path to the folder where the results will be stored.
      we strongly recommed to use a subfolder within
      the default SINGUGINS results folder
      `/root/030_RESULTS/<your_subfolder>`
    specific_sites : list, optional
      List of specific RINEXs site IDs to include in the process.
       Default is an empty list.
    start_epoch : datetime, optional
      Start date and time for the selected RINEXs. Default is May 3, 2000.
    end_epoch : datetime, optional
      End date and time for the selected RINEXs. Default is December 31, 2099.
    bdrnx_folder : str, optional
      Path to the folder containing the RINEX files. Default is "/root/020_BDRNX/".
    spotgins_run_kwargs : dict, optional
      Additional keyword arguments to pass to the `spotgins_run` function.
      Default is an empty dictionary.

    Returns
    -------
    None
    """
    rnxs_lis = opera.rinex_finder(
        bdrnx_folder,
        specific_sites=specific_sites,
        start_epoch=start_epoch,
        end_epoch=end_epoch,
    )

    opera.spotgins_run(
        rnxs_path_inp=rnxs_lis,
        results_folder_inp=results_folder,
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

    args = parser.parse_args()

    # Convert spotgins_run_kwargs from JSON string to dictionary
    import json

    spotgins_run_kwargs = json.loads(args.spotgins_run_kwargs)

    singugins_run(
        results_folder=args.results_folder,
        specific_sites=args.specific_sites,
        start_epoch=args.start_epoch,
        end_epoch=args.end_epoch,
        bdrnx_folder=args.bdrnx_folder,
        spotgins_run_kwargs=spotgins_run_kwargs,
    )


if __name__ == "__main__":
    main()
