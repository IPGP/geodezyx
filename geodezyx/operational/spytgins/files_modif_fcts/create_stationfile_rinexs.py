#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08/05/2025 09:28:19

@author: psakic
"""

from spytgins.files_classes import StationFile
from spytgins import spotgins_utils


def create_stationfile_rinexs(
    rinexs,
    masterfile=None,
    outputfile=None,
    coords_from_master=False,
    ninecharfile=None,
):
    """
    Convert a list of rinex files to a station file

    Parameters
    ----------
    rinexs : list of str
        List of rinex files
    masterfile : str
        Path to the master file
        If not provided, the master file will be searched in the SPOTGINS stations directory
    outputfile : str
        Path to the output file
        If not provided, the output file will be named 'station_file.dat.new'
        in the SPOTGINS stations directory
    coords_from_master : bool
        Get the station coordinates from the master file rather than the station file.
        Default is False
    ninecharfile : str
        Path to the nine-character file, a file that contains 9-char. site names
        If not provided, the nine-character file will be searched in the SPOTGINS stations directory

    Returns
    -------
    str
        Path to the output file
    StationFile
        The station file object
    """

    masterfile, _ = spotgins_utils.get_spotgins_files(masterfile, None)
    sf = StationFile.from_rinexs(rinexs, masterfile, ninecharfile)
    outputfile1 = sf.write(outputfile, coords_from_master=coords_from_master)
    return outputfile1, sf


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert a list of rinex files to a station file"
    )
    parser.add_argument(
        "-r", "--rinexs", nargs="+", help="List of rinex files", required=True
    )
    parser.add_argument(
        "-l",
        "--rinexs_list",
        help="Path to a file containing a list of RINEX files, one per line. "
        "This option is used as an alternative to providing the list directly with the -r option.",
        required=False,
    )
    parser.add_argument(
        "-m",
        "--masterfile",
        help="Path to the master file",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-o",
        "--outputfile",
        help="Path to the output file",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-n",
        "--ninecharfile",
        help="Path to the nine-character file, a file that contains 9-char. site names",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-cm",
        "--coords_from_master",
        help="Get the station coordinates from the master file rather than the station file.",
        action="store_true",
        required=False,
    )

    args = parser.parse_args()

    if args.rinexs_list:
        with open(args.rinexs_list, "r") as f:
            rinexs_list = [line.strip() for line in f.readlines()]
    else:
        rinexs_list = args.rinexs

    create_stationfile_rinexs(
        rinexs=rinexs_list,
        masterfile=args.masterfile,
        outputfile=args.outputfile,
        ninecharfile=args.ninecharfile,
        coords_from_master=args.coords_from_master,
    )


if __name__ == "__main__":
    main()
