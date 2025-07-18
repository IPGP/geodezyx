#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 15/03/2025 11:45:18

@author: psakic
"""

from spytgins.files_classes import StationFile
from spytgins import spotgins_utils

def create_stationfile_sitelogs(sitelogs, masterfile=False, outputfile=None, coords_from_master=False):
    """
    Convert a list of sitelogs to a station file

    Parameters
    ----------
    sitelogs : list of str
        List of sitelog files
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

    Returns
    -------
    str
        Path to the output file
    StationFile
        The station file object
    """
    masterfile, _ = spotgins_utils.get_spotgins_files(masterfile, None)
    sf = StationFile.from_sitelogfiles(sitelogs, masterfile)
    outputfile1 = sf.write(outputfile, coords_from_master=coords_from_master)
    return outputfile1, sf


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Convert a list of sitelogs to a station file')
    parser.add_argument('-s','--sitelogs', nargs='+', help='List of sitelog files', required=True)
    parser.add_argument('-m', '--masterfile', help='Path to the master file', required=False, default=None)
    parser.add_argument('-o', '--outputfile', help='Path to the output file', required=False, default=None)
    parser.add_argument('-cm', '--coords_from_master',
                        help='Get the station coordinates from the master file rather than the station file.',
                        action='store_true')

    args = parser.parse_args()

    create_stationfile_sitelogs(sitelogs=args.sitelogs,
                                masterfile=args.masterfile,
                                outputfile=args.outputfile,
                                coords_from_master=args.coords_from_master)

if __name__ == '__main__':
    main()
