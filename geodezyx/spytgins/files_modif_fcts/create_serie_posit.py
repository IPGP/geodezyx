#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 17/03/2025 19:10:52

@author: psakic
"""
import os.path

from geodezyx import utils as gzyx_utils
from geodezyx import files_rw as gzyx_files_rw
from geodezyx import time_series as gzyx_time_series

from spytgins.files_classes import MasterStationFile
from spytgins import spotgins_utils
import argparse


def create_serie_posit(ippp_files_list, output_dir, ac, master_path=None, data_src="unknown"):
    """
    Create a bunch of SPOTGINS's time series files
    from a list of GINS's solutions IPPP files.

    Several stations can be provided at once.

    Parameters
    ----------
    ippp_files_list : list
        List of GINS's solutions IPPP files to be processed.
    output_dir : str
        Directory where the position files will be saved.
    ac : str
        Analysis center identifier.
    master_path : str, optional
        Path to the master file.
        If not provided, it will use the default file based on the environnement variable $SPOTGINS.
        Default is None.
    data_src : str, optional
        Source of the data. Default is "unknown".

    Returns
    -------
    list
        List of the created position files.
    """

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    d = gzyx_files_rw.read_gins_solution_multi(ippp_files_list)
    master_path, _ = spotgins_utils.get_spotgins_files(master_path, None)
    mas = MasterStationFile.from_file(master_path)
    out = []
    for name, ts in d.items():
        #xyzcol = mas.find_entry(name, "NAME", ['X_position', 'Y_position', 'Z_position'])
        #flhcol = mas.find_entry(name, "NAME", ['X_position', 'Y_position', 'Z_position'])

        xyzcol = ['X_position', 'Y_position', 'Z_position']
        flhcol = ['LATITUDE', 'LONGITUDE', 'HEIGHT']
        
        xyzarr = mas.find_entry(name, "NAME", xyzcol)
        flharr = mas.find_entry(name, "NAME", flhcol)

        xref, yref, zref = xyzarr
        enu = gzyx_time_series.Point(xref, yref, zref, initype="XYZ")
        ts.rm_duplicat_pts("XYZ")
        ts.ENUcalc(enu)
        enu.F, enu.L, enu.H = flharr

        out_export = gzyx_time_series.export_ts_as_spotgins(ts,
                                                            output_dir,
                                                            ac=ac,
                                                            data_src=data_src)
        out.append(out_export)

    return out


def main():
    parser = argparse.ArgumentParser(
        description="Create a series of position files from a list of GINS's solutions IPPP files."
    )
    parser.add_argument('-i', '--ippp_files_folder', nargs='+', required=True,
                        help='Folders where are stored the GINS\'s solutions IPPP files to be processed.')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Directory where the position files will be saved.')
    parser.add_argument('-ac', required=True,
                        help='Analysis center identifier.')
    parser.add_argument('-m', '--master_path', default=None,
                        help='Path to the master file. Default is None.')
    parser.add_argument('-ds', '--data_src', default="unknown",
                        help='Source of the data. Default is "unknown".')

    args = parser.parse_args()

    ippp_files_list = []
    for fol in args.ippp_files_folder:
        ippp_files_list.extend(gzyx_utils.find_recursive(fol, "*IPPP*"))

    created_files = create_serie_posit(
        ippp_files_list=ippp_files_list,
        output_dir=args.output_dir,
        ac=args.ac,
        master_path=args.master_path,
        data_src=args.data_src
    )

if __name__ == "__main__":
    main()
