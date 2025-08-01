#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 16/03/2025 10:49:09

@author: psakic
"""

from  spytgins.files_classes import MasterStationFile
import os
import shutil
from datetime import datetime as dt
import numpy as np
from spytgins import spotgins_utils
import argparse

def update_master_sitelogs(local_sitelogs, ac, spotgins_logdir=None, master_path=None, output_dir=None):
    """
    Update the master station file with new sitelogs.

    Parameters
    ----------
    local_sitelogs : str or list of str
        Path new sitelog files. Can be a path (as str) or a list of paths.
    ac : str
        Acronym of the AC responsible for the stations.
    spotgins_logdir : str, optional
        Path to the SPOTGINS sitelog directory. Defaults to None.
    master_path : str, optional
        Path to the master station file. Defaults to None.
    output_dir : str, optional
        Directory to save the updated master station file. Defaults to None.

    Returns
    -------
    str
        Path to the updated master station file.
    MasterStationFile
        The updated master station file.
    """
    # Get the master station values
    master_path, _ = spotgins_utils.get_spotgins_files(master_path, None)
    master = MasterStationFile.from_file(master_path)
    if not spotgins_logdir:
        spotgins_logdir = os.path.join(spotgins_utils.get_spotgins_dir(), 'metadata', 'logfiles')

    # get the output  directory
    if not output_dir:
        output_dir = os.path.dirname(master_path)

    # List of updated logfiles
    newout_path = os.path.join(output_dir, 'new_logfiles.list')
    newout = open(newout_path, 'w')

    # new logfiles?
    if not type(local_sitelogs) is list:
        local_logs_lis = [local_sitelogs]
    else:
        local_logs_lis = local_sitelogs

    local_logdate_lis = [dt.strptime(l.split('_')[1].split('.')[0], '%Y%m%d') for l in local_logs_lis]

    for (acro9, log) in master.data[master.data['AC'] == ac][['NAME', 'LOGFILE_NAME']].values:

        if np.isnan(log):
            spotgins_logdate = dt(1900, 1, 1)
        else:
            try:
                spotgins_logdate = dt.strptime(log.split('_')[1].split('.')[0], '%Y%m%d')
            except:
                print(f"Unable to extract log date for {acro9}")
                continue

        for local_log, local_logdate in zip(local_logs_lis, local_logdate_lis):

            if not (local_log.lower().startswith(acro9.lower()) and local_logdate > spotgins_logdate):
                continue

            print(f"New logfile for {acro9}: {log} -> {local_log}")

            # Update in master
            index = master.data[master.data['NAME'] == acro9].index[0]
            master.data.loc[index, 'LOGFILE_NAME'] = local_log

            # Copy in SPOTGINS logdir
            shutil.copyfile(f"{local_sitelogs}/{local_log}", f"{spotgins_logdir}/{local_log}")
            newout.write(f"{acro9}\n")

    # writing the updated master station file.
    newout.close()
    master_out_path = os.path.join(output_dir, os.path.basename(master_path) + ".new")
    master.write(outname=master_out_path)

    return master_out_path, master


def main():
    parser = argparse.ArgumentParser(
        description="Update the master station file with new sitelogs."
    )
    parser.add_argument('-l','--local_sitelogs', help='Path of the new sitelog files.',
                        nargs='+', required=True)
    parser.add_argument('-ac', help='Acronym of the AC responsible for the stations.', required=True)
    parser.add_argument('-s', '--spotgins_logdir', help='Path to the SPOTGINS sitelog directory.', default=None)
    parser.add_argument('-m', '--master_path', help='Path to the master station file.', default=None)
    parser.add_argument('-o', '--output_dir', help='Directory to save the updated master station file.', default=None)

    args = parser.parse_args()

    updated_file_path, _ = update_master_sitelogs(
        local_sitelogs=args.local_logdir,
        ac=args.ac,
        spotgins_logdir=args.spotgins_logdir,
        master_path=args.master_path,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
