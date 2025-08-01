#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from  spytgins.files_classes import MasterStationFile, StationFile
from spytgins import spotgins_utils
import os
from datetime import timedelta as td
import argparse

def update_station_file(stafile_new, stafile_path=None, master_path=None, output_dir=None):
    """
    Update the GINS "classic" station file with the new metadata from the logfiles descriibed
    in the SPOTGINS master file.

    Parameters
    ----------
    stafile_new : str or StationFile object
        the GINS "classic" station file with new values.
        Can be a path (as str) or directly a StationFile object.
    stafile_path : str
        Path to the GINS "classic" station file to be updated.
        If not provided, it will use the default file based
        on the environnement variable $SPOTGINS.
        Default is None
    master_path : str
        Path to the SPOTGINS master station file.
        If not provided, it will use the default file based
        on the environnement variable $SPOTGINS.
        Default is None
    output_dir : str
        The directory where the updated station file will be saved.
        If not provided, it will use the directory of the station file to be updated.
        Default is None

    Returns
    -------
    str
        Path to the updated station file.
    StationFile
        The updated station file object.
    """

    master_path, stafile_path = spotgins_utils.get_spotgins_files(master_path, stafile_path)

    if not output_dir:
        output_dir = os.path.dirname(stafile_path)

    mast = MasterStationFile.from_file(master_path)
    sta = StationFile.from_file(stafile_path, master_path)

    # station file with the up-to-date metadata (only)
    if type(stafile_new) is str:
        sta_new = StationFile.from_file(stafile_new, master_path)
    else:
        sta_new = stafile_new

    incons_path = os.path.join(output_dir, 'inconsistencies_stafile.txt')
    incons = open(incons_path, 'w')
    cpt = 0
    updates = []
    for acro9 in sta_new.stations.index:
        cpt += 1
        print('DBUG: update ' + acro9)

        if acro9 in sta.data.keys():
            # metadata in the current station file for the station acro9
            old_data = sta.data[acro9]
            last_change_n = list(old_data.keys())[-1]
            last_change = old_data[last_change_n]
        else:
            print('INFO: ' + acro9 + ' not in old station file. We add it.')
            sta.stations.loc[acro9] = sta_new.stations.loc[acro9]
            sta.data[acro9] = sta_new.data[acro9]
            continue

        # metadata in the new station file for the station acro9
        new_data = sta_new.data[acro9]

        # Comparing current and new station files. Inconsistency before last change?
        start = dt(2000, 1, 1)
        end = last_change['start']
        date_dt = start
        inconst = False

        while date_dt < end:
            # Receiver, antenna and eccentricities in the current station file
            (rec_old, ant_old, ecc_old) = spotgins_utils.extract_from_stafile(sta, date_dt, acro9)
            # Receiver, antenna and eccentricities in the new station file
            (rec_new, ant_new, ecc_new) = spotgins_utils.extract_from_stafile(sta_new, date_dt, acro9)

            if rec_old is not None and rec_new is not None:
                if rec_old != rec_new:
                    # different receiver
                    incons.write(acro9 + ' ' + str(date_dt.date()) + ' REC ' + rec_old + ' ' + rec_new + '\n')
                    inconst = True
            if ant_old is not None and ant_new is not None:
                if ant_old[0] != ant_new[0] or ant_old[1] != ant_new[1]:
                    # different antenna
                    incons.write(acro9 + ' ' + str(date_dt.date()) + ' ANT ' + str(ant_old) + ' ' + str(ant_new) + '\n')
                    inconst = True
            if ecc_old is not None and ecc_new is not None:
                if ecc_old != ecc_new:
                    # different eccentricities
                    incons.write(acro9 + ' ' + str(date_dt.date()) + ' ECC ' + str(ecc_old) + ' ' + str(ecc_new) + '\n')
                    inconst = True

            # 30 days increment for saving processing time
            date_dt += td(days=30)

        if not inconst:
            # pas d'incoherence avant la derniere date. On peut mettre a jour en remplacant TOUT l'actuel par TOUT le nouveau.
            sta.data[acro9] = sta_new.data[acro9]

        else:
            print('WARN: Inconsistencies for ' + acro9 + '. Check summary file inconsistencies_stafile.txt')

    incons.close()

    # Wrinting the updated global station file
    stafile_out_path = os.path.join(output_dir, os.path.basename(stafile_path) + ".updated.new")
    sta.write(outname=stafile_out_path)

    return stafile_out_path, sta


def main():
    parser = argparse.ArgumentParser(
        description="Update the GINS station file with new metadata from logfiles."
    )
    parser.add_argument('-n', '--stafile_new',
                        required=True, help='Path to the new GINS station file.')
    parser.add_argument('-s', '--stafile_path',
                        help='Path to the existing GINS station file to be updated.', default=None)
    parser.add_argument('-m', '--master_path',
                        help='Path to the SPOTGINS master station file.', default=None)
    parser.add_argument('-o', '--output_dir',
                        help='Directory to save the updated station file.', default=None)

    args = parser.parse_args()

    updated_file_path, _ = update_station_file(
        stafile_new=args.stafile_new,
        stafile_path=args.stafile_path,
        master_path=args.master_path,
        output_dir=args.output_dir,
    )

    if updated_file_path:
        print(f"INFO: Updated station file saved to: {updated_file_path}")
    else:
        print("EROR: Failed to update the station file.")


if __name__ == "__main__":
    main()
