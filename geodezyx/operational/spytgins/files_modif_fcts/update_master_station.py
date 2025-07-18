#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from spytgins.files_classes import MasterStationFile
import argparse
import os
import numpy as np
from spytgins import spotgins_utils
from geodezyx import conv


def update_master_station(stations_new, ac, master_path=None, data_source='unknown', output_dir=None, xyz_input=False):
    """
    Update the master station file with new stations
    
    Parameters
    ----------
    
    stations_new : str or pandas.DataFrame
        List of new stations to add to the master station file. 
        Needs to contain 4char OR 9char acronym, long, lat, h and associated logfile, separated with spaces of tabs.
        Can be a path (as str) or directly a pandas.DataFrame with columns ['acro', 'long', 'lat', 'h', 'log'].

    ac : str
        Acronym of the AC responsible for the stations added

    master_path : str 
        The path of the master station file to update.
        If not provided, it will use the default file based on the environnement variable $SPOTGINS.
        Default is None

    data_source : str
        Data source
        Default is 'unknown'
        
    output_dir : str
        The directory where the updated master station file will be saved.
        If not provided, it will use the directory of the master file.
        Default is None

    xyz_input : bool
        If True, the input stations_new dataframe is expected to contain 'x', 'y', 'z' ,
        columns instead of 'long', 'lat', 'h'.
        Default is False.
        
    Returns
    -------
    str
        Path to the updated master station file.
    MasterStationFile
        The updated master station file
    """

    ### get the new stations values
    if type(stations_new) is str:
        if not os.path.isfile(stations_new):
            raise FileNotFoundError(f'Unable to find the new stations list {stations_new}')
        else:
            stations_new = pd.read_csv(stations_new, sep=r'\s+', header=None)
    else:
        stations_new = stations_new

    if xyz_input:
        cols_name = ['acro', 'x', 'y', 'z']
    else:
        cols_name = ['acro', 'long', 'lat', 'h']

    # a log column is provided
    if len(stations_new.columns) == 5:
        cols_name = cols_name + ['log']
        ## set the correct columns names
        stations_new.columns = cols_name
    # no log column is provided
    else:
        stations_new.columns = cols_name
        stations_new['log'] = 'None'

    ## the function works mainly with the 'long', 'lat', 'h' values, we convert it ASAP
    if xyz_input:
        llh_out = conv.XYZ2GEO_vector(stations_new[['x', 'y', 'z']].values)
        stations_new['lat'] = llh_out[:, 0]
        stations_new['long'] = llh_out[:, 1]
        stations_new['h'] = llh_out[:, 2]

    print(stations_new)

    ### get the master station values
    master_path, _ = spotgins_utils.get_spotgins_files(master_path, None)
    master = MasterStationFile.from_file(master_path)

    # get the output  directory
    if not output_dir:
        output_dir = os.path.dirname(master_path)

    for irow, row in stations_new.iterrows():
        # print(f'DBUG: {acro} in master file')

        (acro, long, lat, h, log) = row[['acro', 'long', 'lat', 'h', 'log']].values

        if len(acro) == 4:
            acro9 = spotgins_utils.build_longname(acro, long, lat)
        else:
            acro9 = acro
            acro = acro9[0:4]

        isnew = True

        for (name, latm, longm) in master.data[['NAME', 'LATITUDE', 'LONGITUDE']].values:
            if acro == name[0:4]:
                dist = spotgins_utils.haversine(long, lat, float(longm), float(latm))
                if dist < 0.5:
                    isnew = False
                    print(f'WARN: {acro9} already exists in the master file')

        if isnew:
            print(f'INFO: {acro9} is new and added in the master file')
            id_ = spotgins_utils.build_id(long, lat, master.data['ID'].values)
            (x, y, z) = spotgins_utils.geo2cart(lat * np.pi / 180., long * np.pi / 180., h)
            logpath = 'None' if log == 'None' else log
            master.add_entry(ac, acro9, id_, lat, long, h, x, y, z, logpath, data_source)

    master.write_data_source = True

    master_out_path = os.path.join(output_dir, os.path.basename(master_path) + ".new")
    master.write(outname=master_out_path)

    return master_out_path, master


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', dest='newlist', required=True,
                        help='List of new stations to add to the master station file.'
                             'Needs to contain 4char OR 9char acronym, long, lat, h and associated logfile,'
                             'Accept x y z instead of long lat h if -x is set,'
                             'separated with spaces of tabs')
    parser.add_argument('-ac', dest='ac', required=True,
                        help='Acronym of the AC responsible of the stations added')
    parser.add_argument('-m', dest='master', help='Initial master station file to update')
    parser.add_argument('-ds', dest='ds', help='Data source (default: %(default)s)', default='unknown')
    parser.add_argument('-x', dest='xyz', action='store_true',
                        help='Input file contains x, y, z instead of long, lat, h')

    args = parser.parse_args()

    update_master_station(master_path=args.master,
                          stations_new=args.newlist,
                          ac=args.ac,
                          data_source=args.ds,
                          xyz_input=args.xyz)


if __name__ == '__main__':
    main()
