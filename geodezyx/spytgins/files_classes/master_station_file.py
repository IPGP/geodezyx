#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from datetime import datetime as dt
import numpy as np
import os
from spytgins import spotgins_utils


class MasterStationFile:
    """
    Class to handle SPOTGINS master station files
    """

    def __init__(self):

        self.filename = ''
        self.write_data_source = True
        self.nb_stations = 0
        self.data = pd.DataFrame(
            columns=['AC', 'NAME', 'ID', 'LATITUDE', 'LONGITUDE', 'HEIGHT', 'X_position', 'Y_position', 'Z_position',
                     'LOGFILE_NAME', 'DATA_SOURCE'])
        self.headerlines = ['#STATION MASTER FILE', '#  -COORDINATES OF STATIONS CONSIDERED', '#  -WITH LOGFILES INFOS',
                            '# Last update: ' + str(dt.today()), '#',
                            '#AC    NAME       ID       LATITUDE    LONGITUDE        HEIGHT       X_position       Y_position       Z_position  LOGFILE_NAME            DATA_SOURCE']

    def add_entry(self, ac, name, id, lat, long, h, x, y, z, log, repo='unknown'):
        """
        Add a new entry to a master station file
        """

        self.data.loc[self.nb_stations] = np.array([ac, name, id, lat, long, h, x, y, z, log, repo])

        self.nb_stations += 1

    def find_entry(self, value, col_inp, col_out=None):
        """
        Find an entry (value) in a specific column (col_inp) in the master station file.

        Parameters
        ----------
        value : str
            The value to search for in the specified column.
        col_inp : str
            The column name to search within.
        col_out : str, optional
            The column name from which to return the value. If None, the entire row is returned.

        Returns
        -------
        any
            The value from the specified output column if col_out is provided, otherwise the entire row.
            Returns None if no entry is found.
        """
        site_found = self.data[self.data[col_inp] == value]

        if len(site_found) == 0:
            print(f"WARN: No {value} entry found. Add site to master file first!")
            return None
        elif len(site_found) > 1:
            print(f"WARN: Several {value} entries in master file!")
        else:
            pass
        if not col_out:
            return site_found.values[0]
        else:
            return site_found[col_out].values[0]

    def write(self, outname=None):
        """
        Write the master station file
        """

        if outname is None:
            if spotgins_utils.get_spotgins_dir():
                outname = os.path.join(spotgins_utils.get_spotgins_dir(), 'metadata/stations', self.filename + '.new')
            else:
                outname = self.filename + '.new'

        self.data = self.data.sort_values(['NAME'], axis=0)

        print('writing ' + outname + ', NB CWD: ' + os.getcwd())
        out = open(outname, 'w')

        for hl in self.headerlines:
            out.write(hl + '\n')

        for (ac, name, id, lat, long, h, x, y, z, log, repo) in self.data.values:
            lat = float(lat)
            long = float(long)
            h = float(h)
            x = float(x)
            y = float(y)
            z = float(z)

            # string formatting
            # lat = f'{sign(lat)}{np.abs(lat):09.6f}'
            # long = f'{sign(long)}{np.abs(long):010.6f}'
            # h = f'{sign(h)}{np.abs(h):011.6f}'
            # print(lat)
            lat = f'{lat:9.6f}'
            # print (lat)

            long = f'{long:10.6f}'
            h = f'{h:11.6f}'

            x = f'{x:15.6f}'
            y = f'{y:15.6f}'
            z = f'{z:15.6f}'

            if str(log) == 'nan':
                log = 'None'

            if self.write_data_source:

                out.write(' ' + ac.rjust(4) + '  ' + name + '  ' + str(id) + '  ' + lat.rjust(10) + '  ' + long.rjust(
                    11) + '  ' + h.rjust(12) + '  ' + x + '  ' + y + '  ' + z + '  ' + log.ljust(
                    22) + '  ' + repo + '\n')

            else:
                out.write(' ' + ac.rjust(4) + '  ' + name + '  ' + str(id) + '  ' + lat.rjust(10) + '  ' + long.rjust(
                    11) + '  ' + h.rjust(12) + '  ' + x + '  ' + y + '  ' + z + '  ' + log.ljust(22) + '\n')

        out.close()
        print("INFO: Master station file written to " + outname)
        return outname

    def check_doublons(self):
        """
        Checks duplicates in names a 5 first chars of ids
        """
        # names
        name = ''
        for name_tmp in self.data['NAME'].values.sort():
            if name_tmp == name:
                # duplicates
                print('WARNING: duplicate name: ' + name)
            name_tmp = name

        # ids
        code = ''
        for code_tmp in self.data['ID'].str.sort():
            if code_tmp == code:
                # duplicates
                print('WARNING: duplicate id: ' + code)
            code = code_tmp

    @classmethod
    def from_file(cls, masterstationfile):
        """
        Creates a master station object from an existing file.

        Parameters
        ----------
        cls : type
            The class type to instantiate.
        masterstationfile : str
            Path to the master station file.

        Returns
        -------
        MasterStationFile
            An instance of MasterStationFile populated with data from the file.
        """

        if not os.path.isfile(masterstationfile):
            print('FATAL: Unable to find the file ' + masterstationfile)

        print("INFO: Import master station file: " + masterstationfile)

        master = MasterStationFile()

        master.data = pd.read_csv(masterstationfile, header=None, sep=r'\s+', comment='#',
                                  names=['AC', 'NAME', 'ID', 'LATITUDE', 'LONGITUDE', 'HEIGHT', 'X_position',
                                         'Y_position', 'Z_position', 'LOGFILE_NAME', 'DATA_SOURCE'])

        master.nb_stations = len(master.data)
        master.filename = masterstationfile

        return master

def sign(num):
    """
    Extract the sign fo a float
    :return '+' or '-'
    """
    if num >= 0:
        sign = '+'
    else:
        sign = '-'

    return sign
