#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
StationFile class to handle GINS "classic" station files

Nota Bene
---------
StationFile's 'data' attribute is stored in a dictionary with the following structure:

 'ELLY00USA': {0: {'start': datetime.datetime(2013, 2, 25, 9, 3),
   'end': datetime.datetime(2013, 6, 12, 20, 59, 59),
   'rec': 'TRIMBLE NETR9',
   'ant': ['TRM59800.00', 'SCIT', '5246354615', '0.0'],
   'ecc_u': [0.0083, 0.0],
   'ecc_n': [0.0, 0.0],
   'ecc_e': [0.0, 0.0],
   'source': 'log_file'},
  1: {'start': datetime.datetime(2013, 6, 12, 21, 0),
   'end': datetime.datetime(2016, 7, 19, 14, 45, 59),
   'rec': 'TRIMBLE NETR9',
   'ant': ['TRM59800.00', 'SCIT', '5246354615', '0.0'],
   'ecc_u': [0.0083, 0.0],
   'ecc_n': [0.0, 0.0],
   'ecc_e': [0.0, 0.0],
   'source': 'log_file'},
  2: {'start': datetime.datetime(2016, 7, 19, 14, 46),
   'end': datetime.datetime(2017, 4, 8, 21, 25, 59),
   'rec': 'TRIMBLE NETR9',
   'ant': ['TRM59800.00', 'SCIT', '5246354615', '0.0'],
   'ecc_u': [0.0083, 0.0],
   'ecc_n': [0.0, 0.0],
   'ecc_e': [0.0, 0.0],
   'source': 'log_file'},
 'AV0900USA': {0: {'start': datetime.datetime(2004, 5, 5, 0, 0),
   'end': datetime.datetime(2004, 10, 4, 17, 12, 59),
   'rec': 'TRIMBLE NETRS',
   'ant': ['TRM29659.00', 'SCIT', '0220314089', '0.0'],
   'ecc_u': [0.0083, 0.0],
   'ecc_n': [0.0, 0.0],
   'ecc_e': [0.0, 0.0],
   'source': 'log_file'},
  1: {'start': datetime.datetime(2004, 10, 4, 17, 13),
   'end': datetime.datetime(2007, 5, 15, 22, 53, 59),
   'rec': 'TRIMBLE NETRS',
   'ant': ['TRM29659.00', 'SCIT', '0220314089', '0.0'],
   'ecc_u': [0.0083, 0.0],
   'ecc_n': [0.0, 0.0],
   'ecc_e': [0.0, 0.0],
   'source': 'log_file'},
"""

import os.path
import re

import numpy as np
import pandas as pd
from spytgins.files_classes import OffsetsFile, MasterStationFile
from spytgins import spotgins_utils
from datetime import datetime as dt


class StationFile:
    """
    Class to handle GINS "classic" station files
    """

    def __init__(self, masterfile=None):
        """
        Input masterfile is recommeded but not mandatory for advanced usages
        """
        if not masterfile:
            self.masterfile = None
        else:
            self.masterfile = MasterStationFile.from_file(masterfile)

        self.stations = pd.DataFrame(
            columns=['id', 'X', 'sX', 'Y', 'sY', 'Z', 'sZ', 'VX', 'sVX', 'VY', 'sVY', 'VZ', 'sVZ',
                     'epoch', 'ref_frame'])
        self.label = 'GRG3'
        self.data = {}
        self.version = '1.0'
        self.radius = '0.6378137000000000E+07'
        self.iflat = '0.2982572221010000E+03'
        self.ref_frame = 'ITRF2020'
        self.headerlines = []

    def set_headerlines(self):
        """
        Set the header lines of the station file based
        on the attributes of the StationFile object
        to the headerlines attribute
        """
        self.headerlines = ['#VERSION ' + self.version,
                            '#RADIUS   ' + self.radius,
                            '#1/FLAT   ' + self.iflat,
                            '#' + self.ref_frame]

    def write(self, outname=None, station='all', coords_from_master=False):
        """
        Write the station file
        """

        if outname is None:
            if spotgins_utils.get_spotgins_dir():
                outname = os.path.join(spotgins_utils.get_spotgins_dir(), 'metadata/stations', 'station_file.dat.new')
            else:
                outname = 'station_file.dat.new'

        if station != 'all':
            # on n'ecrit le fichier sta que pour une seule station
            stations_towrite = self.stations[self.stations.index == station]
            if outname is None:
                outname = station + '.STA'
        else:
            stations_towrite = self.stations

        # Experimental if we want to sort the stations
        # stations_towrite.sort_index(inplace=True)
        # stations_towrite.sort_values(by='id', inplace=True)

        out = open(outname, 'w')
        for hl in self.headerlines:
            out.write(hl + '\n')

        for (acro, values) in stations_towrite.iterrows():

            acro9 = self.masterfile.find_entry(values[0], 'ID', 'NAME') if len(acro) == 4 else acro
            # print ('acro9 = ' + acro9
            acro = acro9[0:4]
            (id_, x, s_x, y, s_y, z, s_z, vx, s_vx, vy, s_vy, vz, s_vz, epoch, ref_frame) = values

            if coords_from_master:
                xyz_master = self.masterfile.find_entry(acro9, 'NAME', ['X_position', 'Y_position', 'Z_position'])
                if xyz_master is not None:
                    x, y, z = map(float, xyz_master)
                else:
                    print(f"EROR: no coords for {acro9} found in master file. Fallback to sitelog's ones!")

            x, y, z = float(x), float(y), float(z)
            s_x, s_y, s_z = float(s_x), float(s_y), float(s_z)
            vx, vy, vz = float(vx), float(vy), float(vz)
            s_vx, s_vy, s_vz = float(s_vx), float(s_vy), float(s_vz)

            pline = '   ' + str(int(float(id_))).zfill(5) + '      ' + acro9 + '                       '
            pline += ' mar_xyz' + ('%.6f' % x).rjust(16) + ' ' + '%.6f' % s_x + ('%.6f' % y).rjust(
                16) + ' ' + '%.6f' % s_y + ('%.6f' % z).rjust(16) + ' ' + '%.6f' % s_z
            pline += ' vit_xyz ' + ('%.5f' % vx).rjust(8) + ' ' + '%.4f' % s_vx + ' ' + ('%.5f' % vy).rjust(
                8) + ' ' + '%.4f' % s_vy + ' ' + ('%.5f' % vz).rjust(8) + ' ' + '%.4f' % s_vz
            pline += '  ' + epoch.strftime('%Y%m%d_%H%M%S') + '                 ' + ref_frame

            out.write('\n' + pline + '\n')

            # materiel
            cpt_neq = 0
            for neq in self.data[acro9].keys():
                eq = self.data[acro9][neq]
                end_use = dt(2099, 12, 31, 0, 0, 0) if eq['end'] > dt(2099, 12, 31, 0, 0, 0) else eq['end']
                end = end_use.strftime('%Y%m%d_%H%M%S')
                start = eq['start'].strftime('%Y%m%d_%H%M%S')
                rec = eq['rec']
                (ant, radome, ant_sn, ant_ori) = eq['ant']
                if isinstance(ant_ori, str):
                    ant_ori = ant_ori.split()[-1]  ### avoid parasit values at the begining
                elif np.isnan(ant_ori):
                    ant_ori = 0.
                (ecc_u, secc_u) = eq['ecc_u']
                (ecc_n, secc_n) = eq['ecc_n']
                (ecc_e, secc_e) = eq['ecc_e']
                source = eq['source']

                eline = ' ' + str(int(float(id_))).zfill(5) + str(cpt_neq).zfill(2) + ' ' + self.label + ' ' + acro9
                eline += '  ' + rec.upper().ljust(22)
                eline += 'ecc_une      ' + ('%.6f' % float(ecc_u)).rjust(10) + ' %.6f' % secc_u + '      ' \
                         + ('%.6f' % float(ecc_n)).rjust(10) + ' %.6f' % secc_n + '      ' \
                         + ('%.6f' % float(ecc_e)).rjust(10) + ' %.6f' % secc_e + '            '

                eline += (ant.upper().ljust(16) + radome.upper()[:4] + ' ' + ant_sn.strip()[0:17].ljust(17) +
                          ' ' + '%6.1f' % float(ant_ori) + ' ')
                eline += start + ' ' + end
                eline += ' ' + source

                out.write(eline + '\n')

                cpt_neq += 1

        out.close()
        print("INFO: Station file written to " + outname)
        return outname

    def to_offsets_file(self):
        offsetfile = OffsetsFile.from_station_file(self)
        offsetfile.write_amplitude = True
        offsetfile.write(outfile='offsets_change.dat')

    @classmethod
    def from_file(cls, stationfile, masterfile):
        """
        Create a StationFile object from a GINS station file

        Parameters
        ----------
        stationfile : str
            Path to the station file

        masterfile : str
            Path to the master station file

        Returns
        -------
        sfout : StationFile
            StationFile object

        """
        if not masterfile:
            print("INFO: no masterfile provided, but we use the default one")
            masterfile, _ = spotgins_utils.get_spotgins_files()

        sfout = StationFile(masterfile)
        print("INFO: Import station file " + stationfile)
        data = open(stationfile)
        for line in data:
            # print (line)
            if line[0:1] == '#':
                sfout.headerlines.append(line.strip())
                if '#VERSION' in line:
                    sfout.version = line.strip().split()[1]
                if '#RADIUS' in line:
                    sfout.radius = line.strip().split()[1]
                if '#1/FLAT' in line:
                    sfout.iflat = line.strip().split()[1]

            elif line.strip() == '':
                continue

            elif line[0:3] == '   ':
                # first line of station bloc
                ltab = line.strip().replace('mar_xyz', '').replace('vit_xyz', '').split()
                # domes = ltab[0] + ltab[1]
                id_ = int(float(ltab[0]))
                acro = ltab[1].strip()
                acro9 = sfout.masterfile.find_entry(id_, "ID", 'NAME') if len(acro) == 4 else acro

                # print('acro9 = #' + acro9 + '#')
                epoch = dt.strptime(ltab[14], '%Y%m%d_%H%M%S')
                last_end = dt(1900, 1, 1)
                sfout.stations.loc[acro9] = np.array([id_] + ltab[2:14] + [epoch] + [ltab[15]])
            else:
                # material line
                cpt_mat = line[6:8]
                acro = line[14:23].strip()
                acro9 = sfout.masterfile.find_entry(id_, "ID", 'NAME') if len(acro) == 4 else acro

                sfout.label = line[9:13]
                rec = line[25:46].strip()
                (ecc_u, secc_u, ecc_n, secc_n, ecc_e, secc_e) = [float(val) for val in line[62:127].split()]
                ant = line[141:156].strip()
                dome = line[157:161]
                ant_sn = line[162:179].strip()
                ant_ori = line[180:186].strip()
                if ant_ori == '':
                    ant_ori = '0.0'
                start = dt.strptime(line[187:202], '%Y%m%d_%H%M%S')
                end = dt.strptime(line[203:218], '%Y%m%d_%H%M%S')
                source = line.strip()[218:]

                if acro9 not in sfout.data.keys():
                    neq = 0
                    sfout.data[acro9] = {}

                if start > end:
                    print(f'WARNING: Start and end date issues for {acro9} (s={start} e={end})')
                    # return

                if start < last_end:
                    print(f'WARNING: Start date issues for {acro9} (s={start})')

                    # return
                last_end = end
                sfout.data[acro9][neq] = {'start': start,
                                          'end': end,
                                          'rec': rec,
                                          'ant': [ant, dome, ant_sn, ant_ori],
                                          'ecc_u': [ecc_u, secc_u],
                                          'ecc_n': [ecc_n, secc_n],
                                          'ecc_e': [ecc_e, secc_e],
                                          'source': source}

                neq += 1

        return sfout

    @classmethod
    def from_sitelogfiles(cls, sitelogfiles, masterfile=None):
        """
        Create a StationFile object from a single sitelogfile or list of sitelogfiles

        Parameters
        ----------

        sitelogfiles : list or str
            List of sitelogfiles paths
            Can alo handle a single sitelogfile path

        masterfile : str
            Path to the master station file

        Returns
        -------
        sfout : StationFile
            A StationFile object
        """

        if not masterfile:
            print("INFO: no masterfile provided, but we use the default one")
            masterfile, _ = spotgins_utils.get_spotgins_files()

        import rinexmod
        slgfil_lis = [sitelogfiles] if not isinstance(sitelogfiles, list) else sitelogfiles

        sfout = StationFile(masterfile)
        slinslg_stk = []
        ### read, stack and save the sitelogs as rinexmod MetaData
        for slg in slgfil_lis:
            # print("DBUG: Import sitelog file " + os.path.basename(slg))
            mdaslg = rinexmod.metadata.MetaData(slg)
            data_slg, sta_lin_slg = metadata2stationfile_data(mdaslg)
            sfout.data = {**sfout.data, **data_slg}
            slinslg_stk.append(sta_lin_slg)

        sfout.set_stations_from_lines(slinslg_stk)
        return sfout

    #### GAMIT PART
    @classmethod
    def from_gamit(cls, station_info_path, lfile_path=None,
                   force_fake_coords=False, masterfile=None,
                   ninecharfile=None, rev=False):
        """
        Create a StationFile object from GAMIT station.info & associated files.

        Parameters
        ----------
        station_info_path : str
            Path to the GAMIT station information file.
        lfile_path : str, optional
            Path to the L-file. Default is None.
        force_fake_coords : bool, optional
            If True, force the use of fake coordinates. Default is False.
        masterfile : str, optional
            Path to the SPOTGINS master station file. Default is None.
        ninecharfile : str, optional
            Path to the nine-character file, i.e. a list of 9char. codes.
            Default is None.
        rev : bool, optional
            reverse order of inputs for station.info filled from the
            newest to the oldest change
            The default is False.

        Returns
        -------
        sfout : StationFile
            A StationFile object created from the GAMIT station information.
        """
        import rinexmod

        if not masterfile:
            print("INFO: no masterfile provided, but we use the default one")
            masterfile, _ = spotgins_utils.get_spotgins_files()

        sfout = StationFile(masterfile)

        if not ninecharfile:
            print("INFO: no ninecharfile provided, but we use masterfile as fallback to get the 9char codes")
            ninecharfile = MasterStationFile.from_file(masterfile).data['NAME'].values

        mda_lis = rinexmod.rinexmod_api.gamit2mda_objs(station_info_path,
                                                       lfile_inp=lfile_path,
                                                       force_fake_coords=force_fake_coords,
                                                       ninecharfile_inp=ninecharfile,
                                                       rev=rev)
        slinslg_stk = []
        for mda in mda_lis:
            data_out, sta_lin_out = metadata2stationfile_data(mda)
            sfout.data = {**sfout.data, **data_out}
            slinslg_stk.append(sta_lin_out)

        sfout.set_stations_from_lines(slinslg_stk)
        return sfout

    #### RINEXPART
    @classmethod
    def from_rinexs(cls, rinex_paths, masterfile=None, ninecharfile=None):
        """
        Create a StationFile object from RINEX files.

        Parameters
        ----------
        cls : type
            The class type (StationFile).
        rinex_paths : list or str
            Path(s) to the RINEX files. Can be a single path or a list of paths.
        masterfile : str, optional
            Path to the SPOTGINS master station file. Default is None.
        ninecharfile : str, optional
            Path to the nine-character file, i.e., a list of 9-character codes. Default is None.

        Returns
        -------
        StationFile
            A StationFile object created from the provided RINEX files.
        """
        import rinexmod

        if not masterfile:
            print("INFO: no masterfile provided, but we use the default one")
            masterfile, _ = spotgins_utils.get_spotgins_files()

        sfout = StationFile(masterfile)

        if not ninecharfile:
            print("INFO: no ninecharfile provided, but we use masterfile as fallback to get the 9char codes")
            ninecharfile = MasterStationFile.from_file(masterfile).data['NAME'].values

        mdaobjs_lis = rinexmod.rinexmod_api.rinexs2mda_objs(rinex_paths, ninecharfile_inp=ninecharfile)
        mdaobjs_lis_merged, _ = rinexmod.rinexmod_api.group_mda(mdaobjs_lis)


        slinslg_stk = []
        for mda_site in mdaobjs_lis_merged:
            mda_site.merge_instrus()
            data_out, sta_lin_out = metadata2stationfile_data(mda_site)
            sfout.data = {**sfout.data, **data_out}
            slinslg_stk.append(sta_lin_out)

        sfout.set_stations_from_lines(slinslg_stk)

        return sfout

    ### GENERAL METHODS
    def set_stations_from_lines(self, slinslg_stk):
        """
        Stack stations lines created with metadata2stationfile_data
        to set self.stations attribute.

        Parameters
        ----------
        slinslg_stk : list
            List of station lines created with metadata2stationfile_data.

        Returns
        -------
        None
        """
        ### stack the station DataFrame
        stations = pd.DataFrame(slinslg_stk)
        stations.columns = ["acro"] + list(self.stations.columns)

        # manage the ID, when we work with a masterfile
        if self.masterfile:
            stations["id"] = stations["acro"].apply(lambda x: self.masterfile.find_entry(x, 'NAME', 'ID'))
            stations["id"] = stations["id"].replace(np.nan, 0.)
            stations["id"].astype(int)
        else:
            print("WARN: no masterfile provided, SPOTGINS IDs are forced to 0!")
            stations["id"] = 0

        sum0 = np.sum(stations["id"] == 0)
        if sum0 > 0:
            print(f"WARN: {sum0} sitelog-imported sites have a dummy SPOTGINS ID (00000) bc. not in the master file")

        # finalize the DataFrame
        stations.set_index("acro", inplace=True)

        ### actual set of the stations attribute
        self.stations = stations
        self.set_headerlines()

        return None

    def extract(self, date_dt, acro9):
        """
        Extract the receiver, antenna and eccentricities for a given station and date
        """
        return spotgins_utils.extract_from_stafile(self, date_dt, acro9)


#### ASSOCIATED FUNCTIONS ####

def metadata2stationfile_data(mda_inp):
    """
    Convert a rinexmod's MetaData Object to a
    StationFile's data dictionary and a
    list of StationFile lines.

    Parameters
    ----------
    mda_inp : MetaData
        The rinexmod's MetaData object to be converted.

    Returns
    -------
    dsite : dict
        A dictionary containing the station data.
    slin : list
        A list representing the individual line for the StationFile.stations DataFrame.
    """
    dsite = {mda_inp.site_id9: {}}
    d = dsite[mda_inp.site_id9]
    # slin is the individual line for the StationFile.stations DataFrame
    # columns=['id', 'X', 'sX', 'Y', 'sY', 'Z', 'sZ', 'VX', 'sVX', 'VY', 'sVY', 'VZ', 'sVZ', 'epoch', 'ref_frame']

    slin = [mda_inp.site_id9,  ### this column will be set as the index of the DataFrame below
            00000,  ### will be replaced by the right ID from the master file below
            float(mda_inp.misc_meta["X coordinate (m)"]), 0.,
            float(mda_inp.misc_meta["Y coordinate (m)"]), 0.,
            float(mda_inp.misc_meta["Z coordinate (m)"]), 0.,
            0., 0., 0., 0., 0., 0., dt(2010, 1, 1), 'i20i20']
    for i, ins in enumerate(mda_inp.instrus):

        d[i] = {}
        d[i]['start'] = ins['dates'][0]
        d[i]['end'] = ins['dates'][1]
        d[i]['rec'] = ins['receiver']['Receiver Type']
        d[i]['rec_sn'] = ins['receiver']['Serial Number']
        # specific filter REGEX for alignment
        ali_raw = ins['antenna']['Alignment from True N']
        if isinstance(ali_raw, str):
            ali_regex = re.search(r"^-?\d+(\.\d+)?$", ali_raw)
            align = ali_regex.groups(1) if ali_regex else "0"
        else:
            align = ali_raw
        ant = (ins['antenna']['Antenna Type'][:16],
               ins['antenna']['Antenna Radome Type'],
               ins['antenna']['Serial Number'],
               align)
        d[i]['ant'] = ant
        d[i]['ecc_u'] = [float(ins['antenna']['Marker->ARP Up Ecc. (m)']), 0.]
        d[i]['ecc_n'] = [float(ins['antenna']['Marker->ARP North Ecc(m)']), 0.]
        d[i]['ecc_e'] = [float(ins['antenna']['Marker->ARP East Ecc(m)']), 0.]
        d[i]['source'] = mda_inp.filename

    return dsite, slin
