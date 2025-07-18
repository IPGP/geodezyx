#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from datetime import datetime as dt


class OffsetsFile:

    def __init__(self):

        self.write_amplitude = False

        self.header = ['# Offset file v2.0\n'
                       '# SOLUTION        : SPOTGINS\n'
                       '# SITE            : 9-character site ID\n'
                       '# YYYYMMDDHRMN    : Year, month, day, hour, minute\n'
                       '# TYPE            : One among the following list:\n'
                       '#                   - Antenna_and_Radome_Type_Changed\n'
                       '#                   - Antenna_Code_Changed\n'
                       '#                   - Antenna_Type_Changed\n'
                       '#                   - Elevation_Cutoff_Changed\n'
                       '#                   - Equipment_Site_Change\n'
                       '#                   - Logfile_Changed\n'
                       '#                   - Radome_Type_Changed\n'
                       '#                   - Receiver_Make_and_Model_Changed\n'
                       '#                   - Receiver_Make_Changed\n'
                       '#                   - Receiver_Model_Changed\n'
                       '#                   - Receiver_Setting_Corrected\n'
                       '#                   - Temporary_Step_Product_Change\n'
                       '#                   - Unknown_Event\n'
                       '#                   - Volcanic_Eruption\n'
                       '#                   - Earthquake_Event\n'
                       '#                   - Other_Event\n'
                       '#                   Old deprecated types:\n'
                       '#                   - CHANGE (equipment)\n'
                       '#                   - EQ (earthquake)\n'
                       '#                   - UNK (unknown)\n'
                       '# OFFSET          : Is an offset effectively detected in the station position time series? [Optional, Y/N/U]\n'
                       '#                   - Y : yes\n'
                       '#                   - N : no\n'
                       '#                   - U : unknown\n'
                       '# DN/SN           : Estimated north offset and 1-sigma standard deviation in mm [Optional, -9999999 by default]\n'
                       '# DE/SE           : Estimated east offset and 1-sigma standard deviation in mm [Optional, -9999999 by default]\n'
                       '# DU/SU           : Estimated height offset and 1-sigma standard deviation in mm [Optional, -9999999 by default]\n'
                       '# COMMENT         : Free format comment [Optional]\n'
                       '#\n'
                       '#____SITE YYYYMMDDHRMN ___________________________TYPE OFFSET ______DN ______SN ______DE ______SE ______DU ______SU COMMENT\n']

        self.stations = []
        self.data = {}

    def add_offset(self, acro, date_offset, type, description,
                   dn=np.nan, sn=np.nan,
                   de=np.nan, se=np.nan,
                   du=np.nan, su=np.nan):

        if acro not in self.stations:
            self.stations.append(acro)
            self.data[acro] = {}

        date_jd = round(dt2jd(date_offset), 4)
        date_str = date_offset.strftime('%Y%m%d%H%M')

        self.data[acro][date_jd] = {'datestr': date_str, 'type': type,
                                    'dn': dn, 'sn': sn, 'de': de,
                                    'se': se, 'du': du, 'su': su,
                                    'description': description}

    def remove_offset(self, acro, date_offset):

        date_jd = round(dt2jd(date_offset), 4)

        if date_jd not in self.data[acro].keys():
            print('FATAL: No offset found for ' + acro + ' at date %.4f' % date_jd)
        else:
            del self.data[acro][date_jd]

    def modify_offset(self, acro, date_offset,
                      newdate=None, newtype=None, newdesc=None,
                      newdn=None, newsn=None,
                      newde=None, newse=None,
                      newdu=None, newsu=None):

        date_jd = round(dt2jd(date_offset), 4)

        if date_jd not in self.data[acro].keys():
            print('FATAL: No offset found for ' + acro + ' at date %.4f' % date_jd)
        else:
            if newtype is not None:
                self.data[acro][date_jd]['type'] = newtype
            if newdesc is not None:
                self.data[acro][date_jd]['description'] = newdesc
            if newdn is not None:
                self.data[acro][date_jd]['dn'] = newdn
            if newsn is not None:
                self.data[acro][date_jd]['sn'] = newsn
            if newdn is not None:
                self.data[acro][date_jd]['de'] = newde
            if newse is not None:
                self.data[acro][date_jd]['se'] = newse
            if newdu is not None:
                self.data[acro][date_jd]['du'] = newdu
            if newsu is not None:
                self.data[acro][date_jd]['su'] = newsu

            if newdate is not None:
                print('Modification with newdate ' + str(newdate) + ' at date ' + str(date_offset))
                newdate_jd = dt2jd(newdate)
                # self.add_offset()
                self.data[acro][newdate_jd] = self.data[acro][date_jd]
                self.data[acro][newdate_jd]['datestr'] = newdate.strftime('%Y%m%d%H%M')

                self.remove_offset(acro, date_offset)

    def write(self, outfile='offsets_file.dat', write_amplitude=False):

        dataout = open(outfile, 'w')

        if not self.write_amplitude and write_amplitude:
            self.write_amplitude = True

        # header
        if not self.write_amplitude:
            self.header[-1] = '#____SITE DDDDD.DDDD YYYYMMDDHRMN __TYPE DESCRIPTION\n'

        for hl in self.header:
            write = True
            if 'Estimated' in hl and not self.write_amplitude:
                write = False

            if write:
                dataout.write(hl)

        # offsets data
        for acro in self.stations:
            for date_offset_jd in sorted(self.data[acro].keys()):
                offset_data = self.data[acro][date_offset_jd]
                date_offset_str = offset_data['datestr']
                type = offset_data['type']
                desc = offset_data['description']

                if not self.write_amplitude:
                    dataout.write(
                        acro + ' ' + '%.4f' % date_offset_jd + ' ' + date_offset_str + ' ' + type + ' ' + desc + '\n')
                else:

                    dn = ('%.1f' % offset_data['dn']).rjust(8)
                    # print (offset_data['sn'])

                    if str(dn).strip() == 'nan':
                        dn = '-9999999'

                    if offset_data['sn'] == '-9999999' or str(offset_data['sn'] == 'nan'):
                        sn = '-9999999'
                    else:
                        # print (offset_data['sn'])
                        sn = ('%.1f' % offset_data['sn']).rjust(8)

                    de = ('%.1f' % offset_data['de']).rjust(8)
                    if str(de).strip() == 'nan':
                        de = '-9999999'

                    if (offset_data['se'] == '-9999999') or str(offset_data['se'] == 'nan'):
                        se = '-9999999'
                    else:
                        se = ('%.1f' % offset_data['se']).rjust(8)

                    du = ('%.1f' % offset_data['du']).rjust(8)
                    if str(du).strip() == 'nan':
                        du = '-9999999'
                    if offset_data['su'] == '-9999999' or str(offset_data['su'] == 'nan'):
                        su = '-9999999'
                    else:
                        su = ('%.1f' % offset_data['su']).rjust(8)

                    dataout.write(acro + ' ' + '%.4f' % date_offset_jd + ' ' + date_offset_str + ' ' + type.ljust(
                        6) + ' ' + dn + ' ' + sn + ' ' + de + ' ' + se + ' ' + du + ' ' + su + ' ' + desc + '\n')

        dataout.close()

    @classmethod
    def from_file(cls, offsetfilename):

        offsetfile = OffsetsFile()

        offdata = open(offsetfilename)
        values = False
        for line in offdata:

            if "#____SITE" in line:
                # header line
                if "______DN" in line:
                    values = True
                    break
        offdata.close()
        offdata = open(offsetfilename)
        for line in offdata:
            if line[0:1] != '#':
                cols = line.strip().split()

                # print (cols)
                # print (len(cols))

                if values:
                    (acro, date_jd, date_str, type, dn, sn, de, se, du, su) = line[0:94].split()
                    desc = line[95:].strip()

                    dn = float(dn)
                    de = float(de)
                    du = float(du)
                    if sn == 'nan':
                        sn = np.nan
                    if se == 'nan':
                        se = np.nan
                    if su == 'nan':
                        su = np.nan
                else:
                    (acro, date_jd, date_str, type) = line[0:40].split()
                    desc = line[41:].strip()
                    dn = '-9999999'
                    sn = '-9999999'
                    de = '-9999999'
                    se = '-9999999'
                    du = '-9999999'
                    su = '-9999999'

                date_dt = dt.strptime(date_str, '%Y%m%d%H%M')

                type = type.ljust(6)
                # print (sn)
                offsetfile.add_offset(acro, date_dt, type, desc, dn, sn, de, se, du, su)
        offdata.close()

        return offsetfile

    @classmethod
    def from_station_file(cls, stationfile):  # ,masterfile):

        # sta = StationFile.from_file(stationfile, masterfile)

        offsetfile = OffsetsFile()

        # data
        codes = {0: 'Other site log entry (firmware upgrade, cable change, ...)',
                 1: 'Receiver change',
                 2: 'Antenna change',
                 3: 'Receiver & Antenna change',
                 4: 'Radome change',
                 5: 'Receiver & Radome change',
                 6: 'Antenna & Radome change',
                 7: 'Receiver, Antenna & Radome change'}

        for acro in stationfile.data.keys():
            rec = ''
            rec_sn = ''
            ant = ''
            ant_sn = ''
            dome = ''
            ant_ori = 0.
            for neq in stationfile.data[acro].keys():
                eq = stationfile.data[acro][neq]
                rec_tmp = eq['rec']
                rec_sn_tmp = eq['rec_sn']
                (ant_tmp, dome_tmp, ant_sn_tmp, ant_ori_tmp) = eq['ant']
                ant_ori_tmp = float(ant_ori_tmp)

                ## Antenna orientation must be Numbers (not NaN)
                ant_ori_are_num = not (np.isnan(ant_ori_tmp) or np.isnan(ant_ori))

                if neq != 0:
                    code = 0
                    if (rec != rec_tmp) or (rec_sn != rec_sn_tmp):
                        code += 1
                    if ((ant != ant_tmp) or (ant_sn != ant_sn_tmp) or
                            (ant_ori_are_num and ant_ori != ant_ori_tmp)):
                        code += 2
                    if dome != dome_tmp:
                        code += 4

                    # if code > 0:
                    # change
                    desc = codes[code]
                    chg_date = eq['start']

                    offsetfile.add_offset(acro=acro, date_offset=chg_date,
                                          type='CHANGE', description=desc)

                rec = rec_tmp
                rec_sn = rec_sn_tmp
                ant = ant_tmp
                dome = dome_tmp
                ant_sn = ant_sn_tmp
                ant_ori = ant_ori_tmp

        return offsetfile


def dt2jd(datedt):
    return datedt.toordinal() + datedt.hour / 24. + datedt.minute / 1440. + datedt.second / 86400. + 1721424.5 - 2400000.5
