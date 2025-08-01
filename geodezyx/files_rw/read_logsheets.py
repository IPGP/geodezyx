#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.files_rw contains reading functions to 
import GNSS logsheets/sitelogs.

it can be imported directly with:
from geodezyx import files_rw

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""


########## BEGIN IMPORT ##########
#### External modules
import copy
import datetime as dt
import glob
#### Import the logger
import logging
import os
import re

import dateutil
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

log = logging.getLogger('geodezyx')


##########  END IMPORT  ##########

#  _                     _               _                            _
# | |                   | |             | |                          | |
# | |     ___   __ _ ___| |__   ___  ___| |_ ___   _ __ ___  __ _  __| | ___ _ __
# | |    / _ \ / _` / __| '_ \ / _ \/ _ \ __/ __| | '__/ _ \/ _` |/ _` |/ _ \ '__|
# | |___| (_) | (_| \__ \ | | |  __/  __/ |_\__ \ | | |  __/ (_| | (_| |  __/ |
# |______\___/ \__, |___/_| |_|\___|\___|\__|___/ |_|  \___|\__,_|\__,_|\___|_|
#               __/ |
#              |___/
#



# ===== functions =====

def read_blocks_logsheet(input_file, block_id):
    """
    Read a Logsheet Block

    Parameters
    ----------
    input_file : str
        Logsheet path.
    block_id : int
        the ID of the wished block.
        1.   Site Identification of the GNSS Monument
        2.   Site Location Information 
        3.   GNSS Receiver Information
        4.   GNSS Antenna Information
        5.   Surveyed Local Ties
        6.   Frequency Standard
        7.   Collocation Information
        8.   Meteorological Instrumentation
        etc. etc. ...

    Returns
    -------
    ObjList : List of Block Objects
         List of Block Objects.
         Warning: if the Block contains only one sub-block
         Then the output is a 1-item List
    """
    proto_objects = [None,Site(),Location(),Receiver(),Antenna()]

    blkstr = str(block_id) + '.'
    blkstrnxt = str(block_id+1) + '.'
    blkstrX = str(block_id) + '.x'

    ObjList = []
    inblock = False
    

    for line in open(input_file,'rb'):
        
        ### manage unknown caracters
        line = str(line.decode('utf-8',errors='ignore'))


        if line == '\n':
            continue

        if line.startswith(blkstrnxt) or line.startswith(blkstrX):
            inblock = False
            break

        if line.startswith(blkstr) and (not 'Information' in line or block_id == 2):
            Obj = copy.copy(proto_objects[block_id])
            ObjList.append(Obj)
            inblock = True

        if inblock:
            if 'Approximate Position' in line:
                continue
            prop = line[4:30].strip().replace(' ', '_').replace('/','_').replace('(','').replace(')','').replace('.','').replace(',','_')
            if 'Marker->ARP Up Ecc' in prop:
                prop = 'Up_Ecc'
            elif 'Marker->ARP North Ecc' in prop:
                prop = 'North_Ecc'
            elif 'Marker->ARP East Ecc' in prop:
                prop = 'East_Ecc'
            elif 'Latitude (N is +)' in prop:
                prop = 'Latitude'
            elif 'Longitude (E is +)' in prop:
                prop = 'Longitude'
            elif 'Elevation (m,ellips.)' in prop:
                prop = 'Elevation'
            data = line[31:].strip()
            try:
                setattr(Obj,prop,float(data))
            except:
                setattr(Obj,prop,data)

    return ObjList

def mono_logsheet_read(logsheet_path,return_lists = False):
    """
    Read a sigle LogSheet File to Logsheet block objects

    Parameters
    ----------
    logsheet_path : str
        Logsheet path.
    return_lists : Bool, optional
        The default is False.
        if True, this legacy mode returns each Logsheet Block object in a list, 
        so each of them can be managed immediatly by
        write_station_info_from_datalists function.
        Useful only for such write_<...> functions.

    Returns
    -------
    period_lis : List of Period Block Object
        Period description in Period Block Object.
    sit : Site Block Object
        Site description in Site Block Object.
    loc : Location Block Object
        Location description in Site Block Object.
    """
    
    sit_lis = read_blocks_logsheet(logsheet_path,1)
    loc_lis = read_blocks_logsheet(logsheet_path,2)
    rec_lis = read_blocks_logsheet(logsheet_path,3)
    ant_lis = read_blocks_logsheet(logsheet_path,4)

    ant_install = [e.Date_Installed for e in ant_lis]
    rec_install = [e.Date_Installed for e in rec_lis]

    # merging all install date
    all_install_date = sorted(list(set(ant_install + rec_install )))

    # for each install date, find the Rec/Ant couple
    date_ant_rec_couple_lis = []
    for d in all_install_date:
        # exclude extremal date
        if d > dt.datetime(2098,1,1,tzinfo=dateutil.tz.tzutc()):
            continue

        potential_ant = []
        potential_rec = []

        for a in ant_lis:
            if a.Date_Installed <= d < a.Date_Removed:
                potential_ant.append(a)
        for r in rec_lis:
            if r.Date_Installed <= d < r.Date_Removed:
                potential_rec.append(r)

        site_code = sit_lis[0].Four_Character_ID
        if len(potential_ant) == 0:
            log.warning('%s: missing Antenna, %s, skip ...',site_code,d)
            continue
        if len(potential_rec) == 0:
            log.warning('%s: missing Receiver, %s, skip ...',site_code,d)
            continue
        if len(potential_ant) != 1:
            log.warning('s%: several Antennas found, %s',site_code,d)
        if len(potential_rec) != 1:
            log.warning('s%: several Receivers found, %s',site_code,d)

        date_ant_rec_couple_lis.append((d , potential_ant[0] , potential_rec[0]))

    if len(date_ant_rec_couple_lis) != len(set(date_ant_rec_couple_lis)):
        log.debug('bug 2')

    # construction of a period
    period_lis = []
    for i in range(len(date_ant_rec_couple_lis)):
        d1 = date_ant_rec_couple_lis[i][0]
        if i+1 == len(date_ant_rec_couple_lis):
            d2 = dt.datetime(2099,1,1,tzinfo=dateutil.tz.tzutc())
        else:
            d2 = date_ant_rec_couple_lis[i+1][0]
        a = date_ant_rec_couple_lis[i][1]
        r = date_ant_rec_couple_lis[i][2]

        period_lis.append((d1 , d2 , a , r))

    sit = sit_lis[0]
    loc = loc_lis[0]


    if return_lists:
        return [period_lis] , [sit] , [loc]
    else:
        return period_lis , sit , loc


def multi_logsheet_read(pathin,wildcardin='*log',return_dico=True,
                        output_mode='classic'):

    """
    Read multiple logsheets

    Parameters
    ----------
    pathin : list or str.
        If list: list of the Logsheet paths.
        If str: path of the directory containing the Logsheets.
    wildcardin : str, optional
        a wildcard describing the logsheet pattern.
        used only for pathin in a string path.
        The default is '*log'.
    return_dico : bool, optional
        If False, returns period_lis_lis , stat_lis , loc_lis
        the False mode is useful for station.info generation
        with write_station_info_from_datalists
        The default is True.
    output_mode : str, optional
        Defines the output mode if a dictionary is asked.
        (It is not relevant if list are asked)
        'classic': returns in the period_lis , stat , loc
        'pretty': returns the period_lis, stat, loc but with 
        period_lis as a DataFrame
        'legacy': returns in the dico period_lis_lis , stat_lis , loc_lis 
        (should not be used anymore)
        The default is 'classic'.
        
    Returns
    -------
    stations_dico : dict
        a dictionnary like stations_dico['STAT'] = (period_lis, station, loc).

    """
    
    if utils.is_iterable(pathin):
        logsheet_list = pathin
    else:
        fullpath = os.path.join(pathin,wildcardin)
        logsheet_list = sorted(glob.glob(fullpath))
        
        if not logsheet_list:
            log.error("no logsheets found, exiting ...")
            log.error(fullpath)
            return None
        
    log.info("logsheet_list: ")
    log.info(logsheet_list)

    period_lis_lis = []
    stat_lis       = []
    loc_lis        = []

    stations_dico = dict()

    for ls in logsheet_list:
        try:
            p,s,l = mono_logsheet_read(ls,return_lists=True)
        except:
            log.warning("%s skipped for unknown reason ...",ls)
            log.warning("       logsheet must be checked")
            continue

        period_lis_lis = period_lis_lis + p
        stat_lis       = stat_lis       + s
        loc_lis        = loc_lis        + l
        
        if output_mode == "classic":
            stations_dico[s[0].Four_Character_ID] = (p[0],s[0],l[0])
        elif output_mode == "pretty":
            p_DF = pd.DataFrame(p[0])
            p_DF.columns = ('start','end','ant','rec')
            stations_dico[s[0].Four_Character_ID] = (p_DF,s[0],l[0])
        elif output_mode == "legacy":
            stations_dico[s[0].Four_Character_ID] = (p,s,l)
        else:
            log.error('check output_mode value')
            raise Exception   


    if not return_dico:
        return period_lis_lis , stat_lis , loc_lis
    else:
        return stations_dico
    
# ============================= OBJECTS ======================================

class Event(object):
    def __init__(self):
        self.__Date_Installed            =  dt.datetime(1980,1,1)
        self.__Date_Removed              =  dt.datetime(2099,1,1)

    @property
    def Date_Installed(self):
        return self.__Date_Installed

    @Date_Installed.setter
    def Date_Installed(self,indate):
        if isinstance(indate,dt.datetime):
            self.__Date_Installed = indate
        elif 'CCYY-MM-DD' in indate:
            self.__Date_Installed = dt.datetime(1980,1,1,
                                                tzinfo=dateutil.tz.tzutc())
        elif 'XXXX' in indate:
            self.__Date_Installed = dt.datetime(1980,1,1,
                                                tzinfo=dateutil.tz.tzutc())
        else:
            self.__Date_Installed = dateutil.parser.parse(indate).replace(hour=0,
            minute=0, second=0, microsecond=0,tzinfo=dateutil.tz.tzutc())

    @property
    def Date_Removed(self):
        return self.__Date_Removed

    @Date_Removed.setter
    def Date_Removed(self,indate):
        if isinstance(indate,dt.datetime):
            self.__Date_Removed = indate
        elif 'CCYY-MM-DD' in indate:
            self.__Date_Removed = dt.datetime(2099,1,1,
                                              tzinfo=dateutil.tz.tzutc())
        elif 'XXXX' in indate:
            self.__Date_Removed = dt.datetime(2099,1,1,
                                              tzinfo=dateutil.tz.tzutc())
        else:
            self.__Date_Removed = dateutil.parser.parse(indate).replace(hour=0,
            minute=0, second=0, microsecond=0,tzinfo=dateutil.tz.tzutc())

class Receiver(Event):
    def __init__(self):
        Event.__init__(self)
        self.Receiver_Type             =  'XXXX'
        self.Satellite_System          =  'XXXX'
        self.Serial_Number             =  'XXXX'
        self.Firmware_Version          =  'XXXX'
        self.Elevation_Cutoff_Setting  =  'XXXX'
        self.Temperature_Stabiliz      =  'XXXX'
        self.Additional_Information    =  'XXXX'
#        self.Date_Installed            = 'XXXX'
#        self.Date_Removed              = 'XXXX'

    def __repr__(self):
        return "Receiver Type  :{},\nDate Installed:{},\nDate Removed  :{}".format(self.Receiver_Type,
                                 self.Date_Installed,
                                 self.Date_Removed)


    def FirmwareSmart(self):
        aaa = str(self.Firmware_Version).split()
        for a in aaa:
            try:
                return float(a)
            except:
                continue
        return 0.


class Antenna(Event):
    def __init__(self):
        Event.__init__(self)
        self.Antenna_Type             =  'XXXX'
        self.Serial_Number            =  'XXXX'
        self.Antenna_Reference_Point  =  'XXXX'
        self.Up_Ecc                   =  0.
        self.North_Ecc                =  0.
        self.East_Ecc                 =  0.
        self.Alignment_from_True_N    =  0.
        self.__Antenna_Radome_Type    =  'XXXX'
        self.Radome_Serial_Number     =  'XXXX'
        self.Antenna_Cable_Type       =  'XXXX'
        self.Antenna_Cable_Length     =  'XXXX'
        self.Additional_Information   =  'XXXX'

    @property
    def Antenna_Radome_Type(self):
        return self.__Antenna_Radome_Type

    @Antenna_Radome_Type.setter
    def Antenna_Radome_Type(self,inradome):
        if inradome == '':
            self.__Antenna_Radome_Type = 'NONE'
        else:
            self.__Antenna_Radome_Type = inradome[0:4].upper()

    def __repr__(self):
        return "Antenna Type  :{},\nDate Installed:{},\nDate Removed  :{}".format(self.Antenna_Type,self.Date_Installed,self.Date_Removed)

    def AntTypSmart(self):
        #Elimination of the Radome type
        if len(self.Antenna_Type.split()) > 1:
            return self.Antenna_Type.split()[0]
        else:
            return self.Antenna_Type

    def ARPSmart(self):
        #return 'DH'+self.Antenna_Reference_Point
        # see http://rses.anu.edu.au/geodynamics/gps/papers/gamit/apdx2_ps.pdf
        #return 'DHPAB'
        return 'DHARP'


class Site(object):
    def __init__(self):
        self.Site_Name               =  'XXXX'
        self.Four_Character_ID       =  'XXXX'
        self.Monument_Inscription    =  'XXXX'
        self.__IERS_DOMES_Number     =  'XXXX'
        self.CDP_Number              =  'XXXX'
        self.Monument_Description    =  'XXXX'
        self.Height_of_the_Monument  =  'XXXX'
        self.Monument_Foundation     =  'XXXX'
        self.Foundation_Depth        =  'XXXX'
        self.Marker_Description      =  'XXXX'
        self.Date_Installed          =  'XXXX'
        self.Geologic_Characteristic =  'XXXX'
        self.Bedrock_Type            =  'XXXX'
        self.Bedrock_Condition       =  'XXXX'
        self.Fracture_Spacing        =  'XXXX'
        self.Fault_zones_nearby      =  'XXXX'
        self.Distance_activity       =  'XXXX'
        self.Additional_Information  =  'XXXX'

    @property
    def IERS_DOMES_Number(self):
        return  self.__IERS_DOMES_Number
    @IERS_DOMES_Number.setter
    def IERS_DOMES_Number(self,indomes):
        #if indomes == '':
        #    self.__IERS_DOMES_Number = '99999M001'
        #else:
        #    self.__IERS_DOMES_Number = indomes

        if re.search('[0-9]{4}M[0-9]{3}', indomes.strip()):
            self.__IERS_DOMES_Number = indomes.strip()
        else:
            self.__IERS_DOMES_Number = '99999M001'

    def __repr__(self):
        return "Site Name   :{},\nCode        :{},\nInstallation:{}".format(self.Site_Name,
                                 self.Four_Character_ID,
                                 self.Date_Installed)


class Location(object):
    def __init__(self):
        self.City_or_Town            = 'XXXX'
        self.State_or_Province       = 'XXXX'
        self.Country                 = 'XXXX'
        self.Tectonic_Plate          = 'XXXX'
        self.X_coordinate_m          = 0.
        self.Y_coordinate_m          = 0.
        self.Z_coordinate_m          = 0.
        self.Latitude                = 0.
        self.Longitude               = 0.
        self.Elevation_m_ellips      = 0.
        self.Additional_Information  = 'XXXX'

        self.X_velocity              = 0.
        self.Y_velocity              = 0.
        self.Z_velocity              = 0.

        self.X_coordinate_sigma      = 0.008
        self.Y_coordinate_sigma      = 0.008
        self.Z_coordinate_sigma      = 0.008

        self.X_velocity_sigma        = 0.008
        self.Y_velocity_sigma        = 0.008
        self.Z_velocity_sigma        = 0.008

        self.Reference_epoch         = dt.datetime(2005,1,1)

    def __repr__(self):
        return "X:{},Y:{},Z:{},\nRef. epoc.:{}".format(self.X_coordinate_m,
                                                     self.Y_coordinate_m,
                                                     self.Z_coordinate_m,
                                                     self.Reference_epoch)

    def export_as_string(self):
        Str_list = []
        Str_list.append('City or Town             : {:}'.format(self.City_or_Town))
        Str_list.append('State or Province        : {:}'.format(self.State_or_Province))
        Str_list.append('Country                  : {:}'.format(self.Country))
        Str_list.append('Tectonic Plate           : {:}'.format(self.Tectonic_Plate))
        Str_list.append('Approximate Position (ITRF)')

        X = self.X_coordinate_m
        Y = self.Y_coordinate_m
        Z = self.Z_coordinate_m

        lat , lon , h = conv.xyz2geo(X, Y, Z)
        lat_deg , lat_min , lat_sec = conv.deg2deg_dec2dms(lat)
        lon_deg , lon_min , lon_sec = conv.deg2deg_dec2dms(lon)

        Str_list.append('X coordinate (m)       : {:}'.format(X))
        Str_list.append('Y coordinate (m)       : {:}'.format(Y))
        Str_list.append('Z coordinate (m)       : {:}'.format(Z))

        Str_list.append('Latitude (N is +)      : {:+3d}{:2d}{:5.2d}'.format((lat_deg,lat_min,lat_sec)))
        Str_list.append('Longitude (E is +)     : {:+3d}{:2d}{:5.2d}'.format((lat_deg,lat_min,lat_sec)))
        Str_list.append('Elevation (m,ellips.)  : {:7.1f}'.format(h))
        Str_list.append('Additional Information   : N/A')

        out_str = "\n".join(Str_list)

        return out_str
