#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 16:07:11 2019

@author: psakicki, mansur, chaiyap

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics

Copyright (C) 2019 Pierre Sakic (GFZ, pierre.sakic@gfz-postdam.de)
GitHub repository :
https://github.com/PierreS1/GeodeZYX-Toolbox-Lite

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""


import datetime as dt
import numpy as np
import dateutil.parser
import re
import glob
import scipy
from scipy.interpolate import interp1d
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.figure import Figure
import csv
import pandas as pd
from collections import Counter
import copy
import warnings
import time
from random import randrange
import operator
from natsort import natsorted, ns
import tabulate
from matplotlib.widgets import MultiCursor
import matplotlib
import linecache
import math
import softs_runner

################################################
import geodezyx.conv.time as convt
import geodezyx.utils.utils as utils


#import genefun as genefun
#import geo_files_converter_lib as gfc
#import geo_trop as gtro
#import geodetik as geok     maybe we should leave as comment until define the proper paths


def read_bull_B(path):

    if not genefun.is_iterable(path):
        path = [path]

    path = sorted(path)

    DFstk = []

    for path_solo in path:
        S = genefun.extract_text_between_elements(path_solo,"1 - DAILY FINAL VALUES" ,
                                         "2 - DAILY FINAL VALUES" )

        L = S.replace('\t','\n').split("\n")

        L2 = []
        for e in L:
            if len(e) > 0:
                if e[0] !=  " ":
                    L2.append(e)
        L3 = []
        for e in L2:
            L4 = []
            for ee in e.split():
                L4.append(float(ee))
            L3.append(L4)

        DF = pd.DataFrame(np.vstack(L3))
        DFstk.append(DF)

    DFout = pd.concat(DFstk)
    DFout.columns = ["year","month","day","MJD","x","y","UT1-UTC","dX","dY",
                  "x err","y err","UT1 err","X err","Y err"]

    return DFout





def read_clk(file_path_in, returns_pandas = True, interval=None):
    """
    General description

    Parameters
    ----------
    file_path_in :  str
        Path of the file in the local machine.

    returns_pandas :  bool
        Define if pandas function will be used to put the data on tables

    interval :  int
        Defines which interval should be used in the data tables. The interval is always in minutes unit.

    Returns
    -------
    out1 : float or int oandas table
        Returns a panda table format with the data extracted from the file.

    out2 :  list
       In case of the pandas function are not in use, the return will be a list with the data extract from the file.


    """

    file = open(file_path_in,'r')
    fil = file.readlines()
    file.close()


    types_m, name_m, year_m, month_m, day_m, hr_m, minute_m, seconds_m, epoch_m, offset_m, rms_m = [],[],[],[],[],[],[],[],[],[],[]

    Clk_read  = []

        #####IGNORES HEADER

    le = 0
    i = 0
    count = 0
    for le in range(len(fil)):
     linhaatual = linecache.getline(file_path_in, le)
     header_end = (linhaatual[60:73])
     count +=1
     if header_end =="END OF HEADER":
        i = count
#############################################################
    while i <= len(fil):
          linhaatual = linecache.getline(file_path_in, i)
          types = (linhaatual[0:2]) # Column refers to AS or AR
          name = (linhaatual[3:7])# Name of the station or satellite
          year = int((linhaatual[8:12]))
          month = int((linhaatual[13:15]))
          day = int((linhaatual[16:18]))
          hr= int((linhaatual[19:22]))# hour
          minute = int((linhaatual[22:25]))
          seconds = float((linhaatual[25:33]))
          offset = float(((linhaatual[40:59])))# Clock offset
          ##### check if there is a value for thr rms
          if (linhaatual[65:70]) == '     ' or len(linhaatual[65:70])==0:
              rms = 0
          else:
              rms = float(linhaatual[61:79])
          ############
          epoch = dt.datetime(year,month,day,hr,minute,int(seconds),int(((int(seconds)-seconds)*10**-6))) # epoch as date.time function

        ############ Data in a defined interval
          if interval:
              if (float(minute)%interval) == 0 and float(seconds) == 0:

                  if returns_pandas:
                            Clk_dados = [types,name,year,month,day,hr,minute,seconds,epoch,offset,rms]
                            Clk_read.append(Clk_dados)
                  else:
                            types_m.append(types)
                            name_m.append(name)
                            year_m.append(year)
                            month_m.append(month)
                            day_m.append(day)
                            hr_m.append(hr)
                            minute_m.append(minute)
                            seconds_m.append(seconds)
                            epoch_m.append(epoch)
                            offset_m.append(offset)
                            rms_m.append(rms)
         #########################################################
        ########################################################### standart file epochs
          else:
              if returns_pandas:
                            Clk_dados = [types,name,year,month,day,hr,minute,seconds,epoch,offset,rms]
                            Clk_read.append(Clk_dados)
              else:
                            types_m.append(types)
                            name_m.append(name)
                            year_m.append(year)
                            month_m.append(month)
                            day_m.append(day)
                            hr_m.append(hr)
                            minute_m.append(minute)
                            seconds_m.append(seconds)
                            epoch_m.append(epoch)
                            offset_m.append(offset)
                            rms_m.append(rms)


          i +=1

  ################################################################################################
   ############################################ put on pandas table format
    if returns_pandas:
     Clk_readed = pd.DataFrame(Clk_read, columns=['type','name','year','month','day','h','minutes','seconds','epoch','offset','rms'])
     Clk_readed['epoch'] =  pd.to_datetime(Clk_readed[[ 'year' ,'month' ,'day','h','minutes','seconds']])
     Clk_readed.path = file_path_in

     return Clk_readed

          ###############################
    else:
           print(" ...")
           return  types_m, name_m, year_m, month_m, day_m, hr_m, minute_m, seconds_m, offset_m, rms_m




############################# END FUNCTION READ CLK 30







def read_epos_sta_coords_mono(filein,return_df=True):
    """
    read an EPOS's YYYY_DDD_XX_sta_coordinates coordinates files
    and return a list of Points objects

    ... TBC ...
    """
    F = open(filein)

    Points_list_stk = []
    Lines_4_DF_stk = []

    for l in F:
        fields = l.split()
        if l[0] != " ":
            continue
        if "SITE" in fields[0]:
            namestat = fields[8]
            numstat  = int(fields[2])
            tecto_plate = fields[4]
            MJD_ref  = int(fields[5])
            MJD_strt = int(fields[6])
            MJD_end  = int(fields[7])
            MJD_mid = np.mean([MJD_strt , MJD_end])
            T = geok.MJD2dt(MJD_mid)                #### SHOULD BE MODIFIED THE GEOK FUNCTION DOES NO LONGER EXIST

        if "POS_VEL:XYZ" in fields[0]:
            X  = float(fields[4])
            Y  = float(fields[5])
            Z  = float(fields[6])
            Vx = float(fields[7])
            Vy = float(fields[8])
            Vz = float(fields[9])

        if "SIG_PV_XYZ" in fields[0]:
            sX  = float(fields[4].replace("D","E"))
            sY  = float(fields[5].replace("D","E"))
            sZ  = float(fields[6].replace("D","E"))
            sVx = float(fields[7])
            sVy = float(fields[8])
            sVz = float(fields[9])

            #### Last useful line for the point, store it
            point = Point(X,Y,Z,T,"XYZ",sX,sY,sZ,name=namestat)
            point.anex["Vx"] = sVx
            point.anex["Vy"] = sVy
            point.anex["Vz"] = sVz
            Points_list_stk.append(point)

            #### And store for the DataFrame
            tup_4_DF = (namestat,numstat,tecto_plate,
                        MJD_ref,MJD_strt,MJD_end,
                        X,Y,Z,sX,sY,sZ,
                        Vx,Vy,Vz,sVx,sVy,sVz)



            Lines_4_DF_stk.append(tup_4_DF)


        columns = ("site","site_num","tecto_plate",
                   "MJD_ref","MJD_start","MJD_end",
                   "x","y","z","sx","sy","sz",
                   "Vx","Vy","Vz","sVx","sVy","sVz")

        DFout = pd.DataFrame(Lines_4_DF_stk,
                     columns=columns)


    if return_df:
        return DFout
    else:
        return Points_list_stk



def read_epos_sta_coords_multi(filein_list,return_dict = True):

    filein_list  = sorted(filein_list)
    Points_list  = []
    statname_stk = []

    for fil in filein_list:
        Points_daily_list = read_epos_sta_coords_mono(fil,return_df=False)
        Points_list   = Points_list + Points_daily_list
        statname_stk  = statname_stk + [e.name for e in Points_daily_list]

    statname_uniq = sorted(list(set(statname_stk)))

    ts_dict = dict()

    for point in Points_list:
        if not point.name in ts_dict.keys():
            ts_dict[point.name] = TimeSeriePoint(stat=point.name)
        ts_dict[point.name].add_point(point)

    if return_dict:
        return ts_dict
    else:
        ts_list = []
        for k , val in ts_dict.items():
            ts_list.append(val)
        return ts_list


def read_epos_sta_kinematics(filein):
    """
    read an EPOS kinematic solutions
    """

    F = open(filein)
    Lines_4_DF_stk = []
    for l in F:
        fields = l.split()
        if l[0] != "K" and l[0] != "U" and l[0] != "X":
            continue
        if l[0] == "K" or l[0] == "U" or l[0] == "X":
            namstat = fields[2]
            numstat = int(fields[1])
            MJD_epo = float(fields[3])
            numobs = int(fields[4])

            X = float(fields[6])
            Y = float(fields[7])
            Z = float(fields[8])
            sX = float(fields[10])
            sY = float(fields[11])
            sZ = float(fields[12])

            N = float(fields[14])
            E = float(fields[15])
            U = float(fields[16])
            sN = float(fields[18])
            sE = float(fields[19])
            sU = float(fields[20])

            tup_4_df = (namstat,numstat,MJD_epo,numobs,X,Y,Z,sX,sY,sZ,
                        N,E,U,sN,sE,sU)
            Lines_4_DF_stk.append(tup_4_df)

    columns = ("site","site_num",
                   "MJD_epo","numobs",
                   "x","y","z","sx","sy","sz",
                   "N","E","U","sN","sE","sU")

    DFout = pd.DataFrame(Lines_4_DF_stk,
                     columns=columns)
    return DFout




def read_erp(caminho_arq,ac):
    """
    General description
    Units: ('MJD','X-P (arcsec)', 'Y-P (arcsec)', 'UT1UTC (E-7S)','LOD (E-7S/D)','S-X (E-6" arcsec)','S-Y (E-6" arcsec)',
    'S-UT (E-7S)','S-LD (E-7S/D)','NR (E-6" arcsec)', 'NF (E-6" arcsec)', 'NT (E-6" arcsec)',
    'X-RT (arcsec/D)','Y-RT (arcsec/D)','S-XR (E-6" arcsec/D)','S-YR (E-6" arcsec/D)', 'C-XY', 'C-XT',
    'C-YT', 'DPSI', 'DEPS','S-DP','S-DE')

    Parameters
    ----------
    file_path_in :  str
        Path of the file in the local machine.

    which AC :  str
        The analisys center that will be used


    Returns
    -------
    out1 :  pandas table
        Returns a panda table format with the data extracted from the file.

    Obs: this function need to be improved to read also the values of UTC and etc. So far reads only the Pole rates. 

    """

    #### FIND DELIVERY DATE
    name = os.path.basename(caminho_arq)

    if len(name) == 12:
        dt_delivery = convt.sp3name2dt(caminho_arq)
    elif len(name) == 38:
        dt_delivery = convt.sp3name_v3_2dt(caminho_arq)
    else:
        dt_delivery = convt.posix2dt(0)


    le = open(caminho_arq, 'r')
    letudo = le.readlines()
    le.close()
    tamanho = len(letudo) 

    

    numeros = ['0','1','2','3','4','5','6','7','8','9']


    ERP=[]


    if caminho_arq[-3:] in ('snx','ssc'):
        file = open(caminho_arq)
        Lines =  file.readlines()
        XPO_stk  = []
        XPO_std_stk = []
        YPO_stk  = []
        YPO_std_stk = []
        LOD_stk  = []
        LOD_std_stk = []
        MJD_stk = []
        marker = False

        for i in range(len(Lines)):

            if len(Lines[i].strip()) == 0:
                continue
            else:

                if Lines[i].split()[0] == '+SOLUTION/ESTIMATE':
                    marker = True

                if Lines[i].split()[0] == '-SOLUTION/ESTIMATE':
                    marker = False

                if utils.contains_word(Lines[i],'XPO') and marker:
                    Doy = (Lines[i][30:33])
                    Year = (Lines[i][27:29])
                    Pref_year = '20'
                    Year = int(Pref_year+Year)
                    Date = convt.doy2dt(Year,Doy)
                    XPO = float(Lines[i][47:68])*(10**-3)
                    XPO_std = float(Lines[i][69:80])*(10**-3)
                    XPO_stk.append(XPO)
                    XPO_std_stk.append(XPO_std)
                    MJD_stk.append(convt.jd_to_mjd(convt.date_to_jd(Date.year,Date.month,Date.day)))

                if utils.contains_word(Lines[i],'YPO') and marker:
                    Doy = (Lines[i][30:33])
                    Year = str(Lines[i][27:29])
                    Pref_year = '20'
                    Year = int(Pref_year+Year)
                    Date = convt.doy2dt(Year,Doy)
                    YPO = float(Lines[i][47:68])*(10**-3)
                    YPO_std = float(Lines[i][69:80])*(10**-3)
                    YPO_stk.append(YPO)
                    YPO_std_stk.append(YPO_std)
                    MJD_stk.append(convt.jd_to_mjd(convt.date_to_jd(Date.year,Date.month,Date.day)))

                if utils.contains_word(Lines[i],'LOD') and marker:
                    Doy = (Lines[i][30:33])
                    Year = str(Lines[i][27:29])
                    Pref_year = '20'
                    Year = int(Pref_year+Year)
                    Date = convt.doy2dt(Year,Doy)
                    LOD = float(Lines[i][47:68])*(10**+4)
                    LOD_std = float(Lines[i][69:80])*(10**+4)
                    LOD_stk.append(LOD)
                    LOD_std_stk.append(LOD_std)
                    MJD_stk.append(convt.jd_to_mjd(convt.date_to_jd(Date.year,Date.month,Date.day)))

        MJD = list(set(MJD_stk))
        if len(LOD_stk) == 0:
                LOD_stk = ['0']*len(MJD)
                LOD_std_stk = ['0']*len(MJD)

        for i in range(len(MJD)):

            ERP_data = [ac, MJD[i], XPO_stk[i], YPO_stk[i], 0, LOD_stk[i], XPO_std_stk[i], YPO_std_stk[i],
                     0, LOD_std_stk[i], 0, 0, 0, 0, 0, 0, 0, dt_delivery]
            ERP.append(ERP_data)



    if ac in ('COD','cod','com', 'cof', 'grg', 'mit', 'sio'):
        for i in range(tamanho+1):
            linhaatual = linecache.getline(caminho_arq, i)
            if linhaatual[0:1] in numeros:
                ERP_data = linhaatual.split()
                for j in range(len(ERP_data)):
                    ERP_data[j] = float(ERP_data[j])
                ERP_data.insert(0,ac)
                ERP_data[2] =  ERP_data[2]*(10**-6)
                ERP_data[3] =  ERP_data[3]*(10**-6)
                ERP_data[13] = ERP_data[13]*(10**-6)
                ERP_data[14] = ERP_data[14]*(10**-6)
                del ERP_data[17:]
                ERP_data.append(dt_delivery)

                ERP.append(ERP_data)



    if ac in ('wum','grg','esa', 'mit', 'ngs', 'sio'):
        for i in range(tamanho+1):
            linhaatual = linecache.getline(caminho_arq, i)
            if linhaatual[0:1] in numeros:
                ERP_data = linhaatual.split()
                for j in range(len(ERP_data)):
                    ERP_data[j] = float(ERP_data[j])
                ERP_data.insert(0,ac)
                ERP_data[2] =  ERP_data[2]*(10**-6)
                ERP_data[3] =  ERP_data[3]*(10**-6)
                ERP_data[13] = ERP_data[13]*(10**-6)
                ERP_data[14] = ERP_data[14]*(10**-6)
                ERP_data.append(dt_delivery)

                ERP.append(ERP_data)

#
    if ac in ('gbm', 'gfz'):
        for i in range(tamanho+1):
            linhaatual = linecache.getline(caminho_arq, i)
            if linhaatual[0:1] in numeros:
                ERP_data = linhaatual.split()
                for j in range(len(ERP_data)):
                    ERP_data[j] = float(ERP_data[j])
                ERP_data.insert(0,ac)
                ERP_data[2] =  ERP_data[2]*(10**-6)
                ERP_data[3] =  ERP_data[3]*(10**-6)
                ERP_data[13] = ERP_data[13]*(10**-6)
                ERP_data[14] = ERP_data[14]*(10**-6)
                ERP_data.append(dt_delivery)

                ERP.append(ERP_data)


    header = []
    if ac in ('emr'):
        for i in range(tamanho+1):
            linhaatual = linecache.getline(caminho_arq, i)
            if linhaatual == 'EOP  SOLUTION':
                del ERP_data[:]
                header = ['EOP  SOLUTION']
            if linhaatual[0:1] in numeros and 'EOP  SOLUTION' in header:
                ERP_data = linhaatual.split()
                for j in range(len(ERP_data)):
                    ERP_data[j] = float(ERP_data[j])
                ERP_data.insert(0,ac)
                ERP_data[2] =  ERP_data[2]*(10**-6)
                ERP_data[3] =  ERP_data[3]*(10**-6)
                ERP_data[13] = ERP_data[13]*(10**-6)
                ERP_data[14] = ERP_data[14]*(10**-6)
                del ERP_data[17:]
                ERP_data.append(dt_delivery)

                ERP.append(ERP_data)

    Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD',
                                         'NR', 'NF', 'NT',
                                         'X-RT','Y-RT','S-XR','S-YR',
                                         'Delivered_date'])
    return Erp_end







def read_sp3(file_path_in,returns_pandas = True, name = '',
             epoch_as_pd_index = False):
    """
    Read a SP3 file (GNSS Orbits standard file) and return X,Y,Z coordinates
    for each satellite and for each epoch

    Parameters
    ----------
    file_path_in : str
        path of the SP3 file

    returns_pandas : bool
        if True, return a Pandas DataFrame.
        if False, return separated lists.

    name : str
        a manual name for the file

    epoch_as_pd_index : bool
        if True, the index of the output dataframe contains
        if False, it contains generic integer indexs


    Returns
    -------
    df : Pandas DataFrame
        if returns_pandas == True

    epoch_stk ,  Xstk , Ystk , Zstk , Clkstk : lists
        if returns_pandas == False

    """

    AC_name =  os.path.basename(file_path_in)[:3]

    fil = open(file_path_in)

    header = True

    epoch_stk = []
    Xstk , Ystk , Zstk , Clkstk = [],[],[],[]

    data_stk  = []

    for l in fil:
        if l[0] == '*':
            header = False

        if header:
            continue
        if 'EOF' in l:
            continue

        if l[0] == '*':
            epoc   = convt.tup_or_lis2dt(l[1:].strip().split())
        else:
            sat_nat = l[1:2].strip()
            sat_sv  = int(l[2:4].strip())
            sat_sat = l[1:4].strip()

            X   = float(l[4:18])
            Y   = float(l[18:32])
            Z   = float(l[32:46])
            Clk = float(l[46:60])

            if returns_pandas:
                line_data = [epoc,sat_sat,sat_nat,sat_sv,X,Y,Z,Clk,AC_name]
                data_stk.append(line_data)
            else:
                epoch_stk.append(epoc)
                Xstk.append(X)
                Ystk.append(Y)
                Zstk.append(Z)
                Clkstk.append(Clk)


    AC_name_stk = [AC_name] * len(Xstk)

    if returns_pandas:
        df = pd.DataFrame(data_stk, columns=['epoch','sat', 'const', 'sv',
                                             'x','y','z','clk','AC'])
        if epoch_as_pd_index:
            df.set_index('epoch',inplace=True)
        df.filename = os.path.basename(file_path_in)
        df.path = file_path_in

        if name != '':
            df.name = name
        else:
            df.name = os.path.basename(file_path_in)

        return df
    else:
        print("INFO : return list, very beta : no Sat. Vehicule Number info ...")
        return  epoch_stk ,  Xstk , Ystk , Zstk , Clkstk , AC_name_stk




def read_sp3_header(sp3_path):
    """
    Read a SP3 file header and return a Pandas DataFrame
    with sat. PRNs and sigmas contained in the header

    Parameters
    ----------
    sp3_path : str
        path of the SP3 file


    Returns
    -------
    df : Pandas DataFrame
        2 columns "sat", "sigma"

    Note
    -------
    More infos about the sigma
    http://acc.igs.org/orbacc.txt
    """


    F = open(sp3_path)
    ac_name = os.path.basename(sp3_path)[:3]


    Lines = F.readlines()

    Sat_prn_list = []
    Sat_sig_list = []

    for il , l in enumerate(Lines):
        if il == 1:
            date = convt.MJD2dt(int(l.split()[4]))
        if l[:2] == "+ ":
            Sat_prn_list.append(l)
        if l[:2] == "++":
            Sat_sig_list.append(l)
        if l[0] == "*":
            break

    ### PRN part
    Sat_prn_list_clean = []
    for prn_line in Sat_prn_list:
        prn_line_splited = prn_line.split()
        prn_line_splited = [e for e in prn_line_splited if not "+" in e]
        prn_line_splited = [e for e in prn_line_splited if not  e == "0"]
        Sat_prn_list_clean = Sat_prn_list_clean + prn_line_splited

    sat_nbr = int(Sat_prn_list_clean[0])

    Sat_prn_list_clean = Sat_prn_list_clean[1:]

    Sat_prn_string = "".join(Sat_prn_list_clean)

    Sat_prn_list_final = []
    for i in range(sat_nbr):
        Sat_prn_list_final.append(Sat_prn_string[i*3:i*3+3])

    ### Sigma part
    Sat_sig_list_clean = []
    for sig_line in Sat_sig_list:
        sig_line_splited = sig_line.split()
        sig_line_splited = [e for e in sig_line_splited if not "+" in e]
        Sat_sig_list_clean = Sat_sig_list_clean + sig_line_splited

    Sat_sig_list_final = [int(e) for e in Sat_sig_list_clean[:sat_nbr]]


    ### Export part
    AC_name_list = [ac_name] * sat_nbr
    Date_list    = [date] * sat_nbr

    Header_DF = pd.DataFrame(list(zip(AC_name_list,Sat_prn_list_final,
                                      Sat_sig_list_final,Date_list)),
                             columns=["AC","sat","sigma","epoch"])

    return Header_DF


















