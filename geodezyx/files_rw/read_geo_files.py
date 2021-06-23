#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.files_rw contains reading functions to 
import files containing geodetic observations/products.

it can be imported directly with:
from geodezyx import files_rw

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU GPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (GFZ, pierre.sakic@gfz-postdam.de)
GitHub repository :
https://github.com/GeodeZYX/GeodeZYX-Toolbox_v4
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
import gzip
import linecache
import io
import numpy as np
import os 
import pandas as pd
import re


#### geodeZYX modules
from geodezyx import conv
from geodezyx import time_series
from geodezyx import utils

#### Import star style
# from geodezyx import *                   # Import the GeodeZYX modules
# from geodezyx.externlib import *         # Import the external modules
# from geodezyx.megalib.megalib import *   # Import the legacy modules names

##########  END IMPORT  ##########

    
############### Reading Tropospheric ################################################

 #  _______                              _                     ______ _ _           
 # |__   __|                            | |                   |  ____(_) |          
 #    | |_ __ ___  _ __   ___  ___ _ __ | |__   ___ _ __ ___  | |__   _| | ___  ___ 
 #    | | '__/ _ \| '_ \ / _ \/ __| '_ \| '_ \ / _ \ '__/ _ \ |  __| | | |/ _ \/ __|
 #    | | | | (_) | |_) | (_) \__ \ |_) | | | |  __/ | |  __/ | |    | | |  __/\__ \
 #    |_|_|  \___/| .__/ \___/|___/ .__/|_| |_|\___|_|  \___| |_|    |_|_|\___||___/
 #                | |             | |                                               
 #                |_|             |_|        
                
                
def read_snx_trop(snxfile,dataframe_output=True,version=2):
    """
    Read troposphere solutions from Troposphere SINEX
    """
    
    STAT , epoc = [] , []
    tro , stro , tgn , stgn , tge , stge = [] , [] , [] , [] , [] , []
    
    flagtrop = False
    
    for line in open(snxfile,"r",encoding = "ISO-8859-1"):

        if re.compile('TROP/SOLUTION').search(line):
            flagtrop = not flagtrop
            continue
        
        if line[0] == ' ':
            fields = line.split()
        else:
            continue
        
        if flagtrop ==True:
            
            STAT.append(fields[0].upper())
            
            if not ':' in fields[1]:
                epoc.append(conv.convert_partial_year(fields[1]))
            else:
                date_elts_lis = fields[1].split(':')
                if version == 2:
                    yy =  int(date_elts_lis[0])
                else:
                    yy =  int(date_elts_lis[0]) + 2000
                    
                doy = int(date_elts_lis[1])
                sec = int(date_elts_lis[2])
                epoc.append(conv.doy2dt(yy,doy,seconds=sec))
                
            if len(fields) == 8:
                tro.append(np.nan if '*' in fields[2] else fields[2])
                stro.append(np.nan if '*' in fields[3] else fields[3])
                tgn.append(np.nan if '*' in fields[4] else fields[4])
                stgn.append(np.nan if '*' in fields[5] else fields[5])
                tge.append(np.nan if '*' in fields[6] else fields[6])
                stge.append(np.nan if '*' in fields[7] else fields[7])
                            
            elif len(fields) == 4:
                tro.append(np.nan if '*' in fields[2] else fields[2])
                stro.append(np.nan if '*' in fields[3] else fields[3])
                tgn.append(np.nan)
                stgn.append(np.nan)
                tge.append(np.nan)
                stge.append(np.nan)
                
            else:
                tro.append(np.nan)
                stro.append(np.nan)
                tgn.append(np.nan)
                stgn.append(np.nan)
                tge.append(np.nan)
                stge.append(np.nan)
    
    outtuple = \
    list(zip(*sorted(zip(STAT , epoc , tro , stro , tgn , stgn , tge , stge))))
    
    if dataframe_output:
        return Tropsinex_DataFrame(outtuple)
                
def read_gfz_trop(trpfile):
    """
    This function is reading GFZ troposphere sinex into Pandas DataFrame

    Parameters
    ----------
    trpfile : Str
        File name of GFZ troposphere sinex.

    Returns
    -------
    DF : Pandas DataFrame
        Pandas Dataframe of GFZ troposphere sinex.

    """
    fields = []
    flagtrop = False

    for line in open(trpfile,"r",encoding = "ISO-8859-1"):
        if re.compile('TROP/SOLUTION').search(line):
            flagtrop = not flagtrop
            continue

        if flagtrop ==True and line[0] == ' ':
            field = line.split()
            fields.append(field)
        else:
            continue

    DF = pd.DataFrame(fields)
    DF.drop(columns=[0,8,9,10,11,12,18,19,20],inplace=True)
    DF.columns = ['STAT','epoc','year','doy','secofday','ztd_est','ztd_est_std','num_sat','tgn_est','tgn_est_std','tge_est','tge_est_std']
    cols_numeric = ['epoc','ztd_est','ztd_est_std','num_sat','tgn_est','tgn_est_std','tge_est','tge_est_std']
    DF[cols_numeric] = DF[cols_numeric].apply(pd.to_numeric, errors='coerce')
    DF['epoc'] = conv.MJD2dt(DF['epoc'].values)
    DF['epoc'] = DF['epoc'].dt.floor('H')
    return DF


def Tropsinex_DataFrame(read_sinex_result):
     """
      General description

        Parameters
        ----------
        read_sinex_result : 
            List values from read_snx_trop function
                    
        Returns
        -------
        DF_Sinex : 
            troposphere information from SINEX in Dataframe
            
     """
     DF_Sinex = pd.DataFrame.from_records(list(read_sinex_result)).transpose()
     colnam = ['STAT','epoc','tro','stro','tgn','stgn','tge','stge']
     DF_Sinex.columns = colnam
     cols_numeric = ['tro','stro','tgn','stgn','tge','stge']
     DF_Sinex[cols_numeric] = DF_Sinex[cols_numeric].apply(pd.to_numeric, errors='coerce')
     
     return DF_Sinex

def read_bernese_trp(trpfile):
    """
    This function reads tropospheric solution in TRP format from Bernese
    GNSS software
    
    Parameter
    ----------
    trpfile:
        Filename of TRP file from Bernese
    
    Return
    ----------
    DF:
        Tropospheric solutions from Bernese in Dataframe
    
    Notes
    ----------
        Written by Chaiyaporn Kitpracha
    """
    flagtrop = False
    field = []
    for line in open(trpfile,"r",encoding = "ISO-8859-1"):
        if re.compile('STATION NAME').search(line):
            headers = line.split()
            headers.remove('YYYY')
            headers.remove('MM')
            headers.remove('DD')
            headers.remove('HH')
            headers.remove('MM')
            headers.remove('SS')
            headers[3] = 'year'
            headers[4] = 'month'
            headers[5] = 'day'
            headers[6] = 'hour'
            headers[7] = 'minute'
            headers[8] = 'second'
            flagtrop = True
            continue
        
        if flagtrop and not line == '\n':
            fields = line.split()
            field.append(fields)
        else:
            continue
        
    DF = pd.DataFrame(field,columns=headers)
    DF['dt'] = pd.to_datetime(DF[['year','month','day','hour','minute','second']])
    cols_num = ['MOD_U', 'CORR_U', 'SIGMA_U', 'TOTAL_U', 'CORR_N', 'SIGMA_N','CORR_E', 'SIGMA_E']
    DF[cols_num] = DF[cols_num].apply(pd.to_numeric, errors='coerce')
    DF.drop(['year','month','day','hour','minute','second'], axis=1,inplace=True)
    return DF

def read_rinex_met(metfile):
    """
    This function reads RINEX Meteorological files and convert to Pandas DataFrame

    Parameter
    ----------
    metfile:
        Path of RINEX Meteorological file in List/String (e.g. made with glob)

    Return
    ----------
    DF:
        Meteorological data in DataFrame

    Notes
    ----------
        Written by Chaiyaporn Kitpracha
    """
    if utils.is_iterable(metfile):
        merge_df = pd.DataFrame()
        for metfile_m in metfile:
            met_df = read_rinex_met_2(str(metfile_m))
            merge_df = pd.concat([merge_df,met_df])
        return merge_df
    else:
        met_df = read_rinex_met_2(metfile)
        return met_df

def read_rinex_met_2(metfile):
    ln=0
    for line in open(metfile,"r",encoding = "ISO-8859-1"):
        if re.compile('MARKER NAME').search(line):
            marker = line.split()[0]
            marker = marker.upper()
        if re.compile('# / TYPES OF OBSERV').search(line):
            tmp = line.split()
            headers = tmp[1:int(tmp[0])+1]
        if re.compile('TD SENSOR MOD/TYPE/ACC').search(line):
            tmp = line.split()
            temp_unc = float(tmp[-4])
        if re.compile('PR SENSOR MOD/TYPE/ACC').search(line):
            tmp = line.split()
            press_unc = float(tmp[-4])
        if re.compile('HR SENSOR MOD/TYPE/ACC').search(line):
            tmp = line.split()
            humrel_unc = float(tmp[-4])
        if re.compile('END OF HEADER').search(line):
            break
        ln = ln+1

    df = pd.read_csv(metfile,skiprows=range(0,ln+1),delim_whitespace=True,names=['year','month','day','hour','minute','second']+headers)
    df['year'] = df['year'] + 2000 if df['year'].any() <= 79 else df['year'].any() + 1900
    df['STA'] = marker
    df['epoch'] = pd.to_datetime(df[['year','month','day','hour','minute','second']],errors='coerce')
    df.drop(['year','month','day','hour','minute','second'], axis=1,inplace=True)
    if press_unc is not None:
        df['PR_std'] = press_unc
    if temp_unc is not None:
        df['TD_std'] = temp_unc
    if humrel_unc is not None:
        df['HR_std'] = humrel_unc
    df.set_index('epoch',inplace=True)
    return df







 #  ______                _   _                _____                                         _ 
 # |  ____|              | | (_)              / ____|                                       | |
 # | |__ _   _ _ __   ___| |_ _  ___  _ __   | |  __ _ __ __ ___   _____ _   _  __ _ _ __ __| |
 # |  __| | | | '_ \ / __| __| |/ _ \| '_ \  | | |_ | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
 # | |  | |_| | | | | (__| |_| | (_) | | | | | |__| | | | (_| |\ V /  __/ |_| | (_| | | | (_| |
 # |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \_____|_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
 #                                                                        __/ |                
 #                                                                       |___/       



def read_erp_bad(path,return_array=False):
    """
    This function is discontinued, use read_erp1 instead
    """
    F = open(path)
    L = []
    for l in F:
        if "MJD" in l:
            keys = l.split()
        try:
            float(l[0])
        except:
            continue
        l=utils.str_2_float_line(l.strip())
        L.append(l)

    M = np.vstack(L)

    if return_array:
        return M
    else:
        return pd.DataFrame(M)


############################# READ CLK 30
#def read_clk1(file_path_in, returns_pandas = True, interval=None):
#    """
#    General description
#
#    Parameters
#    ----------
#    file_path_in :  str
#        Path of the file in the local machine.
#
#    returns_pandas :  bool
#        Define if pandas function will be used to put the data on tables
#
#    interval :  int
#        Defines which interval should be used in the data tables. The interval is always in minutes unit.
#
#    Returns
#    -------
#    out1 : float or int oandas table
#        Returns a panda table format with the data extracted from the file.
#
#    out2 :  list
#       In case of the pandas function are not in use, the return will be a list with the data extract from the file.
#
#
#    """
#
#    file = open(file_path_in,'r')
#    fil = file.readlines()
#    file.close()
#
#
#    types_m, name_m, year_m, month_m, day_m, hr_m, minute_m, seconds_m, epoch_m, offset_m, rms_m = [],[],[],[],[],[],[],[],[],[],[]
#
#    Clk_read  = []
#
#        #####IGNORES HEADER
#
#    le = 0
#    i = 0
#    count = 0
#    for le in range(len(fil)):
#     linhaatual = linecache.getline(file_path_in, le)
#     header_end = (linhaatual[60:73])
#     count +=1
#     if header_end =="END OF HEADER":
#        i = count
##############################################################
#    while i <= len(fil):
#          linhaatual = linecache.getline(file_path_in, i)
#          types = (linhaatual[0:2]) # Column refers to AS or AR
#          name = (linhaatual[3:7])# Name of the station or satellite
#          year = int((linhaatual[8:12]))
#          month = int((linhaatual[13:15]))
#          day = int((linhaatual[16:18]))
#          hr= int((linhaatual[19:22]))# hour
#          minute = int((linhaatual[22:25]))
#          seconds = float((linhaatual[25:33]))
#          offset = float(((linhaatual[40:59])))# Clock offset
#          ##### check if there is a value for thr rms
#          if (linhaatual[65:70]) == '     ' or len(linhaatual[65:70])==0:
#              rms = 0
#          else:
#              rms = float(linhaatual[61:79])
#          ############
#          epoch = dt.datetime(year,month,day,hr,minute,int(seconds),int(((int(seconds)-seconds)*10**-6))) # epoch as date.time function
#
#        ############ Data in a defined interval
#          if interval:
#              if (float(minute)%interval) == 0 and float(seconds) == 0:
#
#                  if returns_pandas:
#                            Clk_dados = [types,name,year,month,day,hr,minute,seconds,epoch,offset,rms]
#                            Clk_read.append(Clk_dados)
#                  else:
#                            types_m.append(types)
#                            name_m.append(name)
#                            year_m.append(year)
#                            month_m.append(month)
#                            day_m.append(day)
#                            hr_m.append(hr)
#                            minute_m.append(minute)
#                            seconds_m.append(seconds)
#                            epoch_m.append(epoch)
#                            offset_m.append(offset)
#                            rms_m.append(rms)
#         #########################################################
#        ########################################################### standart file epochs
#          else:
#              if returns_pandas:
#                            Clk_dados = [types,name,year,month,day,hr,minute,seconds,epoch,offset,rms]
#                            Clk_read.append(Clk_dados)
#              else:
#                            types_m.append(types)
#                            name_m.append(name)
#                            year_m.append(year)
#                            month_m.append(month)
#                            day_m.append(day)
#                            hr_m.append(hr)
#                            minute_m.append(minute)
#                            seconds_m.append(seconds)
#                            epoch_m.append(epoch)
#                            offset_m.append(offset)
#                            rms_m.append(rms)
#
#
#          i +=1
#
#  ################################################################################################
#   ############################################ put on pandas table format
#    if returns_pandas:
#     Clk_readed = pd.DataFrame(Clk_read, columns=['type','name','year','month','day','hr','minutes','seconds','epoch','offset','rms'])
#     Clk_readed.path = file_path_in
#
#     return Clk_readed
#
#          ###############################
#    else:
#           print(" ...")
#           return  types_m, name_m, year_m, month_m, day_m, hr_m, minute_m, seconds_m, offset_m, rms_m
#



############################# END FUNCTION READ CLK 30




################################################READ ERP GUS
# def read_erp1(caminho_arq,ac):

    # le = open(caminho_arq, 'r')
    # letudo = le.readlines()
    # le.close()
    # tamanho = len(letudo) #usado para saber quantas linhas tem o arquivo



    # para = tamanho #onde o arquivo para de ser lido

    # numeros = ['0','1','2','3','4','5','6','7','8','9']
    # le = 0
    # numlin = 0 #numero de linhas da matriz de epocas
    # numcol = 16 #numero de colunas que a matriz final deve ter


    # while le <= para:
     # linhaatual = linecache.getline(caminho_arq, le)
     # if linhaatual[0:1] in numeros:
         # numlin +=1
     # le +=1

    # dt = np.dtype(object)
    # ERP = np.empty((numlin,numcol), dtype=dt)
    # n = 0
    # j = 0
    # l = 0
    # g = 0
# #####################################
    # if ac == 'wum':
        # while l <=para:
         # linhaatual = linecache.getline(caminho_arq, l)
         # if linhaatual[0:1] in numeros:
             # ERP[j,0] = float(linhaatual[0:8])
             # ERP[j,1] = round((float(linhaatual[11:17])*(10**-6)),6)
             # ERP[j,2] = float(linhaatual[18:25])*(10**-6)
             # ERP[j,4] = float(linhaatual[26:34])*(10**-7)
             # ERP[j,3] = float(linhaatual[35:41])*(10**-7)
             # ERP[j,5] = float(linhaatual[44:47])*(10**-6)
             # ERP[j,6] = float(linhaatual[50:53])*(10**-6)
             # ERP[j,7] = float(linhaatual[57:61])*(10**-7)
             # ERP[j,8] = float(linhaatual[63:69])*(10**-7)
             # ERP[j,9] = float(linhaatual[69:73])
             # ERP[j,10] = float(linhaatual[74:76])
             # ERP[j,11] = float(linhaatual[76:80])
             # ERP[j,12] = float(linhaatual[81:86])*(10**-6)
             # ERP[j,13] = float(linhaatual[88:93])*(10**-6)
             # ERP[j,14] = float(linhaatual[96:101])*(10**-6)
             # ERP[j,15] = float(linhaatual[102:107])*(10**-6)

             # j +=1
         # l +=1

        # Erp_end = pd.DataFrame(ERP, columns=['MJD','X-P', 'Y-P', 'UT1UTC','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
                                             # 'X-RT','Y-RT','S-XR','S-YR'])
    # #    Erp_end.set_index('MJD',inplace=True)
        # return Erp_end
# #######################################
    # if ac == 'COD':
        # while n <=para:
         # linhaatual = linecache.getline(caminho_arq, n)
         # if linhaatual[0:1] in numeros:
             # ERP[j,0] = float(linhaatual[0:8])
             # ERP[j,1] = round((float(linhaatual[12:17])*(10**-6)),6)
             # ERP[j,2] = float(linhaatual[20:27])*(10**-6)
             # ERP[j,4] = float(linhaatual[38:41])*(10**-7)
             # ERP[j,3] = float(linhaatual[29:35])*(10**-7)
             # ERP[j,5] = float(linhaatual[47:49])*(10**-6)
             # ERP[j,6] = float(linhaatual[53:55])*(10**-6)
             # ERP[j,7] = float(linhaatual[59:61])*(10**-7)
             # ERP[j,8] = float(linhaatual[65:66])*(10**-7)
             # ERP[j,9] = float(linhaatual[67:70])
             # ERP[j,10] = float(linhaatual[71:74])
             # ERP[j,11] = float(linhaatual[75:77])
             # ERP[j,12] = float(linhaatual[81:85])*(10**-6)
             # ERP[j,13] = float(linhaatual[88:92])*(10**-6)
             # ERP[j,14] = float(linhaatual[94:98])*(10**-6)
             # ERP[j,15] = float(linhaatual[100:104])*(10**-6)

             # j +=1
         # n +=1

        # Erp_end = pd.DataFrame(ERP, columns=['MJD','X-P', 'Y-P', 'UT1UTC','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
                                             # 'X-RT','Y-RT','S-XR','S-YR'])
    # #    Erp_end.set_index('MJD',inplace=True)
        # return Erp_end
# ################################################
    # if ac == 'gbm':
        # while g <=para:
         # linhaatual = linecache.getline(caminho_arq, g)
         # if linhaatual[0:1] in numeros:
             # ERP[j,0] = float(linhaatual[0:9])
             # ERP[j,1] = round((float(linhaatual[12:18])*(10**-6)),6)
             # ERP[j,2] = float(linhaatual[19:25])*(10**-6)
             # ERP[j,4] = float(linhaatual[26:37])*(10**-7)  ###UT1 - TAI!!!!!!!!!!!!
             # ERP[j,3] = float(linhaatual[39:44])*(10**-7)
             # ERP[j,5] = float(linhaatual[46:49])*(10**-6)
             # ERP[j,6] = float(linhaatual[50:54])*(10**-6)
             # ERP[j,7] = float(linhaatual[56:59])*(10**-7)
             # ERP[j,8] = float(linhaatual[60:64])*(10**-7)
             # ERP[j,9] = float(linhaatual[65:70])
             # ERP[j,10] = float(linhaatual[71:72])
             # ERP[j,11] = float(linhaatual[73:76])
             # ERP[j,12] = float(linhaatual[77:82])*(10**-6)
             # ERP[j,13] = float(linhaatual[83:88])*(10**-6)
             # ERP[j,14] = float(linhaatual[92:96])*(10**-6)
             # ERP[j,15] = float(linhaatual[100:104])*(10**-6)

             # j +=1
         # g +=1

        # Erp_end = pd.DataFrame(ERP, columns=['MJD','X-P', 'Y-P', 'UT1UTC','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
                                             # 'X-RT','Y-RT','S-XR','S-YR'])
    # #    Erp_end.set_index('MJD',inplace=True)
        # return Erp_end
        # #out = '/home/mansur/Documents/MGEX/Saida/outERP_'+caminho_arq[28:-4]+'.txt'
        # #np.savetxt(out, ERP, fmt="%02s %04s %05s %05s %05s %05s %05s %10s %10s %10s %10s %10s %10s %10s %10s %10s")
################################################### END READ ERP GUS
##########################################################################################################################

def read_clk_old(file_path_in, returns_pandas = True, interval=None,old_naming=True):
    """
    Read an IGS clk file
    Slow and complex, use read_clk instead
    
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
        if old_naming:
            name_list = ['offset','rms']
        else:
            name_list = ['bias','sigma']
            
        Clk_readed = pd.DataFrame(Clk_read, columns= ['type','name','year','month','day','h','minutes','seconds','epoch'] + name_list)
        Clk_readed['epoch'] =  pd.to_datetime(Clk_readed[[ 'year' ,'month' ,'day','h','minutes','seconds']])
        Clk_readed.path = file_path_in
    
        return Clk_readed

          ###############################
    else:
           print(" ...")
           return  types_m, name_m, year_m, month_m, day_m, hr_m, minute_m, seconds_m, offset_m, rms_m




############################# END FUNCTION READ CLK 30





# def read_erp(caminho_arq,ac):
    # """
    # General description
    # Units: ('MJD','X-P (arcsec)', 'Y-P (arcsec)', 'UT1UTC (E-7S)','LOD (E-7S/D)','S-X (E-6" arcsec)','S-Y (E-6" arcsec)',
    # 'S-UT (E-7S)','S-LD (E-7S/D)','NR (E-6" arcsec)', 'NF (E-6" arcsec)', 'NT (E-6" arcsec)',
    # 'X-RT (arcsec/D)','Y-RT (arcsec/D)','S-XR (E-6" arcsec/D)','S-YR (E-6" arcsec/D)', 'C-XY', 'C-XT',
    # 'C-YT', 'DPSI', 'DEPS','S-DP','S-DE')

    # Parameters
    # ----------
    # file_path_in :  str
        # Path of the file in the local machine.

    # which AC :  str
        # The analisys center that will be used


    # Returns
    # -------
    # out1 :  pandas table
        # Returns a panda table format with the data extracted from the file.

    # Obs: this function need to be improved to read also the values of UTC and etc. So far reads only the Pole rates. 

    # """

    # #### FIND DELIVERY DATE
    # name = os.path.basename(caminho_arq)
    
    # if len(name) == 12:
        # dt_delivery = conv.sp3name2dt(caminho_arq)
    # elif len(name) == 38:
        # dt_delivery = conv.sp3name_v3_2dt(caminho_arq)
    # else:
        # dt_delivery = conv.posix2dt(0)


    # le = open(caminho_arq, 'r')
    # letudo = le.readlines()
    # le.close()
    # tamanho = len(letudo) 

    

    # numeros = ['0','1','2','3','4','5','6','7','8','9']


    # ERP=[]


    # if caminho_arq[-3:] in ('snx','ssc'):
        # file = open(caminho_arq)
        # Lines =  file.readlines()
        # XPO_stk  = []
        # XPO_std_stk = []
        # YPO_stk  = []
        # YPO_std_stk = []
        # LOD_stk  = []
        # LOD_std_stk = []
        # MJD_stk = []
        # marker = False

        # for i in range(len(Lines)):

            # if len(Lines[i].strip()) == 0:
                # continue
            # else:

                # if Lines[i].split()[0] == '+SOLUTION/ESTIMATE':
                    # marker = True

                # if Lines[i].split()[0] == '-SOLUTION/ESTIMATE':
                    # marker = False

                # if utils.contains_word(Lines[i],'XPO') and marker:
                    # Doy = (Lines[i][30:33])
                    # Year = (Lines[i][27:29])
                    # Pref_year = '20'
                    # Year = int(Pref_year+Year)
                    # Date = conv.doy2dt(Year,Doy)
                    # XPO = float(Lines[i][47:68])*(10**-3)
                    # XPO_std = float(Lines[i][69:80])*(10**-3)
                    # XPO_stk.append(XPO)
                    # XPO_std_stk.append(XPO_std)
                    # MJD_stk.append(conv.jd_to_mjd(conv.date_to_jd(Date.year,Date.month,Date.day)))
                    
                # if utils.contains_word(Lines[i],'YPO') and marker:
                    # Doy = (Lines[i][30:33])
                    # Year = str(Lines[i][27:29])
                    # Pref_year = '20'
                    # Year = int(Pref_year+Year)
                    # Date = conv.doy2dt(Year,Doy)
                    # YPO = float(Lines[i][47:68])*(10**-3)
                    # YPO_std = float(Lines[i][69:80])*(10**-3)
                    # YPO_stk.append(YPO)
                    # YPO_std_stk.append(YPO_std)
                    # MJD_stk.append(conv.jd_to_mjd(conv.date_to_jd(Date.year,Date.month,Date.day)))
                    
                # if utils.contains_word(Lines[i],'LOD') and marker:
                    # Doy = (Lines[i][30:33])
                    # Year = str(Lines[i][27:29])
                    # Pref_year = '20'
                    # Year = int(Pref_year+Year)
                    # Date = conv.doy2dt(Year,Doy)
                    # LOD = float(Lines[i][47:68])*(10**+4)
                    # LOD_std = float(Lines[i][69:80])*(10**+4)
                    # LOD_stk.append(LOD)
                    # LOD_std_stk.append(LOD_std)
                    # MJD_stk.append(conv.jd_to_mjd(conv.date_to_jd(Date.year,Date.month,Date.day)))
                    
        # MJD = list(set(MJD_stk))
        # if len(LOD_stk) == 0:
                # LOD_stk = ['0']*len(MJD)
                # LOD_std_stk = ['0']*len(MJD)

        # for i in range(len(MJD)):

            # ERP_data = [ac, MJD[i], XPO_stk[i], YPO_stk[i], 0, LOD_stk[i], XPO_std_stk[i], YPO_std_stk[i],
                     # 0, LOD_std_stk[i], 0, 0, 0, 0, 0, 0, 0, dt_delivery]
            # print(ERP_stk)
            # ERP.append(ERP_data)



    # if ac in ('COD','cod','com', 'cof', 'grg', 'mit', 'sio'):
        # for i in range(tamanho+1):
            # linhaatual = linecache.getline(caminho_arq, i)
            # if linhaatual[0:1] in numeros:
                # ERP_data = linhaatual.split()
                # for j in range(len(ERP_data)):
                    # ERP_data[j] = float(ERP_data[j])
                # ERP_data.insert(0,ac)
                # ERP_data[2] =  ERP_data[2]*(10**-6)
                # ERP_data[3] =  ERP_data[3]*(10**-6)
                # ERP_data[13] = ERP_data[13]*(10**-6)
                # ERP_data[14] = ERP_data[14]*(10**-6)
                # del ERP_data[17:]
                # ERP_data.append(dt_delivery)

                # ERP.append(ERP_data)



    # if ac in ('wum','grg','esa', 'mit', 'ngs', 'sio'):
        # for i in range(tamanho+1):
            # linhaatual = linecache.getline(caminho_arq, i)
            # if linhaatual[0:1] in numeros:
                # ERP_data = linhaatual.split()
                # for j in range(len(ERP_data)):
                    # ERP_data[j] = float(ERP_data[j])
                # ERP_data.insert(0,ac)
                # ERP_data[2] =  ERP_data[2]*(10**-6)
                # ERP_data[3] =  ERP_data[3]*(10**-6)
                # ERP_data[13] = ERP_data[13]*(10**-6)
                # ERP_data[14] = ERP_data[14]*(10**-6)
                # ERP_data.append(dt_delivery)

                # ERP.append(ERP_data)

# #
    # if ac in ('gbm', 'gfz'):
        # for i in range(tamanho+1):
            # linhaatual = linecache.getline(caminho_arq, i)
            # if linhaatual[0:1] in numeros:
                # ERP_data = linhaatual.split()
                # for j in range(len(ERP_data)):
                    # ERP_data[j] = float(ERP_data[j])
                # ERP_data.insert(0,ac)
                # ERP_data[2] =  ERP_data[2]*(10**-6)
                # ERP_data[3] =  ERP_data[3]*(10**-6)
                # ERP_data[13] = ERP_data[13]*(10**-6)
                # ERP_data[14] = ERP_data[14]*(10**-6)
                # ERP_data.append(dt_delivery)

                # ERP.append(ERP_data)


    # header = []
    # if ac in ('emr'):
        # for i in range(tamanho+1):
            # linhaatual = linecache.getline(caminho_arq, i)
            # if linhaatual == 'EOP  SOLUTION':
                # del ERP_data[:]
                # header = ['EOP  SOLUTION']
            # if linhaatual[0:1] in numeros and 'EOP  SOLUTION' in header:
                # ERP_data = linhaatual.split()
                # for j in range(len(ERP_data)):
                    # ERP_data[j] = float(ERP_data[j])
                # ERP_data.insert(0,ac)
                # ERP_data[2] =  ERP_data[2]*(10**-6)
                # ERP_data[3] =  ERP_data[3]*(10**-6)
                # ERP_data[13] = ERP_data[13]*(10**-6)
                # ERP_data[14] = ERP_data[14]*(10**-6)
                # del ERP_data[17:]
                # ERP_data.append(dt_delivery)

                # ERP.append(ERP_data)

    # Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD',
                                         # 'NR', 'NF', 'NT',
                                         # 'X-RT','Y-RT','S-XR','S-YR',
                                         # 'Delivered_date'])
    # return Erp_end


def list_files(dire,file = None):
    """
    written for MGEX combi by GM
    should be discontinued
    """
    List_arq = os.listdir(dire)

    path = []

    for l in range(len(List_arq)):
        atual = List_arq[l]
        if file == 'sp3':
            if atual[-3:] == 'SP3' or atual[-3:] == 'sp3':
                path.append(os.path.join(dire,List_arq[l]))

        if file == 'clk':
            if atual[-3:] == 'CLK' or atual[-3:] == 'clk':
                path.append(os.path.join(dire,List_arq[l]))

        if file == 'erp':
            if atual[-3:] == 'ERP' or atual[-3:] == 'erp':
                path.append(os.path.join(dire,List_arq[l]))

        if file == 'eph':
            if atual[-3:] == 'EPH' or atual[-3:] == 'eph':
                path.append(os.path.join(dire,List_arq[l]))


        if file == 'snx':
            if atual[-3:] == 'SNX' or atual[-3:] == 'snx':
                path.append(os.path.join(dire,List_arq[l]))


    return path


 
