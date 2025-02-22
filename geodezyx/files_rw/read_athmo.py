#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 10:52:09 2021

@author: psakic
"""

import re

########## BEGIN IMPORT ##########
#### External modules
import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
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
                epoc.append(conv.year_decimal2dt(fields[1]))
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
    
    Parameters
    ----------
    trpfile:
        Filename of TRP file from Bernese
    
    Returns
    -------
    DF:
        Tropospheric solutions from Bernese in Dataframe
    
    Notes
    -----
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

    Parameters
    ----------
    metfile:
        Path of RINEX Meteorological file in List/String (e.g. made with glob)

    Returns
    -------
    DF:
        Meteorological data in DataFrame

    Notes
    -----
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
