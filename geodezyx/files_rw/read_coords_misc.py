#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 11:00:45 2021

@author: psakicki
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
#### Import the logger
import logging

import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

log = logging.getLogger('geodezyx')


# def read_spotgins_materfile(materfile_path):
#     df = pd.read_csv(materfile_path, header=None, sep=r'\s+', comment='#',
#                 names=['AC', 'NAME', 'ID', 'LATITUDE', 'LONGITUDE', 'HEIGHT', 'X_position',
#                        'Y_position', 'Z_position', 'LOGFILE_NAME', 'DATA_SOURCE'])
#     return df


def read_spotgins_quick(p):
    df = pd.read_csv(p , comment='#', header=None, sep=r'\s+')
    

    with open(p) as f:
        header_line = [line for line in f if line.startswith('#')][-1]
    df.columns = header_line.strip().split()
    df = df[df.columns[:4]]
    df.columns = ['mjd',"E","N","U"]
    df['date'] = conv.mjd2dt(df['mjd'])    
    df.set_index('date', inplace=True)
    df = df.drop("mjd",axis=1)

    return df


def diff_spotgins_quick(df1,df2):
    def _rndidx(d):
        duse = d.index.to_series().dt.round("1h").values
        dout = d.set_index(duse)
        return dout
    
    def _jjul(d):
        d = _rndidx(d)
        duse = d.index.to_series()
        duse = conv.numpy_dt2dt(duse)
        jjul = conv.dt2jjul_cnes(duse, True)
                
        dout = d.set_index(jjul) # + jjul[1]/86400)
        return dout 

    lbda = _rndidx
            
    dfout = (lbda(df1) - lbda(df2)).dropna()
    #print(dfout.to_string())
    return dfout

def stations_in_epos_sta_coords_file_mono(coords_file_path):
    """
    Gives stations in a EPOS coords. file (YYYY_DDD_sta_coordinates)

    Parameters
    ----------
    coords_file_path : str
        path of the EPOS coords. file.

    Returns
    -------
    epoch : datetime
        the main mean epoch in the EPOS coords. file.
    stats_list : list
        list of 4 char station list.
    """

    site_line_list = utils.grep(coords_file_path , " SITE            m")

    stats_list = []
    mean_mjd_list = []
    for l in site_line_list:
        stat = l.split()[8].lower()
        stats_list.append(stat)
        mean_mjd = np.mean([float(l.split()[6]) , float(l.split()[7])])
        mean_mjd_list.append(mean_mjd)

    mjd_final = utils.most_common(mean_mjd_list)
    epoch = conv.mjd2dt(mjd_final)
    return epoch , stats_list


def stations_in_sinex_mono(sinex_path):
    """
    Gives stations list in a SINEX file

    Parameters
    ----------
    sinex_path : str
         path of the SINEX file.

    Returns
    -------
    epoch : datetime
        the main mean epoch in the SINEX.
    stats_list : list
        list of 4 char station list.

    """

    extract = utils.extract_text_between_elements(sinex_path,'+SITE/ID','-SITE/ID')
    extract = extract.split('\n')
    extract2 = []
    for e in extract:
        if e != '' and e[0] == ' ' and e != '\n':
            extract2.append(e)

    stats_list = [e.split()[0].lower() for e in extract2]

    extract = utils.extract_text_between_elements(sinex_path,'+SOLUTION/EPOCHS','-SOLUTION/EPOCHS')
    extract = extract.split('\n')
    extract2 = []
    for e in extract:
        if e != '' and e[0] == ' ' and e != '\n':
            extract2.append(e)

    epoch = conv.datestr_sinex_2_dt(utils.most_common([e.split()[-1] for e in extract2]))

    return epoch , stats_list


def stations_in_coords_file_multi(files_path_list,files_type = "sinex"):
    """
    Gives stations list in a SINEX or EPOS coords. file

    Parameters
    ----------
    files_path_list : iterable (list)
        path of the SINEX/coords files list (e.g. made with glob).
    files_type : str, optional
        "sinex" or "EPOS_sta_coords" depending on the type. The default is "sinex".

    Returns
    -------
    datadico : dict.
        a dico with stat 4 char. code as key and epoch list as values (for timeline_plotter)

    """
    

    if   files_type == "sinex":
        extract_fct = stations_in_sinex_mono
    elif files_type == "EPOS_sta_coords":
        extract_fct = stations_in_epos_sta_coords_file_mono
    else:
        log.error("check station file type !!")
        return None

    epoch_list_stk , stats_list_stk = [] , []

    for file_path in files_path_list:
        try:
            epoch , stats_list = extract_fct(file_path)
            epoch_list = [epoch] * len(stats_list)
            epoch_list_stk = epoch_list_stk + epoch_list
            stats_list_stk = stats_list_stk + stats_list
        except:
            log.error('something wrong with' , file_path)
            log.error('TIPS: good files type ? ("sinex" or "EPOS_sta_coords"): ', files_type)

            continue

    datadico = dict()
    for e,s in zip(epoch_list_stk,stats_list_stk):
        if not s in datadico.keys():
            datadico[s] = []
        else:
            datadico[s].append(e)

    for v in datadico.values():
        v.sort()

    return datadico


def stations_in_sinex_multi(sinex_path_list):
    """
    Gives stations list in a SINEX

    Just a wrapper of stations_in_sinex_or_EPOS_sta_coords_file_multi !!!
    This other function should be used !!!

    Args :
        sinex_path_list : path of the SINEX files list (e.g. made with glob)
    Returns :
        datadico : a dico with stat 4 char. code as key and epoch list as values (for timeline_plotter)
    """

    return stations_in_coords_file_multi(sinex_path_list,"sinex")


def sinex_bench_antenna_df_2_disconts(df_antenna_in, stat, return_full=False):
    df_antenna_work = df_antenna_in[df_antenna_in["Code"] == stat]
    start_list     = conv.datestr_sinex_2_dt(df_antenna_work["_Data_Start"])
    end_list       = conv.datestr_sinex_2_dt(df_antenna_work["_Data_End__"])
    if return_full:
        return start_list,end_list
    else:
        clean_list = sorted(list(set(start_list + end_list)))
        clean_list = [e for e in clean_list if e != dt.datetime(1970, 1, 1, 0, 0)]
        return clean_list
