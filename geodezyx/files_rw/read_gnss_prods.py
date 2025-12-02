#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 10:24:19 2021

@author: psakic
"""

########## BEGIN IMPORT ##########
#### External modules
import gzip
#### Import the logger
import logging
import os

import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
# from geodezyx import time_series
from geodezyx import utils
from geodezyx.reffram import orb_df_velocity_calc

log = logging.getLogger('geodezyx')

 #   _____ _            _       __ _ _           
 #  / ____| |          | |     / _(_) |          
 # | |    | | ___   ___| | __ | |_ _| | ___  ___ 
 # | |    | |/ _ \ / __| |/ / |  _| | |/ _ \/ __|
 # | |____| | (_) | (__|   <  | | | | |  __/\__ \
 #  \_____|_|\___/ \___|_|\_\ |_| |_|_|\___||___/
                                               
                                               

def read_clk(file_path_inp,names_4char=False):
    """
    Read an IGS clk file

    Parameters
    ----------
    file_path_inp :  str
        Path of the file in the local machine.
        can handle gzip-compressed file (with .gz/.GZ extension) 
    names_4char : bool
        Force the station names to have 4 charaters
        (new naming convention is longer)

    Returns
    -------
    clk_df : pandas DataFrame
        Returns a panda table format with the data extracted from the file.
        
    Note
    ----
    Bias is given in seconds
    """
    head_line = utils.grep(file_path_inp,"END OF HEADER",
                          only_first_occur=True,line_number=True)
    
    clk_df = pd.read_csv(file_path_inp,skiprows=head_line[0]+1,header=None,
                        delim_whitespace = True,
                        names=['type', 'name', 'year', 'month', 'day', 'hour',
                               'minute', 'second',"n_values",'bias', 'sigma'])
    
    
    # Special case when D-0n instead of E-0n (e.g. IAC), 
    # then values are not converted to float and kept as generic objects...
    if clk_df['bias'].dtype == "O" and clk_df['bias'].str.match(r".*D(\+|-)\d\d$").any():
        clk_df['bias'] = clk_df['bias'].str.replace("D","E").astype(float)

    # remove EOF line    
    clk_df = clk_df[clk_df["type"] != "EOF"]

    # convert date to int
    #for d in ['year', 'month', 'day', 'hour','minute']:
    #    clk_df[d] = clk_df[d].astype(int)

    clk_df["year"] = clk_df["year"].astype(int)

    clk_df["ac"] = os.path.basename(file_path_inp)[:3]
    clk_df["name"] = clk_df["name"].str.upper()
    if names_4char:
        clk_df['name'] = clk_df['name'].str[:4]

    clk_df['epoch'] = pd.to_datetime(clk_df[['year', 'month', 'day',
                                           'hour','minute', 'second']])
    clk_df.path = file_path_inp

    return clk_df


def read_clk_from_sp3(file_path_or_orb_df_inp):
    """
    Get the clock values from a SP3 file of an Orbit DataFrame,
    formated as a Clock DataFrame

    Parameters
    ----------
    file_path_or_orb_df_inp : str or DataFrame
        the input SP3 file path or an Orbit DataFrame.

    Returns
    -------
    clk_df_sp3 : DataFrame
        Clock DataFrame.

    """
    if type(file_path_or_orb_df_inp) is str:
        orb_df = read_sp3(file_path_or_orb_df_inp)
    else:
        orb_df = file_path_or_orb_df_inp

    nlines = len(orb_df)

    ### replace 999999.999999 with NaN
    orb_df.loc[np.isclose(orb_df.clk,999999.999999),'clk'] = np.nan

    clk_df_sp3 = pd.DataFrame([["AS"] * nlines,
                                orb_df.sat,
                                orb_df.epoch.dt.year,
                                orb_df.epoch.dt.month,
                                orb_df.epoch.dt.day,
                                orb_df.epoch.dt.hour,
                                orb_df.epoch.dt.minute,
                                orb_df.epoch.dt.second,
                                [1] * nlines,
                                orb_df.clk * 10**-6,
                                [np.nan] * nlines,
                                orb_df.ac,
                                orb_df.epoch]).T

    clk_df_sp3 = clk_df_sp3.infer_objects()

    
    clk_df_sp3.columns = ['type', 'name', 'year', 'month',
                         'day', 'hour', 'minute', 'second',
                         'n_values', 'bias', 'sigma', 'ac', 'epoch']
    
    return clk_df_sp3


def clk_decimate(file_inp,file_out,step=300):
    """
    Decimate a .clk file 

    Parameters
    ----------
    file_inp : str
        path of the input clk file.
    file_out : str
        path of the output clk file.
    step : int, optional
        decimation step in sec. The default is 300.

    Returns
    -------
    file_out : str
        path of the output clk file.

    """


    fin = open(file_inp)

    good_line = True
    outline = []

    for l in fin:
        good_line = True
        if l[0:2] in ("AR","AS"):
            epoc   = conv.tup_or_lis2dt(l[8:34].strip().split())
            if np.mod(int(epoc.minute*60 + epoc.second) , step) == 0:
                good_line = True
            else:
                good_line = False

        if good_line:
            outline.append(l)

    with open(file_out,"w+") as Fout:
        for l in outline:
            Fout.write(l)

    return file_out


def clk_diff(file_1,file_2):
    """
    Compute the clock differences between two .clk file

    Parameters
    ----------
    file_1 : str or pandas DataFrame
        file 1 path.
    file_2 : str or pandas DataFrame
        file 2 path.

    Returns
    -------
    DFdiff : DataFrame
        clock differences epoch byx epoch between the two files.

    """
    if type(file_1) is str:
        df1 = read_clk(file_1)
    else:
        df1 = file_1

    if type(file_2) is str:
        df2 = read_clk(file_2)
    else:
        df2 = file_2

    name_list = sorted(set(df1['name']).intersection(set(df2['name'])))

    epoc_common_stk = []
    name_list_stk   = []
    type_list_stk   = []
    d_clk_stk        = []

    for nam in name_list:
        df1nam_bool = (df1['name'] == nam)
        df2nam_bool = (df2['name'] == nam)

        df1nam = df1[df1nam_bool]
        df2nam = df2[df2nam_bool]

        epoc_common = set(df1nam['epoc']).intersection(set(df2nam['epoc']))

        df1a = df1nam[ df1nam['epoc'].isin(epoc_common) ]
        df2a = df2nam[ df2nam['epoc'].isin(epoc_common) ]

        d_clk = np.array(df1a['clk_bias']) - np.array(df2a['clk_bias'])

        type_list_stk   += list(df1a['type'])
        name_list_stk   += [nam] * len(d_clk)
        epoc_common_stk += list(epoc_common)
        d_clk_stk        += list(d_clk)

    diff_df = pd.DataFrame((list(zip(name_list_stk,type_list_stk,
                                   epoc_common_stk,d_clk_stk))),
                              columns = ['name','type','date','diff'])
    return diff_df

 #   ____       _     _ _       _______ _____ ____     __ _ _           
 #  / __ \     | |   (_) |     / / ____|  __ \___ \   / _(_) |          
 # | |  | |_ __| |__  _| |_   / / (___ | |__) |__) | | |_ _| | ___  ___ 
 # | |  | | '__| '_ \| | __| / / \___ \|  ___/|__ <  |  _| | |/ _ \/ __|
 # | |__| | |  | |_) | | |_ / /  ____) | |    ___) | | | | | |  __/\__ \
 #  \____/|_|  |_.__/|_|\__/_/  |_____/|_|   |____/  |_| |_|_|\___||___/
                                                                      
def read_sp3(file_path_inp,returns_pandas = True, name = '',
             epoch_as_pd_index = False,km_conv_coef=1,
             skip_null_epoch=True,
             new_col_names=True,
             pos_vel_col=False):
    """
    Read a SP3 file (GNSS Orbits standard file) and return X,Y,Z coordinates
    for each satellite and for each epoch

    Parameters
    ----------
    file_path_inp : str
        path of the SP3 file
        can handle gzip-compressed file (with .gz/.GZ or .Z extension) 

    returns_pandas : bool
        if True, return a Pandas DataFrame.
        if False, return separated lists.

    name : str
        a manual name for the file

    epoch_as_pd_index : bool
        if True, the index of the output dataframe contains
        if False, it contains generic integer indexs
        
    km_conv_coef : float
        a conversion coefficient to change the units
        to get meters : 10**3
        to get milimeters : 10**6
        
    skip_null_epoch :bool
        Do not write an epoch if all sats are null (filtering)
        
    new_col_names : bool
        A legacy option to have (or not) consistent column names with read_rinex
        Default is False
        
    Returns
    -------
    df : Pandas DataFrame
        if returns_pandas == True

    epoch_stk ,  Xstk , Ystk , Zstk , Clkstk : lists
        if returns_pandas == False
        
    Note
    ----
    SP3 coordinates are usually given in ITRF, i.e. an ECEF system
    
    If the SP3 contains velocity records, this option adds a 'pv' column
    containing 'p' for position or 'v' for velocity.
    
    Warning: velocity 'v' records are in dm/s per default, and the
    same km_conv_coef coefficient as the position records will be applied!

    """
    
    ac_name =  os.path.basename(file_path_inp)[:3]

    if file_path_inp[-2:] in ("gz","GZ"):
        f = gzip.open(file_path_inp, "r+")
        lines = [e.decode('utf-8') for e in f]
    elif file_path_inp[-2:] in (".Z",):
        import ncompress
        fh = open(file_path_inp, 'rb')
        f = ncompress.decompress(fh).decode("utf-8")
        lines = f.split('\n')
    else:
        f = open(file_path_inp,'r')
        lines = f.readlines()

    header = True

    #### List/DF initialization
    epoch_stk = []
    xstk , ystk , zstk , clkstk = [],[],[],[]
    #Typestk     = []
    data_stk    = []
    ac_name_stk = []
   
    # Why this here ?!? (PSa 202104)
    # if returns_pandas:
    #     df = pd.DataFrame(data_stk, columns=['epoch','sat', 'const', 'sv','type',
    #                                        'x','y','z','clk','AC'])
        
    if not new_col_names:
        log.warning("you use old column names (not conventional) set new_col_names as True and adapt your code")
        ### DeprecationWarning
        col_names=['epoch','sat','const','sv',
                   'type','x','y','z','clk','AC']
    else:
        col_names=['epoch','prn','sys','prni',
                   'rec','x','y','z','clk','ac']

    #### read the Header as a 1st check
    header = read_sp3_header(lines,ac_name)
    if header.empty:
        log.warning("The SP3 looks empty: ",file_path_inp)
        if returns_pandas:
            df = pd.DataFrame([], columns=col_names)
            return df
        else:
            return  epoch_stk ,  xstk , ystk , zstk , clkstk , ac_name_stk

    for l in lines:
        if l[0] == '*':
            header = False

        if header:
            continue
        if 'EOF' in l:
            break

        if l[0] == '*':
            epoc   = conv.tup_or_lis2dt(l[1:].strip().split())
            
        elif len(l.strip()) == 0:
            continue
            
        else:
            sat_nat = l[1:2].strip()
            
            sat_sv  = int(l[2:4].strip())
            sat_sat = l[1:4].strip()
	    
            # QnD mode, must be imprved to detect nonfloat values
            if '*' in l[4:18] or not (l[4:18] and l[4:18].strip()):
                x = np.nan
            else:
                x   = float(l[4:18]) * km_conv_coef
            if '*' in l[18:32] or not (l[18:32] and l[18:32].strip()):
                y = np.nan
            else:
                y   = float(l[18:32])* km_conv_coef
            if '*' in l[32:46] or not (l[32:46] and l[32:46].strip()):
                z = np.nan
            else:
                z   = float(l[32:46])* km_conv_coef
            if '*' in l[46:60] or not (l[46:60] and l[46:60].strip()):
                clk = np.nan
            else:
                clk = float(l[46:60])
            
            typ = l[0]

            if returns_pandas:
                line_data = [epoc,sat_sat,sat_nat,sat_sv,typ,x,y,z,clk,ac_name]
                data_stk.append(line_data)
            else:
                epoch_stk.append(epoc)
                xstk.append(x)
                ystk.append(y)
                zstk.append(z)
                clkstk.append(clk)


    ac_name_stk = [ac_name] * len(xstk)

    if returns_pandas:
        df = pd.DataFrame(data_stk, columns=col_names)
        
        if skip_null_epoch:
            df = sp3_data_frame_zero_epoch_filter(df)

        if epoch_as_pd_index:
            df.set_index('epoch',inplace=True)
        df.filename = os.path.basename(file_path_inp)
        df.path = file_path_inp

        if name != '':
            df.name = name
        else:
            df.name = os.path.basename(file_path_inp)

        return df
    else:
        log.info("return list, very beta : no Sat. Vehicule Number info ...")
        return  epoch_stk ,  xstk , ystk , zstk , clkstk , ac_name_stk

def read_sp3_header(sp3_inp,ac_name_inp=None):
    """
    Read a SP3 file header and return a Pandas DataFrame
    with sat. PRNs and sigmas contained in the header

    Parameters
    ----------
    sp3_inp : str or list
        path of the SP3 file
        can handle gzip-compressed file (with .gz/.GZ extension) 

        can also handle the sp3 content as a list of strings
        (useful when read_sp3_header is used as a subfunction of read_sp3)
        
    ac_name_inp : str
        force the AC name
        (necessary when read_sp3_header is used as a subfunction of read_sp3)
        
    Returns
    -------
    header_df : Pandas DataFrame
        2 columns "sat", "sigma"

    Note
    -------
    More infos about the sigma
    http://acc.igs.org/orbacc.txt
    """

    if type(sp3_inp) is list:
        ### case when read_sp3_header is used as a subfunction of read_sp3
        lines = sp3_inp
    elif sp3_inp[-2:] in ("gz","GZ"):
        F = gzip.open(sp3_inp, "r+")
        lines = [e.decode('utf-8') for e in F]
    else:
        F = open(sp3_inp,'r+')
        lines = F.readlines()
    
    if not ac_name_inp:
        ac_name = os.path.basename(sp3_inp)[:3]
    else:
        ac_name = ac_name_inp

    sat_prn_list = []
    sat_sig_list = []

    for il , l in enumerate(lines):
        if il == 1:
            date = conv.mjd2dt(int(l.split()[4]))
        if l[:2] == "+ ":
            sat_prn_list.append(l)
        if l[:2] == "++":
            sat_sig_list.append(l)
        if l[0] == "*":
            break

    ### PRN part
    sat_prn_list_clean = []
    for prn_line in sat_prn_list:
        prn_line_splited = prn_line.split()
        prn_line_splited = [e for e in prn_line_splited if not "+" in e]
        prn_line_splited = [e for e in prn_line_splited if not  e == "0"]
        sat_prn_list_clean = sat_prn_list_clean + prn_line_splited
    try:
        sat_nbr = int(sat_prn_list_clean[0])
    except:
        sat_nbr = 0 

    sat_prn_list_clean = sat_prn_list_clean[1:]

    sat_prn_string = "".join(sat_prn_list_clean)

    sat_prn_list_final = []
    for i in range(sat_nbr):
        sat_prn_list_final.append(sat_prn_string[i*3:i*3+3])

    ### Sigma part
    sat_sig_list_clean = []
    for sig_line in sat_sig_list:
        sig_line_splited = sig_line.split()
        sig_line_splited = [e for e in sig_line_splited if not "+" in e]
        sat_sig_list_clean = sat_sig_list_clean + sig_line_splited

    sat_sig_list_final = [int(e) for e in sat_sig_list_clean[:sat_nbr]]


    ### Export part
    ac_name_list = [ac_name] * sat_nbr
    date_list    = [date] * sat_nbr

    header_df = pd.DataFrame(list(zip(ac_name_list,sat_prn_list_final,
                                      sat_sig_list_final,date_list)),
                             columns=["ac","prn","sigma","epoch"])

    return header_df


def sp3_data_frame_zero_epoch_filter(orb_df_inp):
    """
    Filter an Orbit DataFrame (from a SP3) by removing the null epochs

    Parameters
    ----------
    orb_df_inp : DataFrame
        Orbit DataFrame (from a SP3).

    Returns
    -------
    orb_df_out : DataFrame
        Filtered Orbit DataFrame.

    """

    d_fgrp = orb_df_inp[["epoch", "x", "y", "z"]].groupby("epoch")
    df_sum = d_fgrp.agg(np.sum).sum(axis=1)
    epochs = df_sum[np.isclose(df_sum,0)].index
    
    orb_df_out = orb_df_inp[np.logical_not(orb_df_inp["epoch"].isin(epochs))]
    
    return orb_df_out


def sp3_decimate(file_inp,file_out,step=15):
    """
    Decimate a SP3 file 

    Parameters
    ----------
    file_inp : str
        path of the input SP3 file.
    file_out : str
        path of the output SP3 file.
    step : int, optional
        decimation step in minutes. The default is 15.

    Returns
    -------
    file_out : str
        path of the output SP3 file.
    """


    fin = open(file_inp)

    good_line = True
    outline = []
    n_good_lines = 0 
    for l in fin:
        if l[0] == "*":
            epoc   = conv.tup_or_lis2dt(l[1:].strip().split()) 
            if np.mod(epoc.minute , step) == 0:
                good_line = True
                n_good_lines += 1
            else:
                good_line = False

        if good_line:
            outline.append(l)

    ### replace nb epochs
    line0     = outline[0]
    nlines_orig = outline[0].split()[6]
    nlines_ok = "{:7}".format(n_good_lines).strip().zfill(3)
    line0b = line0.replace(nlines_orig,nlines_ok)
    outline[0] = line0b    


    ### replace step
    line1     = outline[1]
    step_orig = outline[1].split()[3]
    step_ok = "{:14.8f}".format(step * 60).strip()
    line1b = line1.replace(step_orig,step_ok)
    outline[1] = line1b


    with open(file_out,"w+") as Fout:
        for l in outline:
            Fout.write(l)

    return file_out

    
