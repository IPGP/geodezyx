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

log = logging.getLogger('geodezyx')

 #   _____ _            _       __ _ _           
 #  / ____| |          | |     / _(_) |          
 # | |    | | ___   ___| | __ | |_ _| | ___  ___ 
 # | |    | |/ _ \ / __| |/ / |  _| | |/ _ \/ __|
 # | |____| | (_) | (__|   <  | | | | |  __/\__ \
 #  \_____|_|\___/ \___|_|\_\ |_| |_|_|\___||___/
                                               
                                               

def read_clk(file_path_in,names_4char=False):
    """
    Read an IGS clk file

    Parameters
    ----------
    file_path_in :  str
        Path of the file in the local machine.
        can handle gzip-compressed file (with .gz/.GZ extension) 
    names_4char : bool
        Force the station names to have 4 charaters
        (new naming convention is longer)

    Returns
    -------
    DFclk : pandas DataFrame
        Returns a panda table format with the data extracted from the file.
        
    Note
    ----
    Bias is given in seconds
    """
    HeadLine = utils.grep(file_path_in,"END OF HEADER",
                          only_first_occur=True,line_number=True)
    
    DFclk = pd.read_csv(file_path_in,skiprows=HeadLine[0]+1,header=None,
                        delim_whitespace = True,
                        names=['type', 'name', 'year', 'month', 'day', 'hour',
                               'minute', 'second',"n_values",'bias', 'sigma'])
    
    
    # Special case when D-0n instead of E-0n (e.g. IAC), 
    # then values are not converted to float and kept as generic objects...
    if DFclk['bias'].dtype == "O" and DFclk['bias'].str.match(r".*D(\+|-)\d\d$").any():
        DFclk['bias'] = DFclk['bias'].str.replace("D","E").astype(float)
        
    # remove EOF line    
    DFclk = DFclk[DFclk["type"] != "EOF"]
    
    # convert date to int
    #for d in ['year', 'month', 'day', 'hour','minute']:
    #    DFclk[d] = DFclk[d].astype(int)
        
    DFclk["year"] = DFclk["year"].astype(int)
    
    DFclk["ac"] = os.path.basename(file_path_in)[:3] 
    DFclk["name"] = DFclk["name"].str.upper()
    if names_4char:
        DFclk['name'] = DFclk['name'].str[:4]
            
    DFclk['epoch'] = pd.to_datetime(DFclk[['year', 'month', 'day',
                                           'hour','minute', 'second']])
    DFclk.path = file_path_in
    
    return DFclk


def read_clk_from_sp3(file_path_or_DForb_in):
    """
    Get the clock values from a SP3 file of an Orbit DataFrame,
    formated as a Clock DataFrame

    Parameters
    ----------
    file_path_or_DForb_in : str or DataFrame
        the input SP3 file path or an Orbit DataFrame.

    Returns
    -------
    DFclk_sp3 : DataFrame
        Clock DataFrame.

    """
    if type(file_path_or_DForb_in) is str:
        DForb = read_sp3(file_path_or_DForb_in)
    else:
        DForb = file_path_or_DForb_in
        
    nlines = len(DForb)

    ### replace 999999.999999 with NaN
    DForb.loc[np.isclose(DForb.clk,999999.999999),'clk'] = np.nan

    DFclk_sp3 = pd.DataFrame([["AS"] * nlines,
                                DForb.sat,
                                DForb.epoch.dt.year,
                                DForb.epoch.dt.month,
                                DForb.epoch.dt.day,
                                DForb.epoch.dt.hour,
                                DForb.epoch.dt.minute,
                                DForb.epoch.dt.second,
                                [1] * nlines,
                                DForb.clk * 10**-6,
                                [np.nan] * nlines,
                                DForb.ac,
                                DForb.epoch]).T
    
    DFclk_sp3 = DFclk_sp3.infer_objects()
    
    
    DFclk_sp3.columns = ['type', 'name', 'year', 'month', 
                         'day', 'hour', 'minute', 'second',
                         'n_values', 'bias', 'sigma', 'ac', 'epoch']
    
    return DFclk_sp3


def clk_decimate(file_in,file_out,step=300):
    """
    Decimate a .clk file 

    Parameters
    ----------
    file_in : str
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


    Fin = open(file_in)

    good_line = True
    outline = []

    for l in Fin:
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
    file_1 : str
        file 1 path.
    file_2 : str
        file 2 path.

    Returns
    -------
    DFdiff : DataFrame
        clock differences epoch byx epoch between the two files.

    """
    if type(file_1) is str:
        DF1 = read_clk(file_1)
    else:
        DF1 = file_1

    if type(file_2) is str:
        DF2 = read_clk(file_2)
    else:
        DF2 = file_2

    name_list = sorted(set(DF1['name']).intersection(set(DF2['name'])))

    epoc_common_stk = []
    name_list_stk   = []
    type_list_stk   = []
    dClk_stk        = []

    for nam in name_list:
        DF1nam_bool = (DF1['name'] == nam)
        DF2nam_bool = (DF2['name'] == nam)

        DF1nam = DF1[DF1nam_bool]
        DF2nam = DF2[DF2nam_bool]

        epoc_common = set(DF1nam['epoc']).intersection(set(DF2nam['epoc']))

        DF1a = DF1nam[ DF1nam['epoc'].isin(epoc_common) ]
        DF2a = DF2nam[ DF2nam['epoc'].isin(epoc_common) ]

        dClk = np.array(DF1a['clk_bias']) - np.array(DF2a['clk_bias'])

        type_list_stk   += list(DF1a['type'])
        name_list_stk   += [nam] * len(dClk)
        epoc_common_stk += list(epoc_common)
        dClk_stk        += list(dClk)

    DFdiff = pd.DataFrame((list(zip(name_list_stk,type_list_stk,
                                   epoc_common_stk,dClk_stk))),
                              columns = ['name','type','date','diff'])
    return DFdiff

 #   ____       _     _ _       _______ _____ ____     __ _ _           
 #  / __ \     | |   (_) |     / / ____|  __ \___ \   / _(_) |          
 # | |  | |_ __| |__  _| |_   / / (___ | |__) |__) | | |_ _| | ___  ___ 
 # | |  | | '__| '_ \| | __| / / \___ \|  ___/|__ <  |  _| | |/ _ \/ __|
 # | |__| | |  | |_) | | |_ / /  ____) | |    ___) | | | | | |  __/\__ \
 #  \____/|_|  |_.__/|_|\__/_/  |_____/|_|   |____/  |_| |_|_|\___||___/
                                                                      
def read_sp3(file_path_in,returns_pandas = True, name = '',
             epoch_as_pd_index = False,km_conv_coef=1,
             skip_null_epoch=True,
             new_col_names=True,
             pos_vel_col=False):
    """
    Read a SP3 file (GNSS Orbits standard file) and return X,Y,Z coordinates
    for each satellite and for each epoch

    Parameters
    ----------
    file_path_in : str
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
    
    AC_name =  os.path.basename(file_path_in)[:3]
    
    if file_path_in[-2:] in ("gz","GZ"):
        F = gzip.open(file_path_in, "r+")
        Lines = [e.decode('utf-8') for e in F]
    elif file_path_in[-2:] in (".Z",):
        import ncompress
        fh = open(file_path_in, 'rb')
        F = ncompress.decompress(fh).decode("utf-8") 
        Lines = F.split('\n')
    else:
        F = open(file_path_in,'r')
        Lines = F.readlines()

    header = True

    #### List/DF initialization
    epoch_stk = []
    Xstk , Ystk , Zstk , Clkstk = [],[],[],[]
    #Typestk     = []
    data_stk    = []
    AC_name_stk = []
   
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
    Header = read_sp3_header(Lines,AC_name)
    if Header.empty:
        log.warning("The SP3 looks empty: ",file_path_in)
        if returns_pandas:            
            df = pd.DataFrame([], columns=col_names)
            return df
        else:
            return  epoch_stk ,  Xstk , Ystk , Zstk , Clkstk , AC_name_stk

    for l in Lines:
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
                X = np.nan
            else:
                X   = float(l[4:18]) * km_conv_coef               
            if '*' in l[18:32] or not (l[18:32] and l[18:32].strip()):
                Y = np.nan
            else:
                Y   = float(l[18:32])* km_conv_coef                
            if '*' in l[32:46] or not (l[32:46] and l[32:46].strip()):
                Z = np.nan
            else:
                Z   = float(l[32:46])* km_conv_coef
            if '*' in l[46:60] or not (l[46:60] and l[46:60].strip()):
                Clk = np.nan
            else:
                Clk = float(l[46:60])
            
            typ = l[0]

            if returns_pandas:
                line_data = [epoc,sat_sat,sat_nat,sat_sv,typ,X,Y,Z,Clk,AC_name]
                data_stk.append(line_data)
            else:
                epoch_stk.append(epoc)
                Xstk.append(X)
                Ystk.append(Y)
                Zstk.append(Z)
                Clkstk.append(Clk)


    AC_name_stk = [AC_name] * len(Xstk)

    if returns_pandas:
        df = pd.DataFrame(data_stk, columns=col_names)
        
        if skip_null_epoch:
            df = sp3_DataFrame_zero_epoch_filter(df)

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
        log.info("return list, very beta : no Sat. Vehicule Number info ...")
        return  epoch_stk ,  Xstk , Ystk , Zstk , Clkstk , AC_name_stk

def read_sp3_header(sp3_in,ac_name_in=None):
    """
    Read a SP3 file header and return a Pandas DataFrame
    with sat. PRNs and sigmas contained in the header

    Parameters
    ----------
    sp3_in : str or list
        path of the SP3 file
        can handle gzip-compressed file (with .gz/.GZ extension) 

        can also handle the sp3 content as a list of strings
        (useful when read_sp3_header is used as a subfunction of read_sp3)
        
    ac_name_in : str
        force the AC name
        (necessary when read_sp3_header is used as a subfunction of read_sp3)
        
    Returns
    -------
    Header_DF : Pandas DataFrame
        2 columns "sat", "sigma"

    Note
    -------
    More infos about the sigma
    http://acc.igs.org/orbacc.txt
    """

    if type(sp3_in) is list: 
        ### case when read_sp3_header is used as a subfunction of read_sp3
        Lines = sp3_in
    elif sp3_in[-2:] in ("gz","GZ"):
        F = gzip.open(sp3_in, "r+")
        Lines = [e.decode('utf-8') for e in F]
    else:
        F = open(sp3_in,'r+')
        Lines = F.readlines()
    
    if not ac_name_in:
        ac_name = os.path.basename(sp3_in)[:3]
    else:
        ac_name = ac_name_in

    Sat_prn_list = []
    Sat_sig_list = []

    for il , l in enumerate(Lines):
        if il == 1:
            date = conv.mjd2dt(int(l.split()[4]))
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
    try:
        sat_nbr = int(Sat_prn_list_clean[0])
    except:
        sat_nbr = 0 

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
                             columns=["ac","prn","sigma","epoch"])

    return Header_DF


def sp3_DataFrame_zero_epoch_filter(DFsp3):
    '''
    Filter an Orbit DataFrame (from a SP3) by removing the null epochs

    Parameters
    ----------
    DFsp3 : DataFrame
        Orbit DataFrame (from a SP3).

    Returns
    -------
    DFsp3_out : DataFrame
        Filtered Orbit DataFrame.

    '''

    DFgrp = DFsp3[["epoch","x","y","z"]].groupby("epoch")
    DFsum = DFgrp.agg(np.sum).sum(axis=1)
    Epochs = DFsum[np.isclose(DFsum,0)].index
    
    DFsp3_out = DFsp3[np.logical_not(DFsp3["epoch"].isin(Epochs))]
    
    return DFsp3_out


def sp3_decimate(file_in,file_out,step=15):
    """
    Decimate a SP3 file 

    Parameters
    ----------
    file_in : str
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


    Fin = open(file_in)

    good_line = True
    outline = []
    n_good_lines = 0 
    for l in Fin:
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

    
