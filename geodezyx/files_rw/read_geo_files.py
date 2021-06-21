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


def read_bull_B(file_path_in):
    """
    Read an Bulletin B file

    Parameters
    ----------
    file_path_in :  str
        Path of the file in the local machine.
    Returns
    -------
    DFout : pandas DataFrame
        Returns a panda table format with the data extracted from the file.
    """

    if not utils.is_iterable(file_path_in):
        file_path_in = [file_path_in]

    file_path_in = sorted(file_path_in)

    DFstk = []

    for path_solo in file_path_in:
        S = utils.extract_text_between_elements(path_solo,
                                                "1 - DAILY FINAL VALUES" ,
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
                                DForb.AC,
                                DForb.epoch]).T
    
    DFclk_sp3 = DFclk_sp3.infer_objects()
    
    
    DFclk_sp3.columns = ['type', 'name', 'year', 'month', 
                         'day', 'hour', 'minute', 'second',
                         'n_values', 'bias', 'sigma', 'ac', 'epoch']
    
    return DFclk_sp3


def read_sp3(file_path_in,returns_pandas = True, name = '',
             epoch_as_pd_index = False,km_conv_coef=1,
             skip_null_epoch=True):
    """
    Read a SP3 file (GNSS Orbits standard file) and return X,Y,Z coordinates
    for each satellite and for each epoch

    Parameters
    ----------
    file_path_in : str
        path of the SP3 file
        can handle gzip-compressed file (with .gz/.GZ extension) 

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

    Returns
    -------
    df : Pandas DataFrame
        if returns_pandas == True

    epoch_stk ,  Xstk , Ystk , Zstk , Clkstk : lists
        if returns_pandas == False

    """
    
    AC_name =  os.path.basename(file_path_in)[:3]
    
    if file_path_in[-2:] in ("gz","GZ"):
        F = gzip.open(file_path_in, "r+")
        Lines = [e.decode('utf-8') for e in F]
    else:
        F = open(file_path_in,'r+')
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

    #### read the Header as a 1st check
    Header = read_sp3_header(Lines,AC_name)
    if Header.empty:
        print("WARN:read_sp3: The SP3 looks empty: ",file_path_in)
        if returns_pandas:
            df = pd.DataFrame([], columns=['epoch','sat', 'const',
                                           'sv','type',
                                           'x','y','z',
                                           'clk','AC'])
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
        df = pd.DataFrame(data_stk, columns=['epoch','sat', 'const',
                                             'sv','type',
                                             'x','y','z',
                                             'clk','AC'])
        
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
        print("INFO : return list, very beta : no Sat. Vehicule Number info ...")
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
            date = conv.MJD2dt(int(l.split()[4]))
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
                             columns=["AC","sat","sigma","epoch"])

    return Header_DF


def sp3_DataFrame_zero_epoch_filter(DFsp3):
    """
    Filter an Orbit DataFrame (from a SP3) by removing the null epochs

    Parameters
    ----------
    DFsp3 : DataFrame
        Orbit DataFrame (from a SP3).

    Returns
    -------
    DFsp3_out : DataFrame
        Filtered Orbit DataFrame.

    """
    DFgrp = DFsp3[["epoch","x","y","z"]].groupby("epoch")
    DFsum = DFgrp.agg(np.sum).sum(axis=1)
    Epochs = DFsum[np.isclose(DFsum,0)].index
    
    DFsp3_out = DFsp3[np.logical_not(DFsp3["epoch"].isin(Epochs))]
    
    return DFsp3_out


def read_erp_multi(path_list, 
                   return_array=False,
                   smart_mode=True,
                   ac=None):
    """
    Input :
        path_list : a list of ERP files
        smart_mode : keep only the latest value (True is recommended)
    """
    path_list = sorted(path_list)
    Lstk = []
    for path in path_list:
        L = read_erp(path,ac)
        Lstk.append(L)

    M = np.vstack(Lstk)

    if smart_mode:
        Msmart = np.array([])
        for ilm , lm in enumerate(np.flipud(M)):
            if ilm == 0:
                Msmart = np.array([lm])
            elif lm[0] in Msmart[:,0]:
                continue
            else:
                Msmart = np.vstack((Msmart,lm))

        M = np.flipud(Msmart)

    if return_array:
        return M
    else:
        return pd.DataFrame(M)



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
    
def AC_equiv_vals(AC1,AC2):
    """
    Find common values between 2 Orbit DataFrame (from SP3)

    Parameters
    ----------
    AC1 : DataFrame
        Orbit DataFrame 1.
    AC2 : DataFrame
        Orbit DataFrame 2.

    Returns
    -------
    AC1_ok : DataFrame
        Cleaned Orbit DataFrame 1.
    AC2_ok : DataFrame
        Cleaned Orbit DataFrame 2.
    ACmerged : DataFrame
        Merged orbit DataFrame.

    """
    ### 1) Merge the 2 DF to find common lines
    ACmerged = time_series.merge(AC1 , AC2 , how='inner', on=['epoch', 'sat'])

    ### 2) Extract merged epoch & sv
    common_epoch = ACmerged["epoch"]
    common_sat   = ACmerged["sat"]
    ### 3) Create a boolean line based on common epoch / sv
    common_sat_epoch_AC1 = (AC1["epoch"].isin(common_epoch)) & (AC1["sat"].isin(common_sat))
    common_sat_epoch_AC2 = (AC2["epoch"].isin(common_epoch)) & (AC2["sat"].isin(common_sat))
    ### 4) Get epoch and sv in the combined sol which correspond to the SP3
    AC1new = AC1[common_sat_epoch_AC1].copy()
    AC2new = AC2[common_sat_epoch_AC2].copy()
    ### 5) A sort to compare the same things
    AC1new.sort_values(by=['sat','sv'],inplace=True)
    AC2new.sort_values(by=['sat','sv'],inplace=True)

    ### Check for > 99999 vals
    AC1_bad_bool_9  = (AC1new["x"] > 9999) & (AC1new["y"] > 9999) & (AC1new["z"] > 9999)
    AC1_bad_bool_9  = np.logical_not(np.array(AC1_bad_bool_9))

    AC2_bad_bool_9  = (AC2new["x"] > 9999) & (AC2new["y"] > 9999) & (AC2new["z"] > 9999)
    AC2_bad_bool_9  = np.logical_not(np.array(AC2_bad_bool_9))

    AC12_bad_bool_9 = np.array(np.logical_and(AC1_bad_bool_9 , AC2_bad_bool_9))

    ### Check for NaN vals
    AC1_bad_bool_nan  = np.isnan(AC1new["x"]) & np.isnan(AC1new["y"]) & np.isnan(AC1new["z"])
    AC1_bad_bool_nan  = np.logical_not(np.array(AC1_bad_bool_nan))

    AC2_bad_bool_nan  = np.isnan(AC2new["x"]) & np.isnan(AC2new["y"]) & np.isnan(AC2new["z"])
    AC2_bad_bool_nan  = np.logical_not(np.array(AC2_bad_bool_nan))

    AC12_bad_bool_nan = np.array(np.logical_and(AC1_bad_bool_nan , AC2_bad_bool_nan))

    AC12_bad_bool = np.logical_and(AC12_bad_bool_9, AC12_bad_bool_nan)

    AC1_ok = AC1new[AC12_bad_bool]
    AC2_ok = AC2new[AC12_bad_bool]

    return AC1_ok , AC2_ok , ACmerged


def read_erp(file_path_in,ac=None):
    """
    
    Read IGS Analysis Center ERP files

    Parameters
    ----------
    file_path_in :  str
        Path of the file in the local machine.

    ac :  str
        The analysis center that will be used. 
        If not precised, will be the first 3 letters of the input name


    Returns
    -------
    out1 :  pandas table
        Returns a panda table format with the data extracted from the file.
        

    Note
    ----
    Columns name
    
    ('MJD','X-P (arcsec)', 'Y-P (arcsec)', 'UT1UTC (E-7S)','LOD (E-7S/D)','S-X (E-6" arcsec)','S-Y (E-6" arcsec)',
    'S-UT (E-7S)','S-LD (E-7S/D)','NR (E-6" arcsec)', 'NF (E-6" arcsec)', 'NT (E-6" arcsec)',
    'X-RT (arcsec/D)','Y-RT (arcsec/D)','S-XR (E-6" arcsec/D)','S-YR (E-6" arcsec/D)', 'C-XY', 'C-XT',
    'C-YT', 'DPSI', 'DEPS','S-DP','S-DE')

    """
    
    caminho_arq = file_path_in

    #### FIND DELIVERY DATE
    name = os.path.basename(caminho_arq)

    if not ac:
        ac = name[:3]

    if len(name) == 12:
        dt_delivery = conv.sp3name2dt(caminho_arq)
    elif len(name) == 38:
        dt_delivery = conv.sp3name_v3_2dt(caminho_arq)
    else:
        dt_delivery = conv.posix2dt(0)


    le = open(caminho_arq, 'r')
    letudo = le.readlines()
    le.close()
    tamanho = len(letudo) #usado para saber quantas linhas tem o arquivo

    #para = tamanho #onde o arquivo para de ser lido

    numeros = ['0','1','2','3','4','5','6','7','8','9']
    #le = 0
    #numlin = 0 #numero de linhas da matriz de epocas
    #numcol = 16 #numero de colunas que a matriz final deve ter

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
                    # Doy = (Lines[i][30:33])
                    # Year = (Lines[i][27:29])
                    # Pref_year = '20'
                    # Year = int(Pref_year+Year)
                    # Date = conv.doy2dt(Year,Doy)
                    Date = conv.datestr_sinex_2_dt(Lines[i].split()[5])
                    XPO = float(Lines[i][47:68])*(10**-3)
                    XPO_std = float(Lines[i][69:80])*(10**-3)
                    XPO_stk.append(XPO)
                    XPO_std_stk.append(XPO_std)
                    MJD_stk.append(conv.dt2MJD(Date))
                    #MJD_stk.append(cmg.jd_to_mjd(cmg.date_to_jd(Date.year,Date.month,Date.day)))
                    
                if utils.contains_word(Lines[i],'YPO') and marker:
                    # Doy = (Lines[i][30:33])
                    # Year = (Lines[i][27:29])
                    # Pref_year = '20'
                    # Year = int(Pref_year+Year)
                    # Date = conv.doy2dt(Year,Doy)
                    Date = conv.datestr_sinex_2_dt(Lines[i].split()[5])
                    YPO = float(Lines[i][47:68])*(10**-3)
                    YPO_std = float(Lines[i][69:80])*(10**-3)
                    YPO_stk.append(YPO)
                    YPO_std_stk.append(YPO_std)
                    MJD_stk.append(conv.dt2MJD(Date))
                    #MJD_stk.append(cmg.jd_to_mjd(cmg.date_to_jd(Date.year,Date.month,Date.day)))

                    
                if utils.contains_word(Lines[i],'LOD') and marker:
                    # Doy = (Lines[i][30:33])
                    # Year = (Lines[i][27:29])
                    # Pref_year = '20'
                    # Year = int(Pref_year+Year)
                    # Date = conv.doy2dt(Year,Doy)
                    Date = conv.datestr_sinex_2_dt(Lines[i].split()[5])
                    LOD = float(Lines[i][47:68])*(10**+4)
                    LOD_std = float(Lines[i][69:80])*(10**+4)
                    LOD_stk.append(LOD)
                    LOD_std_stk.append(LOD_std)
                    #MJD_stk.append(cmg.jd_to_mjd(cmg.date_to_jd(Date.year,Date.month,Date.day)))
                    MJD_stk.append(conv.dt2MJD(Date))

        MJD = list(sorted(set(MJD_stk)))
        if len(LOD_stk) == 0:
                LOD_stk = ['0']*len(MJD)
                LOD_std_stk = ['0']*len(MJD)


        for i in range(len(MJD)):

            ERP_data = [ac, MJD[i], XPO_stk[i], YPO_stk[i], 0, LOD_stk[i], XPO_std_stk[i], YPO_std_stk[i],
                     0, LOD_std_stk[i], 0, 0, 0, 0, 0, 0, 0, dt_delivery]

            ERP.append(ERP_data)



    if ac in ('COD','cod','com', 'cof', 'grg', 'mit', 'sio','igs','igr'):
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
                
        linecache.clearcache()
        

#        Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
#                                                 'X-RT','Y-RT','S-XR','S-YR',"Delivery_date"])
#        return Erp_end
#

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
                
        linecache.clearcache()

#        Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
#                                                 'X-RT','Y-RT','S-XR','S-YR'])
#        return Erp_end
#
    if ac in ('gbm', 'gfz','gfr',"p1_","p1r"):
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
        linecache.clearcache()

#        Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
#                                                 'X-RT','Y-RT','S-XR','S-YR'])  ##EH TBM O RATE XY POR DIA??????
#        return Erp_end

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
        linecache.clearcache()



    Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 
                                         'UT1UTC(UT1 -TAI)','LOD',
                                         'S-X','S-Y','S-UT','S-LD',
                                         'NR', 'NF', 'NT',
                                         'X-RT','Y-RT','S-XR','S-YR',
                                         'Delivered_date'])

        
    return Erp_end


### read_erp2 = read_erp
    
    
def read_erp_snx(snx_in):
    """
    Read ERP in SINEX Format

    Parameters
    ----------
    snx_in : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    print("WARN read_erp_snx not IMPLEMENTED !!!!!")
    
    utils.extract_text_between_elements_2(snx_in,
                                          '\+SOLUTION/ESTIMATE',
                                          '-SOLUTION/ESTIMATE',
                                          return_string = True)
    
    

def read_eop_C04(file_path_in):
    """
    read EOP C04 file

    Parameters
    ----------
    file_path_in : str
        path of the EOP C04 file.

    Returns
    -------
    DF : DataFrame
        out DataFrame.

    """
    DF = pd.read_csv(file_path_in,
                     delim_whitespace=True,
                     skip_blank_lines=True,
                     skiprows=13,header=None)
    
    cols = ['yyyy','mm','dd','MJD','x','y','UT1-UTC','LOD','dX','dY',
            'x_Err','y_Err','UT1-UTC_Err','LOD_Err','dX_Err','dY_Err']
    
    DF.columns = cols
    
    DF["epoch"] = conv.MJD2dt(DF.MJD)
    
    return DF

def read_eop_finals(file_path_in,precnut_model = 2000,
                    simplified_EOP_DF="mixed"):
    """
    read EOP finals file


    Parameters
    ----------
    file_path_in : str
        path of the EOP finals file.
    precnut_model : int, optional
        IAU Precession-Nutation Model. Valid values are 1980 and 2000.
        The default is 2000.

    Returns
    -------
    DF : DataFrame
        out DataFrame - Complete content of the finals file.
    DF2 : DataFrame
        out DataFrame - Simplified content of the finals file.
        Based on C04 EOP DataFrame structure.
        If no Bulletin-B values provided, use Bulletin-A.
        
    Note
    ----
    Format description
    https://github.com/marando/phpIERS/blob/master/src/Marando/IERS/default/readme.finals
    """
    
    if precnut_model == 2000:
        x_or_psi = "X"
        y_or_eps = "Y"
    elif precnut_model == 1980:
        x_or_psi = "PSI"
        y_or_eps = "EPS"   
    else:
        raise Exception
    
    cols = ((0,2    ),
    (2,4    ),
    (4,6    ),
    (6,7    ),
    (7,15   ),
    (15,16  ),
    (16,17  ),
    (17,18  ),
    (18,27  ),
    (27,36  ),
    (36,37  ),
    (37,46  ),
    (46,55  ),
    (55,57  ),
    (57,58  ),
    (58,68  ),
    (68,78  ),
    (78,79  ),
    (79,86  ),
    (86,93  ),
    (93,95  ),
    (95,96  ),
    (96,97  ),
    (97,106 ),
    (106,115),
    (115,116),
    (116,125),
    (125,134),
    (134,144),
    (144,154),
    (154,165),
    (165,175),
    (175,185))
        
    col_names = ['yyyy',
    'mm',
    'dd',
    '[blank]',
    'MJD',
    '[blank]',
    'flg_xy',
    '[blank]',
    'x_A',
    'x_A_Err',
    '[blank]',
    'y_A',
    'y_A_Err',
    '[blanks]',
    'flg_UT1-UTC',
    'UT1-UTC_A',
    'UT1-UTC_A_Err',
    '[blank]',
    'LOD_A',
    'LOD_A_Err',
    '[blanks]',
    'flg_nut',
    '[blank]',
    'd' + x_or_psi + '_A',
    'd' + x_or_psi + '_A_Err',
    '[blank]',
    'd' + y_or_eps + '_A',
    'd' + y_or_eps + '_A_Err',
    'x_B',
    'y_B',
    'UT1-UTC_B',
    'd' + x_or_psi + '_B',
    'd' + y_or_eps + '_B']
    
    DF = pd.read_fwf(file_path_in,colspecs=cols,header=None)
    DF.columns = col_names
    ### remove everything with NaN
    DF.dropna(axis=1,how="all",inplace=True)
    ### if not x_A provided, then its a blank line
    DF = DF[np.logical_not(np.isnan(DF['x_A']))]
    
    ### set mjd as int
    DF["MJD"]   = DF["MJD"].astype(np.int64)
    
    #### DF 2 is a simplified version, base on the C04 DF
    DF2 = pd.DataFrame()
    DF2["epoch"] = conv.MJD2dt(DF["MJD"])
    
    if simplified_EOP_DF == 'mixed':
        ### Create epoch and use B-values per default
        DF2[['MJD','x','y','UT1-UTC','LOD','dX','dY']] = DF[['MJD','x_B','y_B','UT1-UTC_B','LOD_A','dX_B','dY_B']]
        
        DF2["Bul"] = "B"
        ### If no B-values, use A-values
        DF2.loc[np.isnan(DF.x_B),'Bul'] = "A"
        DF2.loc[np.isnan(DF.x_B),'x'] = DF.loc[np.isnan(DF.x_B),'x_A']
        DF2.loc[np.isnan(DF.y_B),'y'] = DF.loc[np.isnan(DF.y_B),'y_A']
        DF2.loc[np.isnan(DF['UT1-UTC_B']),'UT1-UTC'] = DF.loc[np.isnan(DF['UT1-UTC_B']),'UT1-UTC_A']
        DF2.loc[np.isnan(DF.dX_B),'dX'] = DF.loc[np.isnan(DF.dX_B),'dX_A']
        DF2.loc[np.isnan(DF.dX_B),'dY'] = DF.loc[np.isnan(DF.dY_B),'dY_A']
                
    elif simplified_EOP_DF == "A":
        DF2[['MJD','x','y','UT1-UTC','LOD','dX','dY']] = DF[['MJD','x_A','y_A','UT1-UTC_A','LOD_A','dX_A','dY_A']]
        DF2["Bul"] = "A"
        DF2.dropna(axis=1,how="all",inplace=True)       

    elif simplified_EOP_DF == "B":
        DF2[['MJD','x','y','UT1-UTC','LOD','dX','dY']] = DF[['MJD','x_B','y_B','UT1-UTC_B','LOD_A','dX_B','dY_B']]
        DF2["Bul"] = "B"
        DF2.dropna(axis=1,how="all",inplace=True)

    else:
        print("ERR: read_eop_C04: check the  simplified_EOP_DF parameter !!!")
        raise Exception
        
    return DF,DF2

########################################################################################################################################
############################################### READ IERS GUS

def read_IERS(file_path_in):
#    file_path_in =  '/dsk/mansur/MGEX_C/ERP_IERS'

    file = open(file_path_in,'r')
    fil = file.readlines()
    file.close()
    return fil


def read_IERS_info(fil, mjd):
    #    mjd = int(mjd)
    X, Y, ut1utc, dx, dy, xerr, yerr, ut1utcerr, dxerr, dyerr, LOD, sig_lod, OMEGA, sig_ome = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
    PoleXY = []

#    mjd = 58230

    for i in range(len(fil)):
    #    linhaatual = linecache.getline(file_path_in, i)
        linhaatual = fil[i]
        if (linhaatual[15:20]) == str(mjd) and (linhaatual[104:105]) != '':
          X = float(linhaatual[23:30])*(10**-3)
          Y = float(linhaatual[31:39])*(10**-3)
          ut1utc = float(linhaatual[40:49])*(10**-3)
          dx = float(linhaatual[51:58])*(10**-3)
          dy = float(linhaatual[59:65])*(10**-3)
          xerr = float(linhaatual[68:74])*(10**-3)
          yerr = float(linhaatual[77:82])*(10**-3)
          ut1utcerr = float(linhaatual[86:93])*(10**-3)
          dxerr = float(linhaatual[94:100])*(10**-3)
          dyerr = float(linhaatual[101:106])*(10**-3)
        elif (linhaatual[15:20]) == str(mjd) and (linhaatual[66:68]) != '':
          LOD = float(linhaatual[22:29])*(10**-3)
          sig_lod = float(linhaatual[30:37])*(10**-3)
          OMEGA = float(linhaatual[38:53])*(10**-3)
          sig_ome = float(linhaatual[56:70])*(10**-3)

    line_data = [mjd,X,Y,ut1utc,dx,dy,xerr,yerr,ut1utcerr,dxerr,dyerr,LOD,sig_lod,OMEGA,sig_ome]
    PoleXY.append(line_data)

    PoleXY_data = pd.DataFrame(PoleXY, columns=['mjd','X','Y','ut1-utc','dx','dy','xerr','yerr','ut1utcerr','dxerr','dyerr','LOD','sig_lod','OMEGA','sig_ome'])
    PoleXY_data.set_index('mjd',inplace=True)
    return PoleXY_data

################################################## END READ IERS GUS

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


### compar_plot
def stations_in_EPOS_sta_coords_file_mono(coords_file_path):
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

    SITE_line_list = utils.grep(coords_file_path , " SITE            m")

    stats_list = []
    mean_mjd_list = []
    for l in SITE_line_list:
        stat = l.split()[8].lower()
        stats_list.append(stat)
        mean_mjd = np.mean([float(l.split()[6]) , float(l.split()[7])])
        mean_mjd_list.append(mean_mjd)

    mjd_final = utils.most_common(mean_mjd_list)
    epoch = conv.MJD2dt(mjd_final)
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
        extract_fct = stations_in_EPOS_sta_coords_file_mono
    else:
        print("ERR : check station file type !!")
        return None

    epoch_list_stk , stats_list_stk = [] , []

    for file_path in files_path_list:
        try:
            epoch , stats_list = extract_fct(file_path)
            epoch_list = [epoch] * len(stats_list)
            epoch_list_stk = epoch_list_stk + epoch_list
            stats_list_stk = stats_list_stk + stats_list
        except:
            print('WARN : something wrong with' , file_path)
            print('TIPS : good files type ? ("sinex" or "EPOS_sta_coords") : ', files_type)

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



def sinex_bench_antenna_DF_2_disconts(DFantenna_in,stat,return_full=False):
    DFantenna_work = DFantenna_in[DFantenna_in["Code"] == stat]
    Start_List     = conv.datestr_sinex_2_dt(DFantenna_work["_Data_Start"])
    End_list       = conv.datestr_sinex_2_dt(DFantenna_work["_Data_End__"])
    if return_full:
        return Start_List,End_list
    else:
        Clean_list = sorted(list(set(Start_List + End_list)))
        Clean_list = [e for e in Clean_list if e != dt.datetime(1970, 1, 1, 0, 0)]
        return Clean_list
    
    
    


 #   _____ ______ ______    ________ _____   ____   _____   ______ _ _           
 #  / ____|  ____|___  /   / /  ____|  __ \ / __ \ / ____| |  ____(_) |          
 # | |  __| |__     / /   / /| |__  | |__) | |  | | (___   | |__   _| | ___  ___ 
 # | | |_ |  __|   / /   / / |  __| |  ___/| |  | |\___ \  |  __| | | |/ _ \/ __|
 # | |__| | |     / /__ / /  | |____| |    | |__| |____) | | |    | | |  __/\__ \
 #  \_____|_|    /_____/_/   |______|_|     \____/|_____/  |_|    |_|_|\___||___/


def read_pdm_res_slr_mono(res_file_in,
                          sol="sol"):
    """
    Read a PDM7 res(idual) file for SLR Validation
    
    Parameters
    ----------
    res_file_in : str
        path of the input res file.
        
    sol : str or lambda fct
        solution name
        if it is a lambda fct, it will grab the sol name from the full residual path
        e.g. : solnam = lambda x: x.split("/")[-5][4:]

    Returns
    -------
    DFout : Pandas DataFrame
        output DataFrame.

    """
    
    dat = conv.sp3name2dt("xxx" + os.path.basename(res_file_in)[:5])
    
    ### get useful values
    L = utils.extract_text_between_elements_2(res_file_in,
                                              "\+residuals",
                                              "\-residuals")
    L = L[3:-1]
    
    output = io.StringIO()
    output.write("".join(L))
    output.seek(0)
    
    ### 
    if utils.is_lambda(sol):
        sol_wrk = sol(res_file_in)
    else:
        sol_wrk = sol

    ### read
    DFout = pd.read_csv(output,header=None,delim_whitespace = True)
    
    ### rename useful columns
    DFout = DFout.rename(columns={0: 'time',
                                  1: 'sta',
                                  2: 'sat',
                                  4: 'res',
                                  5: 'elev',
                                  6: 'azi',
                                  7: 'amb',
                                  8: 'slr',
                                  9: 'dt_sta',
                                  10: 'delay',
                                  11: 'sig_sta'})
    
    DFout["day"] = dat 
    DFout["epoc"] = dat + DFout["time"].apply(lambda x: dt.timedelta(microseconds=x*86400*10*6))
    DFout['sol'] = sol_wrk 
    DFout["sys"] = DFout["sat"].str[0]
    
    return DFout
    
        
def read_pdm_res_slr_multi(Res_file_list_in,sol="sol"):
    """
    Read a PDM7 res(idual) file for SLR Validation
    
    Parameters
    ----------
    Res_file_list_in : list of str
        List of path of the input res files.
        
    sol : str or lambda fct
        solution name
        if it is a lambda fct, it will grab the sol name from the full residual path
        e.g. : solnam = lambda x: x.split("/")[-5][4:]

    Returns
    -------
    DFout : Pandas DataFrame
        output DataFrame
    """
    DFstk = []
    for res_fil in Res_file_list_in:
        DFmono = read_pdm_res_slr_mono(res_fil,sol=sol)
        DFstk.append(DFmono) 
        
    DFmulti = pd.concat(DFstk)
    DFmulti.reset_index(inplace = True,drop=True)

    return DFmulti
   
    
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


 
