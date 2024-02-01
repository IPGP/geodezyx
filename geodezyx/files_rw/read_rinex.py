#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.files_rw contains functions to 
read RINEX files observation files.

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
import datetime as dt
import numpy as np
import pandas as pd
from io import StringIO
import re
import os
import pathlib
from tqdm import tqdm
import hatanaka

#### geodeZYX modules
from geodezyx import operational,utils


#### Import the logger
import logging
log = logging.getLogger(__name__)


def read_rinex2_obs(rnx_in,
                    set_index=None):
    """
    Read a RINEX Observation, version 2 

    Parameters
    ----------
    rnx_in : see below
        input RINEX.
        can be the path of a RINEX file as string or as Path object,
        or directly the RINEX content as a string, bytes, StringIO object or a 
        list of lines
    set_index : str or list of str, optional
        define the columns for the index.
        If None, the output DataFrame is "flat", with integer index
        ["epoch","prn"] for instance set the epoch and the prn as MultiIndex
        The default is None.

    Returns
    -------
    DFrnxobs : Pandas DataFrame / GeodeZYX's RINEX format
    """
    
    #### open block
    try:
        rnx_wrk = hatanaka.decompress(rnx_in)
    except:
        rnx_wrk = rnx_in
        pass

    LINES = utils.open_readlines_smart(rnx_wrk)
    EPOCHS = operational.rinex_read_epoch(rnx_wrk,out_index=True)
    if type(rnx_in) is str or type(rnx_in) is pathlib.Path:
        filename = os.path.basename(rnx_in)
    else:
        filename = "unknown filename"
    
    #### Split header and Observation body
    for il,l in enumerate(LINES):
        if "END OF HEADER" in l:
            i_end_header = il
            break
        
    LINES_header = LINES[:i_end_header+1]
    
    #### get the observables
    Lines_obs = [l for l in LINES_header if '# / TYPES OF OBSERV' in l]
    ## clean SYS / # / OBS TYPES
    Lines_obs = [l[:60] for l in Lines_obs]
    
    ObsAllList_raw = " ".join(Lines_obs).split()
    ObsAllList = ObsAllList_raw[1:]
    ObsAllList = [e for sublist in [(e,e+"_LLI",e+"_SSI") for e in ObsAllList] for e in sublist]
    
    nobs = int(ObsAllList_raw[0])
    nlines_for_obs = int(np.ceil(nobs/5)) ## 5 is the max num of obs in the RIENX specs
    columns_width = nobs*[14,1,1]
    
    DFall_stk = []
    
    #### reading the epochs    
    for iepoc in tqdm(range(len(EPOCHS)),desc="Reading " + filename):
        epoch = EPOCHS[iepoc,0]
        ## define the start/end indices of the epoch block
        iline_start = EPOCHS[iepoc,1] 
        if iepoc == len(EPOCHS)-1:
            iline_end = None
        else:
            iline_end = EPOCHS[iepoc+1,1]
    
        ### get the lines of the epoch block
        Lines_epoc = LINES[iline_start:iline_end]
        
        ### get the satellites for this epoch block
        epoc,nsat,lineconcat,Sats_split,iline_sats_end = _sats_find(Lines_epoc)
        
        ### for each sat, merge the breaked lines
        Lines_obs = Lines_epoc[iline_sats_end+1:iline_end]
        Lines_obs = [e.ljust(80) for e in Lines_obs] # must be exactly 80 char long fut the column trunk !!!!
        Lines_obs = [e.replace("\r","") for e in Lines_obs] # not 100% sure of this one
        Lines_obs = [e.replace("\n","") for e in Lines_obs]
        Lines_obs_merg = [Lines_obs[nlines_for_obs*n:nlines_for_obs*n+nlines_for_obs] for n in range(nsat)]
        Lines_obs_merg = ["".join(e) for e in Lines_obs_merg]
        
        ## read the epoch block using pandas' fixed width reader 
        B = StringIO("\n".join(Lines_obs_merg))
        DFepoch = pd.read_fwf(B,header=None,widths=columns_width)  
        DFepoch.columns = ObsAllList
        DFepoch["prn"] = Sats_split
        DFepoch["sys"] = DFepoch["prn"].str[0]
        DFepoch["epoch"] = epoch
        
        DFall_stk.append(DFepoch)
        
    
    ## final concat and cosmetic (reorder columns, sort)
    DFrnxobs = pd.concat(DFall_stk)
    DFrnxobs = DFrnxobs.reindex(["epoch","sys","prn"] + list(sorted(ObsAllList)),axis=1)
    DFrnxobs.sort_values(["epoch","prn"],inplace=True)
    DFrnxobs.reset_index(drop=True,inplace=True)
    
    
    if set_index:
        DFrnxobs.set_index(set_index,inplace=True)
        DFrnxobs.sort_index(inplace=True)
        
    return DFrnxobs
        

def read_rinex3_obs(rnx_in,
                    set_index=None):
    """
    Read a RINEX Observation, version 3 or 4

    Parameters
    ----------
    rnx_in : see below
        input RINEX.
        can be the path of a RINEX file as string or as Path object,
        or directly the RINEX content as a string, bytes, StringIO object or a 
        list of lines
    set_index : str or list of str, optional
        define the columns for the index.
        If None, the output DataFrame is "flat", with integer index
        ["epoch","prn"] for instance set the epoch and the prn as MultiIndex
        The default is None.

    Returns
    -------
    DFrnxobs : Pandas DataFrame / GeodeZYX's RINEX format
    """
    
    #### open block
    try:
        rnx_wrk = hatanaka.decompress(rnx_in)
    except:
        rnx_wrk = rnx_in
        pass
    
    LINES = utils.open_readlines_smart(rnx_wrk)
    EPOCHS = operational.rinex_read_epoch(rnx_wrk,out_index=True)
    if type(rnx_in) is str or type(rnx_in) is pathlib.Path:
        filename = os.path.basename(rnx_in)
    else:
        filename = "unknown filename"
        
    #### Split header and Observation body
    for il,l in enumerate(LINES):
        if "END OF HEADER" in l:
            i_end_header = il
            break
        
    LINES_header = LINES[:i_end_header+1]
    #LINES_obs = LINES[i_end_header:]
    
    #### get the systems and observations
    Lines_sys = [l for l in LINES_header if 'SYS / # / OBS TYPES' in l]
    ## clean SYS / # / OBS TYPES
    Lines_sys = [l[:60] for l in Lines_sys]
    
    ## manage the 2 lines systems
    for il,l in enumerate(Lines_sys):
        if l[0] == " ":
            Lines_sys[il-1] = Lines_sys[il-1] + l
            Lines_sys.remove(l)
    
    #### store system and observables in a dictionnary
    dict_sys = dict()
    dict_sys_nobs = dict()
    dict_sys_use = dict() # adds the prn, and LLI and SSI indicators
    
    for il,l in enumerate(Lines_sys):
        Sysobs = l.split()
        dict_sys[Sysobs[0]] = Sysobs[2:]
        dict_sys_nobs[Sysobs[0]] = int(Sysobs[1])
        ## adds the LLI and SSI indicators
        dict_sys_use[Sysobs[0]] = [("prn",)] + [(e,e+"_LLI",e+"_SSI") for e in Sysobs[2:]]
        dict_sys_use[Sysobs[0]] = [e for sublist in dict_sys_use[Sysobs[0]] for e in sublist]
        
        if len(Sysobs[2:]) != int(Sysobs[1]):
            print("difference between theorectical and actual obs nbr in the header")  
    
    ## the max number of observable (for the reading)
    ##ObsAllList = list(sorted(set([e for sublist in list(dict_sys_use.values())[1:] for e in sublist]))) 
    nobs_max = max(dict_sys_nobs.values())
    
    DFall_stk = []
    #### reading the epochs
    for iepoc in tqdm(range(len(EPOCHS)),desc="Reading " + filename):
        epoch = EPOCHS[iepoc,0]
        ## define the start/end indices of the epoch block
        iline_start = EPOCHS[iepoc,1] + 1
        if iepoc == len(EPOCHS)-1:
            iline_end = None
        else:
            iline_end = EPOCHS[iepoc+1,1]
        
        Lines_epoc = LINES[iline_start:iline_end]
        ###  Remove CR (Carriage Return) and LF (Line Feed) 
        Lines_epoc = [l.replace('\r', '') for l in Lines_epoc]
        Lines_epoc = [l.replace('\n', '') for l in Lines_epoc]
                
        ## read the epoch block using pandas' fixed width reader 
        B = StringIO("\n".join(Lines_epoc))
        
        columns_width = [3] + nobs_max*[14,1,1]
        DFepoch = pd.read_fwf(B,header=None,widths=columns_width)      
               
        DFepoch_ok_stk = []
        #### assign the correct observable names for each system
        for sys in dict_sys_use.keys():
            #DFepoch_sys = DFepoch[DFepoch[0].str[0] == sys]
            #                   get the sats of the system sys  ||  get the meaningful columns
            DFepoch_sys_clean = DFepoch[DFepoch[0].str[0] == sys].iloc[:,:len(dict_sys_use[sys])]
            DFepoch_sys_clean.columns = dict_sys_use[sys]
            DFepoch_ok_stk.append(DFepoch_sys_clean)
                        
        DFepoch_ok = pd.concat(DFepoch_ok_stk)
        # An epoch column is created to fasten the process
        col_epoch = pd.Series([epoch]*len(DFepoch_ok),name="epoch")
        col_sys = pd.Series(DFepoch_ok['prn'].str[0],name="sys")
        DFepoch_ok = pd.concat([col_epoch,col_sys,DFepoch_ok],axis=1)
                
        DFall_stk.append(DFepoch_ok)
    
    ## final concat and cosmetic (reorder columns, sort)
    DFrnxobs = pd.concat(DFall_stk)
    #Col_names = list(DFrnxobs.columns)
    #Col_names.remove("epoch")
    #Col_names.remove("prn")
    #Col_names = ["epoch","prn"] + sorted(Col_names)
    #DFrnxobs = DFrnxobs.reindex(Col_names,axis=1)
    
    DFrnxobs.reset_index(drop=True,inplace=True)   
    
    if set_index:
        DFrnxobs.set_index(set_index,inplace=True)
        DFrnxobs.sort_index(inplace=True)
        
    return DFrnxobs

############ UTILITY FUNCTIONS


def DFrnx_clean_LLI_SSI(DFrnx_in):
    """
    Remove the Loss of Lock Indicator (LLI) and Signal Strength Indicator (SSI) 
    columns in a DataFrame RINEX
    
    Parameters
    ----------
    DFrnx_in : Pandas DataFrame / GeodeZYX's RINEX format
        A RINEX DataFrame with LLI/SSI columns.

    Returns
    -------
    DFrnx_out : Pandas DataFrame / GeodeZYX's RINEX format
        A RINEX DataFrame without LLI/SSI columns.
    """
    cols = DFrnx_in.columns
    cols_clean = [e for e in cols if not "LLI" in e and not "SSI" in e]
    DFrnx_out = DFrnx_in[cols_clean]
    return DFrnx_out


def observables_dict_per_sys(DFrnx_in):
    """
    Gives the GNSS observables for each GNSS system in a dictionnary


    Parameters
    ----------
    DFrnx_in : Pandas DataFrame / GeodeZYX's RINEX format
        A RINEX DataFrame.

    Returns
    -------
    dict_sys_obs : dict
        A dictionnary with GNSS system as key (G,R,E...).
        And the observalbes for each system as values
        
    Note
    ----
    
    Use dict_sys_obs_clean_LLI_SSI if your want to remove the LLI & SSI values
    """
    
    dict_sys_obs = dict()

    for sys in DFrnx_in["sys"].unique():
        DFsys = DFrnx_in[DFrnx_in["sys"] == sys]
        DFsys_mini = DFrnx_clean_LLI_SSI(DFsys)
        
        ObsSys0 = (DFsys_mini.isna().sum() == len(DFsys)).apply(np.logical_not)
        ObsSys = ObsSys0.index[ObsSys0][3:] ### strating from 3 to clean epoch sys prn
        
        #init_tup = ("epoch","sys","prn") 
        #init_tup = []
        #ObsSys_full = [init_tup] + [(e,e+"_LLI",e+"_SSI") for e in ObsSys] 
        ObsSys_full = [(e,e+"_LLI",e+"_SSI") for e in ObsSys] 
        ObsSys_full = [e for sublist in ObsSys_full for e in sublist]
                    
        dict_sys_obs[sys] = ObsSys_full
        
    return dict_sys_obs
    
def dict_sys_obs_clean_LLI_SSI(dict_sys_obs_in):
    """
    Clean a `dict_sys_obs` (generated by `observables_dict_per_sys`)
    of its LLI and SSI values

    Parameters
    ----------
    dict_sys_obs_in : dict
        A dictionnary with GNSS system as key (G,R,E...).
        And the observalbes for each system as values.

    Returns
    -------
    dict_sys_obs_out : dict
        Same dictionnary cleanned of its LLI and SSI values.

    """
    dict_sys_obs_out = dict()
    
    for sys,obs in dict_sys_obs_in.items():
        obs_clean = [e for e in obs if not "LLI" in e and not "SSI" in e]
        dict_sys_obs_out[sys] = obs_clean
        
    return dict_sys_obs_out

############ INTERNAL FUNCTIONS


def _sats_find(Lines_inp):
    """
    For RINEX2 only
    search for the satellites for each epoch block by reading the 
    EPOCH/SAT record

    Parameters
    ----------
    Lines_inp : List of str
        the lines of one epoch block (EPOCH/SAT + OBSERVATION records).

    Returns
    -------
    bloc_tuple : tuple
        a 5-tuple containing:
            epoc : datetime, the epoch of block
            nsat : int, the number of satellites
            lineconcat : str, the satellites as a concatenated string  
            Sats_split : list of str, the satellites as a list
            il : int, the index of the last line of the EPOCH/SAT record

    """
    iline_bloc = 0
    nlines_bloc = -1
    for il,l in enumerate(Lines_inp):
        ##### 
        re_epoch = '^ {1,2}([0-9]{1,2} * ){5}'
        #re_sat="[A-Z][0-9][0-9]"
        bool_epoch=re.search(re_epoch,l)
        #bool_sat=re.search(re_sat,l)
        
        ### we found an epoch line
        if bool_epoch:
            in_epoch = True
            nsat = int(l[30:32])
            iline_bloc = 0
            nlines_bloc = int(np.ceil(nsat / 12))
            LineBloc = []
            lineconcat = ""
            epoc = operational.read_rnx_epoch_line(l,rnx2=True)
        
        ### we read the sat lines based on the number of sat
        if iline_bloc <= nlines_bloc:
            LineBloc.append(l[32:].strip())
            lineconcat = lineconcat + l[32:].strip()
            iline_bloc += 1
        
        ### we stack everything when the sat block is over
        if iline_bloc == nlines_bloc:
            in_epoch = False
            Sats_split = [lineconcat[3*n:3*n+3] for n in range(nsat)]
            bloc_tuple = (epoc,nsat,lineconcat,Sats_split,il)
            return bloc_tuple



def _line_reader(linein,nobs):
    """
    DISCONTINUED
    
    For RINEX3 Only 
    
    read the content of an observation line. pd.read_fwf does the job
    """
    columns_width = [3] + nobs*[14,1,1]
    columns_cumul = [0] + list(np.cumsum(columns_width))
    
    out_stk = []
    for i in range(len(columns_cumul)-1):
        slic = slice(columns_cumul[i],columns_cumul[i+1])
        out = linein[slic]
        if i == 0:
            out_stk.append(out)
        else:
            try:
                out_stk.append(float(out))
            except:
                out_stk.append(np.nan)
    
    return out_stk
