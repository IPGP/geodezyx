#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.files_rw contains functions to 
read RINEX files observation files.

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
import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import files_rw


#### Import the logger
import logging
log = logging.getLogger(__name__)


def write_rinex3_body(DFrnx_in):
    """
    Write the body (data records) of a RINEX files (Version 3 only)

    Parameters
    ----------
    DFrnx_in : RINEX DataFrame
        Input RINEX DataFrame.

    Returns
    -------
    rnx_str_out : str
        The RINEX data records as string.
    """
    
    obsdic = files_rw.observables_dict_per_sys(DFrnx_in)
    
    DFepoc = DFrnx_in.groupby("epoch")
    
    lines_stk = []
    for epoc, grpepoc in DFepoc:
        ### store the epoch
        lines_stk.append(files_rw.write_epoch_rinex3(epoc, len(grpepoc)))
        DFsys = grpepoc.groupby('sys')
        
        for sys, grpsys in DFsys:        
            grpsys2 = grpsys[["prn"] + obsdic[sys]].copy()
            l = "{:3s}" + "{:14.3f}{:1d}{:1d}" * int(len(obsdic[sys])/3)
            grpsys2 = grpsys2.apply(lambda x:l.format(*x.values),axis=1)
            grpsys2 = grpsys2.apply(_lli_ssi_0_to_blank)
            
        lines_stk = lines_stk + list(grpsys2)
        
    rnx_str_out = "\n".append(lines_stk)
    
    return rnx_str_out
    


def write_epoch_rinex3(epoch,nsats,epoch_flag=0,clk_offset=0.):
    """
    Write the epoch line in the RINEX3 format

    Parameters
    ----------
    epoch : datetime
        the epoch.
    nsats : int
        numbers of satellites for the epoch.
    epoch_flag : int, optional
        Epoch flag:
        0 : OK
        1 : power failure between previous and current epoch
        >1 : Special event (see documentation)
        The default is 0.
    clk_offset : float, optional
        Receiver clock offset correction (seconds, optional).
        The default is 0..

    Returns
    -------
    lout : str
        output string line.
    """
    
    l = "> {:4d} {:2d} {:2d} {:2d} {:2d} {:11.7f}  {:1d}{:3d}      {:15.12f}"
    
    lout = l.format(epoch.year,
                    epoch.month,
                    epoch.day,
                    epoch.hour,
                    epoch.minute,
                    epoch.second + epoch.microsecond * 10**-6,
                    epoch_flag,
                    nsats,
                    clk_offset)
    
    return lout

#### INTERNAL FUNCTION

def _lli_ssi_nan_to_0(DFrnx_in):
    """
    Change the NaN of the SSI and LLI columns for the conversion to string
    Take a RINEX DataFrame as input and 
    Returns a modified RINEX DataFrame as output
    """
    
    for colname in DFrnx_in:
        if colname.endswith("LLI") or colname.endswith("SSI"):
            DFrnx_in[colname] = DFrnx_in[colname].replace(np.nan,0)
            DFrnx_in[colname] = DFrnx_in[colname].astype(int)
            
    return DFrnx_in


def _lli_ssi_0_to_blank(rnxlinestr_in):
    """
    Remove 0 values of the SSI and LLI and replace them with blank (" ") ones.
    Take a RINEX-formated string line as input 
    Returns a modified RINEX-formated string line as output

    can be used with DF.apply where DF is here a DataFrame of 
    RINEX-formated string lines
    """
    llen = len(rnxlinestr_in)
    rnxlinelist = list(rnxlinestr_in)
    ### positions of the LLI and the SSI
    lli_ssi_pos = list(np.arange(17,llen,14+2)) + list(np.arange(18,llen,14+2))
    
    for i in lli_ssi_pos:
        if rnxlinelist[i] == "0":
            rnxlinelist[i] = " "
            
    return("".join(rnxlinelist))
