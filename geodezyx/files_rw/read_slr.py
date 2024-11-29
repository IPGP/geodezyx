#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 10:55:04 2021

@author: psakic
"""

import datetime as dt
import io
import os

import pandas as pd

from geodezyx import conv, utils


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
                                              r"\+residuals",
                                              r"\-residuals")
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