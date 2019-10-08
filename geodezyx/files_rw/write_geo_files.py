#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 18:00:21 2019

@author: psakicki
"""

from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

def write_sndy_light_dat(ts_in,outdir,outprefix):
    """pas fini"""
    fil = open(os.path.join(outdir,outprefix),'w+')
    if isinstance(ts_in,TimeSeriePoint):
        if ts_in.initype() == 'FLH':
            for pt in ts_in.pts:
                lin = ' '.join([str(e) for e in [pt.F , pt.L , pt.H , pt.T , pt.sF , pt.sL , pt.sH ]])
                fil.write(lin + '\n')
    elif isinstance(ts_in,TimeSerieObs):
        if ts_in.typeobs == 'RPY':
            for att in ts_in.obs:
                lin = ' '.join([str(e) for e in [att.R , att.P , att.Y , att.T , att.Q.w , att.Q.x , att.Q.y , att.Q.z ]])
                fil.write(lin + '\n')
    fil.close()
    
    
def write_clk(DFclk_in,clk_file_out,header=""):
    HEAD = header
    Row_str_stk = []
    for irow, row in DFclk_in.iterrows():
        row_str_proto = "{:2} {:4} {:4d} {:2d} {:2d} {:2d} {:2d} {:9.6f} {:2d}   {:19.12e}"
        row_str = row_str_proto.format(row["type"],row["name"],row["year"],
                                       row["month"],row["day"],row["h"],row["minutes"],
                                       row["seconds"],1,row["offset"])
        Row_str_stk.append(row_str)
        
    CORPSE = "\n".join(Row_str_stk)
       
    OUT = HEAD + CORPSE

    with open(clk_file_out,"w+") as Fout:
        Fout.write(OUT)
        Fout.close()
        
    return OUT
