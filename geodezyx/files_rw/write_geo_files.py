#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 18:00:21 2019

@author: psakicki
"""

#### geodeZYX modules
from geodezyx import conv
from geodezyx import operational
from geodezyx import utils

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names
##########  END IMPORT  ##########




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
    
def write_sp3(SP3_DF_in , outpath):
    """
    Write DOCSTRING
    """
    ################## MAIN DATA
    LinesStk = []


    SP3_DF_in.sort_values(["epoch","sat"],inplace=True)

    EpochList  = SP3_DF_in["epoch"].unique()
    SatList    = sorted(SP3_DF_in["sat"].unique())
    SatListSet = set(SatList)

    for epoc in EpochList:
        SP3epoc   = pd.DataFrame(SP3_DF_in[SP3_DF_in["epoch"] == epoc])
        ## Missing Sat
        MissingSats = SatListSet.difference(set(SP3epoc["sat"]))
        
        for miss_sat in MissingSats:
            miss_line = SP3epoc.iloc[0].copy()
            miss_line["sat"]   = miss_sat
            miss_line["const"] = miss_sat[0]
            miss_line["x"]     = 0.000000
            miss_line["y"]     = 0.000000
            miss_line["z"]     = 0.000000
            miss_line["clk"]   = 999999.999999
            
            SP3epoc = SP3epoc.append(miss_line)

        SP3epoc.sort_values("sat",inplace=True)
        timestamp = conv.dt2sp3_timestamp(conv.numpy_datetime2dt(epoc)) + "\n"

        LinesStk.append(timestamp)

        linefmt = "P{:}{:14.6f}{:14.6f}{:14.6f}{:14.6f}\n"


        for ilin , lin in SP3epoc.iterrows():
            line_out = linefmt.format(lin["sat"],lin["x"],lin["y"],lin["z"],lin["clk"])

            LinesStk.append(line_out)



    ################## HEADER
    ######### SATELLITE LIST

    Satline_stk   = []
    Sigmaline_stk = []

    for i in range(5):
        SatLine = SatList[17*i:17*(i+1)]
        if len(SatLine) < 17:
            complem = " 00" * (17 - len(SatLine))
        else:
            complem = ""

        if i == 0:
            nbsat4line = len(SatList)
        else:
            nbsat4line = ''

        satline = "+  {:3}   ".format(nbsat4line) + "".join(SatLine) + complem + "\n"
        sigmaline = "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"

        Satline_stk.append(satline)
        Sigmaline_stk.append(sigmaline)


    ######### 2 First LINES
    start_dt = conv.numpy_datetime2dt(EpochList.min())

    header_line1 = "#cP" + conv.dt2sp3_timestamp(start_dt,False) + "     {:3}".format(len(EpochList)) + "   u+U IGSXX FIT  XXX\n"

    delta_epoch = int(utils.most_common(np.diff(EpochList) * 10**-9))
    MJD  = conv.dt2MJD(start_dt)
    MJD_int = int(np.floor(MJD))
    MJD_dec = MJD - MJD_int
    gps_wwww , gps_sec = conv.dt2gpstime(start_dt,False,"gps")

    header_line2 = "## {:4} {:15.8f} {:14.8f} {:5} {:15.13f}\n".format(gps_wwww,gps_sec,delta_epoch,MJD_int,MJD_dec)


    ######### HEADER BOTTOM
    header_bottom = """%c G  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
%f  1.2500000  1.025000000  0.00000000000  0.000000000000000
%f  0.0000000  0.000000000  0.00000000000  0.000000000000000
%i    0    0    0    0      0      0      0      0         0
%i    0    0    0    0      0      0      0      0         0
/* PCV:IGSXX_XXXX OL/AL:FESXXXX  NONE     YN CLK:CoN ORB:CoN
/*     GeodeZYX Toolbox Output
/*
/*
"""


    ################## FINAL STACK

    FinalLinesStk = []

    FinalLinesStk.append(header_line1)
    FinalLinesStk.append(header_line2)
    FinalLinesStk = FinalLinesStk + Satline_stk + Sigmaline_stk
    FinalLinesStk.append(header_bottom)
    FinalLinesStk = FinalLinesStk + LinesStk + ["EOF"]

    FinalStr = "".join(FinalLinesStk)

    F = open(outpath,"w+")
    F.write(FinalStr)
    
    
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
