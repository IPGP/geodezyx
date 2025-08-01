#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 20:09:27 2023

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules

P="/home/psakicki/GFZ_WORK/IPGP_WORK/REVOSIMA/0000_Pressure_Mayotte/010_from_Treden/RawData/transfer_3167556_files_8a99b7f8/204657_20210409_1130_data.txt"

DF = pd.read_table(P,sep=",")
#DF = DF.infer_objects()


#pd.to_datetime(DF.Time)


import pyrsktools

prsk = "/home/psakicki/GFZ_WORK/IPGP_WORK/REVOSIMA/0000_Pressure_Mayotte/011_RawData/010_A0A_RBR/204657_20210919_1458.rsk"
rsk = pyrsktools.RSK(prsk, readHiddenChannels=True)
rsk.open()



if 0:
    rsk._reader._query("data")
    
    t1 = np.datetime64("2021-03-03")
    t2 = np.datetime64("2021-06-04")
    
    rsk.readdata(t1,t2)


# DF = pd.DataFrame(rsk.data)


# outputDir="/home/psakicki/aaa_FOURBI/exportRSK"
# rsk.RSK2CSV(outputDir=outputDir)