#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 17:22:21 2019

@author: psakicki
"""

from megalib import *
import geoclass as gcls
from io import BytesIO,StringIO


###### Define export paths
### Exported plot folder
plots_path = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/1806_Jordan/plots_timeseries"
### Exported ASCII coordinate folders
export_path= "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/1806_Jordan/RAW_NEU"
# Make sure the export dir is created
gf.create_dir(plots_path)
gf.create_dir(export_path)

### Exported simple latitude Longitude folder
export_latlonfile_path= "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/1806_Jordan/RAW_NEU/Jordan_latlon.txt"

###### Read a simple coordinate file
### For instance, a coordinate file from the Nevada Geodetic Laboratory (Blewitt et al.,2018).
### Let's take ZIMM, which can be downloaded here
### http://geodesy.unr.edu/gps_timeseries/tenv3/IGS08/ZIMM.IGS08.txyz2
path_coords = "/home/psakicki/Downloads/ZIMM.IGS08.txyz2"
TS = gcls.read_nevada(path_coords,"xyz")

#So far, Can manage several software with the following functions :
#read_rtklib(filein)
#read_gipsy_apps(filein)
#read_tdp(filein)
#read_track(filein)
#read_sonardyne_posi(filein)
#read_sonardyne_attitude(filein)
#read_gins(filein,'kine')
#read_gins_solution(filein)
#read_gins_solution_multi()
#read_gins_double_diff(filein)
#read_gins_double_diff_multi(filelistin)
#read_qinsy(filein,2014,0o4,0o4)
#read_IGS_coords(filein,initype='auto')
#read_nrcan_csv(filein , associated_ps_file = '', statname = '')
#read_hector_neu()
#read_bull_B()
#read_nrcan_pos()


### We can store TSs in a list
TSlist = []
TSlist.append(TS)


#### Determine mean position, and then transfert coordinate into a 
#### topocentric frame centered on this mean position
TS.ENUcalc_from_mean_posi()

#### get this mean position
TS.mean_posi()

### Access to the main attributes

# get the values as list
X,Y,Z,T,sX,sY,sZ = TS.to_list("XYZ") # Cartesian
F,L,H,T,sF,sL,sH = TS.to_list("FLH") # Geographic
E,N,U,T,sE,sN,sU = TS.to_list("ENU") # Topocentric
# get the station name
TS.stat
# get the list of values as point objects
TS.pts
# get a random point object
pt = TS.aleapt()
# get values of a Point object 
pt.X,pt.Y,pt.Z # Cartesian
pt.F,pt.L,pt.H  # Geographic
pt.E,pt.N,pt.U, # Topocentric
pt.T,pt.Tdt # T = time as POSIX, Tdt = time as Datetime


### Add a point manually and add it to the TimeSerie
#point = Point(X,Y,Z,T,'XYZ',sX,sY,sZ)
#tsout.add_point(point)



################################################
### Read discontinuities from the ITRF file for instance
# ftp://itrf.ign.fr/pub/itrf/itrf2014/ITRF2014-soln-gnss.snx
p_discont="/home/psakicki/GFZ_WORK/GGRSP_TESTs/1902_Tests_Quality_GTRF/ITRF2014-soln-gnss.snx"

# Trick to remove the comments after " -"
discont_clean_str = "\n".join([e.split(" -")[0] for e in open(p_discont,"r").readlines()])
DFantenna = gfc.read_sinex_versatile(StringIO(discont_clean_str),
                                     "SOLUTION/DISCONTINUITY",convert_date_2_dt = True)

################################################



for ts in TSlist:

    print("#### Working with", ts.stat)
          
    #### Security if the TS contains only one point
    if ts.nbpts == 1:
        print(ts.stat,"excluded bc ts.nbpts == 1")
        continue

    if 1:
        #### Exclude points with a too high sigma
        std_val_lim = 0.4

        ts2 = gcls.std_dev_cleaner(ts,std_val_lim,"ENU",verbose=False)
        if ts2.nbpts == 0:
            print(ts.stat,"excluded bc ts2.nbpts == 0")
            continue
    
    if 1:
        #### Exclude outliers based on the MAD strategy
        ts3 = gcls.mad_cleaner(ts2,coortype="ENU",
                               detrend_first=True,verbose=True)


    if 1: ### Apply a Time Window on the TS, here we strat in 1997
        ts3 = gcls.time_win(ts3,[(dt.datetime(1997,1,1),
                                  dt.datetime(2099,12,1))],
                                 mode='keep')

    if 1: ### Add discontinuities to the TS
        # Find the discontinuities for the station
        Discont_stat = DFantenna[DFantenna[0] == ts.stat][5]
        Discont_stat = Discont_stat[Discont_stat > dt.datetime(1980,1,1)]
        # Add them in the Time Serie
        ts3.set_discont(Discont_stat)

    if 1: ### PLOT THE TS
        Fig = plt.figure()
        ts3.plot("ENU",fig=Fig)
        ts3.plot_discont()
        gcls.export_ts_plot(ts3,plots_path,close_fig_after_export=False)
     
   
    if 0: # export the TS as a HECTOR/MIDAS ENU file
        gcls.export_ts_as_neu(ts3, export_path , "")
        gcls.export_ts_as_midas_tenu(ts3, export_path , "")