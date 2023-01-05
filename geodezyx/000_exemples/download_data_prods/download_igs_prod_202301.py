#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 18:06:21 2023

@author: psakicki
"""

#### Import star style

from geodezyx import operational

#####################################################################
##### define the wanted products i.e. the file extension
# upper/lower case does not matter, it will be adjusted internally
# depending on the short/long naming
Prod_types  = ("sp3","clk","bia","obx")

#####################################################################
##### define the wanted analysis centers with its 3 char. abreviation
# upper/lower case does not matter, it will be adjusted internally
# depending on the short/long naming
AC_names    = ("IGS")

#####################################################################
##### define the wanted period
# start and end are python datetimes
# start_end_date_easy is a front end function when you want 
# to provide year/doy

start,end = operational.start_end_date_easy(2000, 1,
                                            2021, 365)

#####################################################################
##### define the local destination folder of the downloaded products  
archive_dir = '/scratch/calcgnss/prods_gnss/igs_repro3'

#####################################################################
##### define the IGS data center
# (if needed for  advenced users, the full DC list is within the function)
# we highly recommend to use only the cddis
data_center = 'cddis'

#####################################################################
##### choose the reprocessing campaign
#  3 for repro3 etc...
# set 0 for the routine products
repro = 3

#####################################################################
##### choose Multi-GNSS (MGEX) or "legacy" GPS-Only products
mgex = False


#### create a shorter alias name for the download function
download_fct = operational.multi_downloader_orbs_clks_2

Prods_out_list = download_fct(archive_dir,
                              start,end,
                              AC_names,
                              Prod_types,
                              archtype ='week',
                              archive_center=data_center,
                              parallel_download=1,
                              repro=repro,
                              mgex=mgex)
