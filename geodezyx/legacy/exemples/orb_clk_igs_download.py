#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 14:48:56 2018

@author: psakicki
"""

from megalib import *
import softs_runner
import itertools
import geo_files_converter_lib

#### Destination of the products files
orbit_mgex_dir = "/home/USER/ARCHIVE_PROD_IGS-REGULAR/"

#### Products files type (see doc of fct multi_downloader_orbs_clks())
files_type_list  = ["sp3","clk","erp","bia","snx","ssc","sum"]

#### 1st STEP : Selection of the ACs


######### IGS-REFFRAME
#repro=0
#calc_center_list = ["igsYYP"]


######### REPRO 2
#repro=2
#calc_center_list = ["co2","cf2","em2","es2","gf2","gr2","ig2","mi2","ng2","si2","jp2","ul2"]

######### REGULAR
repro=0
calc_center_list = ["cod","cof","emr","esa","gfz","grg","igr","igs","jpl","mit","ngs","sio"]

#### 2nd STEP : Selects where ACs products are coming from
# in a dictonnary : select the good server as the dict key, and the goods ACs
# in an associated list
# (see doc of fct multi_downloader_orbs_clks())
# IMPORTANT : ALL ACs in calc_center_list must be in the following archive_dict
archive_dict = dict()
archive_dict["cddis"] = calc_center_list

bool_uncompress      = True  # uncompression of the files
bool_long2short_name = False # convert long name to short name

#### 3rd STEP : Time range selection
# start, end are datetime

week_begin = 1770
week_end   = 2030

start = geok.gpstime2dt(week_begin,0)
end   = geok.gpstime2dt(week_end,6)

parallel_download = 4

############### END OF THE PARAMETERS    #######################
############### BEGINING OF THE DOWNLOAD #######################

iter_list = itertools.product(files_type_list , calc_center_list)

for fil , clc in iter_list:

    # Find the corresponding Archive Center for the wished Analysis Center
    archive_center = genefun.dic_key_for_vals_list_finder(archive_dict , clc)

    if not archive_center:
        print("ERR : AC " , clc , "is not in the archive dictionary !!!" )
        print("      Define an archive for this AC in the archive dictionary")
        raise Exception


    # DOWNLOAD
    orblis = softs_runner.multi_downloader_orbs_clks(orbit_mgex_dir , start , end ,
                                                     sp3clk=fil,parallel_download=parallel_download,
                                                     archtype='wkwwww',
                                                     calc_center = clc,
                                                     archive_center = archive_center,
                                                     repro=repro,
                                                     force_weekly_file=False)


    # Post-download actions : uncompression and renaming
    if bool_uncompress:
        orblis = [geo_files_converter_lib.unzip_gz_Z(e) for e in orblis]

    if bool_long2short_name and (clc in archive_dict["cddis_mgex_longname"]):
        orblis = [softs_runner.orbclk_long2short_name(e,center_id_last_letter="m",rm_longname_file=False) for e in orblis]
