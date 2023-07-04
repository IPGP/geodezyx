#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 17:12:25 2023

@author: psakicki
"""

import datetime as dt
import subprocess
import urllib
from geodezyx import utils, conv, operational
from xml.etree import ElementTree as et
import os
from threading import Thread
import time

### Project name
project_name="test_calc_PF_03_testJJ"

### Paths where the RINEXs/Products/results are/will be stored
#  Inputs
cfg_files_root_dir="/home/ovsgnss/010_SOFTS/GeodeZYX-Toolbox_v4/geodezyx/000_exemples/groops_frontend/configfiles/040_prototype4/"
rinex_root_path = "/vol/ovpf/miroir_ovpf/DonneesAcquisition/geodesie/GPSData/"
vmf_tropo_root_dir = "/scratch/calcgnss/prods_tropo_vmf3"
prods_gnss_root_dir = "/scratch/calcgnss/prods_gnss/"
# Outputs
log_root_dir="/home/ovsgnss/020_CALC/groops_process/031_groops_frontend_logs"


### RINEX names and dates
specific_sites=["sneg"]
#specific_sites=["knkl"]
specific_sites=["borg","fjag"]

start_epoch=dt.datetime(2023,5,5)
end_epoch=dt.datetime(2023,5,5)



### IGS AC definition 
#igs_ac_10char = "IGS0OPSFIN"
#igs_ac_10char = "IGS2R03FIN"
#igs_ac_10char = "COD0R03FIN"
#igs_ac_10char = "GRG6RE3FIN"
#igs_ac_10char = "COD0MGXFIN"
#igs_ac_10char = "GRG0MGXFIN"
#igs_ac_10char = "GRG0OPSFIN"
#igs_ac_10char = "GFZ0MGXRAP"
#igs_ac_10char = "TUG0R03FIN"

igs_ac_10char = "COD0OPSFIN"

### Config file names for the different steps
cfg_files_dict = dict()

cfg_files_dict["convTropo"] = "010_convTropo_vmf3Daily_v04a.xml"

cfg_files_dict["convProds_clock"] = "021_convProds_clock_v04a.xml"
cfg_files_dict["convProds_bias"]  = "022_convProds_bias_v04a.xml"
cfg_files_dict["convProds_orbit"] = "023_convProds_orbit_v04a.xml"
cfg_files_dict["convProds_alt_attitude"] = "024_convProds_alt_attitude_v04a.xml"
cfg_files_dict["convProds_attitude"] = "025_convProds_attitude_v04a.xml"

cfg_files_dict["convSite_sitelog"] = "030_convSite_sitelog_v04a.xml"
cfg_files_dict["convSite_rnxObs"] = "031_convSite_rnxObs_v04a.xml"

cfg_files_dict["gnssProcessing"] = "040_gnssProcessing_v04a.xml"




Rinexs = operational.rinex_finder(rinex_root_path,
                                  specific_sites=specific_sites,
                                  start_epoch=start_epoch,
                                  end_epoch=end_epoch)

print(Rinexs)


for rnx in Rinexs:
    operational.groops_ppp_full_runner(rnx,
                           project_name,
                           igs_ac_10char,
                           cfg_files_dict,
                           log_root_dir,
                           vmf_tropo_root_dir,
                           prods_gnss_root_dir,
                           cfg_files_root_dir)
