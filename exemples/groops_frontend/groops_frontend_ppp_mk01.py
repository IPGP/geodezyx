#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 17:12:25 2023

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules

rinex_path = "/vol/ovpf/miroir_ovpf/DonneesAcquisition/geodesie/GPSData/2023/008/fjag0080.23d.Z"
rinex_root_path = "/vol/ovpf/miroir_ovpf/DonneesAcquisition/geodesie/GPSData/"

igs_ac_10char = "IGS0OPSFIN"
igs_ac_10char = "IGS2R03FIN"
igs_ac_10char = "COD0R03FIN"
igs_ac_10char = "GRG6RE3FIN"
igs_ac_10char = "COD0MGXFIN"
igs_ac_10char = "GRG0MGXFIN"
igs_ac_10char = "GRG0OPSFIN"
igs_ac_10char = "GFZ0MGXRAP"
igs_ac_10char = "COD0OPSFIN"
igs_ac_10char = "TUG0R03FIN"


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

#specific_sites=["borg","fjag"]
specific_sites=["sneg"]
start_epoch=dt.datetime(2019,1,1)
end_epoch=dt.datetime(2021,12,31)

Rinexs = operational.rinex_finder(rinex_root_path,
                                  specific_sites=specific_sites,
                                  start_epoch=start_epoch,
                                  end_epoch=end_epoch)

print(Rinexs)
project_name="test_calc_PF_03"
vmf_tropo_root_dir = "/scratch/calcgnss/prods_tropo_vmf3"
prods_gnss_root_dir = "/scratch/calcgnss/prods_gnss/"
cfg_files_root_dir="/opt/softs_gnss/groops/config/040_prototype4/"
log_root_dir="/home/ovsgnss/020_CALC/groops_process/031_groops_frontend_logs"

for rnx in Rinexs:
    groops_ppp_full_runner(rnx,
                           project_name,
