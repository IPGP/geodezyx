#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 22:19:04 2024

@author: psakicki
"""

#### Import the logger
import logging
import os

from geodezyx import operational, utils

log = logging.getLogger('geodezyx')

cfg_template_path = "/home/ovsgnss/.PRIDE_PPPAR_BIN/config_template"
rnx_path = '/scratch/calcgnss/temp_stuffs/2402_tests_PF_pride/SNEG/SNEG00REU_R_20232000000_01D_30S_MO.crx.gz' 
prod_parent_dir = '/scratch/calcgnss/prods_gnss_pride-pppar/prods'  
#prod_ac_name = "GRG0OPSFIN"
#prod_ac_name = "COD0MGXFIN"
prod_ac_name = "WUM0MGXFIN" 
tmp_dir = "/scratch/calcgnss/pride-pppar_process/tmp" 
cfg_dir = "/scratch/calcgnss/pride-pppar_process/cfg" 
run_dir = "/scratch/calcgnss/pride-pppar_process/run" 

#rnx_list = utils.find_recursive('/vol/ovsg/acqui/GPSOVSG/rinex/2024'  ,'*psa1*d.Z*')
#rnx_list = utils.find_recursive('/home/ovsgnss/050_DATA_GNSS/baiededix/OVPF/2024','*borg*d.Z*')
multi_process = 4

########################

cfg_template_path = os.environ['HOME'] + "/.PRIDE_PPPAR_BIN/config_template"
p = '/home/psakicki/aaa_FOURBI/RINEXexemple'
rnx_path_list = utils.find_recursive(p,'*HOUE*gz')

rnx_path_list = ["/home/psakicki/aaa_FOURBI/FOAG00REU_R_20240762000_01H_01S_MO.crx.gz"]
rnx_path_list = ["/home/psakicki/aaa_FOURBI/RINEXexemple/ABMF00GLP_R_20242250000_01D_30S_MO.crx.gz",
                 "/home/psakicki/aaa_FOURBI/RINEXexemple/ABPO00MDG_R_20242250000_01D_30S_MO.crx.gz"]

parent_dir = "/home/psakicki/aaa_FOURBI/pride_pppar_run"

prod_parent_dir = '/home/psakicki/aaa_FOURBI/pride_pppar_prods'  
#prod_ac_name = "GRG0OPSFIN"
#prod_ac_name = "COD0MGXFIN"
prod_ac_name = "WUM0MGXRAP" 
tmp_dir = parent_dir + "/tmp" 
cfg_dir = parent_dir + "/cfg" 
run_dir = parent_dir + "/run" 

mode = "S"

operational.pride_pppar_runner(rnx_path_list, 
                               cfg_template_path,
                               prod_ac_name, 
                               prod_parent_dir, 
                               tmp_dir, cfg_dir, run_dir,
                               mode=mode,
                               force=True)