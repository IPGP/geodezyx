#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 22:19:04 2024

@author: psakicki
"""

from geodezyx import operational, utils, conv
import datetime as dt
import numpy as np
import os
import multiprocessing as mp

#### Import the logger
import logging
log = logging.getLogger('geodezyx')

########################

cfg_template_path = os.environ['HOME'] + "/.PRIDE_PPPAR_BIN/config_template"
parent_dir = os.environ['HOME'] +  "/034_PROCESS_pride-pppar/2410_gwada2410a"
prod_parent_dir = os.environ['HOME'] + '/060_PRODS_GNSS/pride_pppar_prods'  

hours_ago = 2
past_epoch = dt.datetime.utcnow() -  dt.timedelta(seconds=3600*hours_ago) 
pattern = conv.statname_dt2rinexname_long('XXXX',past_epoch)[9:21]
rnx_parent_dir = '/home/ovsgnss/050_DATA_GNSS/data_ovs_glass/GL' 
rnx_path_list = utils.find_recursive(rnx_parent_dir, '*' + '_2023104*' + '*gz') 
print(past_epoch)
print(rnx_path_list) 

prod_ac_name = "WUM0MGXFIN" 
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


