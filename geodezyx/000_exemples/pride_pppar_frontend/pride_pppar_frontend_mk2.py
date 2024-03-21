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
log = logging.getLogger(__name__)

########################

cfg_template_path = os.environ['HOME'] + "/.PRIDE_PPPAR_BIN/config_template"
parent_dir = os.environ['HOME'] +  "/pride_pppar_run"
prod_parent_dir = os.environ['HOME'] + '/pride_pppar_prods'  

past_1hour = dt.datetime.now() -  dt.timedelta(seconds=3600*1) 
pattern = conv.statname_dt2rinexname_long('XXXX',past_1hour)[9:21]
rnx_parent_dir = '/home/sysop/autorino_workflow_tests/rnx' 
rnx_path_list = utils.find_recursive(rnx_parent_dir, '*' + pattern + '*gz')
print(past_1hour)
print(rnx_path_list) 

prod_ac_name = "WUM0MGXULT" 
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
