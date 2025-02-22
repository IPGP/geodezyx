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
parent_dir = os.environ['HOME'] +  "/034_PROCESS_pride-pppar/2410_gwada2410b"
prod_parent_dir = os.environ['HOME'] + '/060_PRODS_GNSS/prods_gnss_pride-pppar/prods'  

prod_ac_name = "COD0MGXFIN" 
prod_ac_name = "WUM0MGX(FIN|RAP)" 
tmp_dir = parent_dir + "/tmp" 
cfg_dir = parent_dir + "/cfg" 
run_dir = parent_dir + "/run" 
mode = "S"
multi_process = 12

site_list=['ABG0',
        'AGAL',
        'AMC0',
        'BULG',
        'CBEZ',
        'CRA2',
        'F562',
        'F8O2',
        'FNG0',
        'HOUZ',
        'LEN0',
        'MAD0',
        'PAR1',
        'PSA1',
        'SOUF',
        'STG0',
        'TAR1']

hours_ago = 2
past_epoch = dt.datetime.utcnow() -  dt.timedelta(seconds=3600*hours_ago) 
pattern = conv.statname_dt2rinexname_long('XXXX',past_epoch)[9:21]

for site in site_list:
    rnx_parent_dir = '/home/ovsgnss/050_DATA_GNSS/data_ovs_glass/GL' 
    rnx_path_list = utils.find_recursive(rnx_parent_dir, '*' + site + '*_2023*' + '*gz') 
    print(past_epoch)
    print(rnx_path_list) 


    operational.pride_pppar_runner(rnx_path_list, 
                                   cfg_template_path,
                                   prod_ac_name, 
                                   prod_parent_dir, 
                                   tmp_dir, cfg_dir, run_dir,
                                   mode=mode,
                                   multi_process = multi_process,
                                   dl_prods=0,
                                   dl_prods_only=0,
                                   force=False)


