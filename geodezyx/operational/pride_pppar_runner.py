#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 18:01:49 2024

@author: psakic
"""

#### Import star style
from geodezyx import operational
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import utils

#import datetime as dt
import os 
from subprocess import Popen, PIPE

def pride_pppar_runner(rnx_path,
                       cfg_template_path,
                       prod_ac_name,
                       prod_parent_dir,
                       tmp_dir,
                       cfg_dir,
                       run_dir,
                       cfg_prefix='pride_pppar_cfg_1a',
                       bin_dir=None):
    
    if not bin_dir:
        bin_dir = os.path.join(os.environ['HOME'],'.PRIDE_PPPAR_BIN')

    srt,end = operational.rinex_start_end(rnx_path)
    doy,year = conv.dt2doy_year(srt)
    
    tmp_dir_use = os.path.join(tmp_dir, year, doy)
    cfg_dir_use = os.path.join(cfg_dir, year, doy)
    run_dir_use = os.path.join(run_dir, year, doy)
    
    utils.create_dir(tmp_dir_use)
    utils.create_dir(cfg_dir_use)
    utils.create_dir(run_dir_use)

    ########### DOWNLOAD PRODUCTS
    for data_center in ('ign','whu'):
        operational.download_products_gnss(prod_parent_dir,
                                           srt,end,
                                           AC_names = ("WUM.*FIN",),
                                           prod_types = ("sp3","clk","bia",
                                                         "obx","erp"),
                                           remove_patterns=("ULA",),
                                           archtype ='year/doy',
                                           new_name_conv = True,
                                           parallel_download=4,
                                           archive_center=data_center,
                                           mgex=True,
                                           repro=0,
                                           sorted_mode=False,
                                           return_also_uncompressed_files=True,
                                           ftp_download=False,
                                           dow_manu=False)
        
    ########### GENERATE CONFIG FILE
    def _find_unzip_prod(prod):
        """
        find the right products for a given day, unzip it in the temp dir,
        return the path of the 2 files
        """
    
        prod_ori =  operational.find_IGS_products_files(prod_parent_dir,
                                                    [prod],
                                                    [prod_ac_name], 
                                                    srt,
                                                    severe=False)[0]
        prod_out = files_rw.unzip_gz_Z(prod_ori,out_gzip_dir=tmp_dir)
        
        return prod_out, prod_ori
        
    
    def _change_value_in_cfg(lines_file_inp,key_inp,val_inp):
        """
        change the values in the readed with realine configfile
        """
        for il,l in enumerate(lines_file_inp):
            f = l.split("=")
            if key_inp in f[0]:
                f[1] = val_inp
                lines_file_inp[il] = "= ".join(f) + '\n'    
    
    with open(cfg_template_path) as fil:
        cfg_lines = fil.readlines()
    
    ### find and unzip the right products
    sp3_path,_ = _find_unzip_prod("SP3")
    bia_path,_ = _find_unzip_prod("BIA")
    clk_path,_ = _find_unzip_prod("CLK")
    obx_path,_ = _find_unzip_prod("OBX")
    erp_path,_ = _find_unzip_prod("ERP")
    
    ## alias for basename
    bnm = os.path.basename
    
    ### change the values in the config file template
    _change_value_in_cfg(cfg_lines, "Product directory", tmp_dir)
    _change_value_in_cfg(cfg_lines, "Satellite orbit", bnm(sp3_path))
    _change_value_in_cfg(cfg_lines, "Satellite clock", bnm(clk_path))
    _change_value_in_cfg(cfg_lines, "ERP", bnm(erp_path))
    _change_value_in_cfg(cfg_lines, "Quaternions", bnm(obx_path))
    _change_value_in_cfg(cfg_lines, "Code/phase bias", bnm(bia_path))
    
    date_prod_midfix = year + doy + '0000'

    cfg_name = cfg_prefix + '_' + prod_ac_name + '_' + date_prod_midfix
    cfg_path = os.path.join(cfg_dir,cfg_name)
    
    ### write the good config file
    with open(cfg_path,'w+') as fout:
        for l in cfg_lines:
            fout.write(l)    
    
    ##### RUN THE STUFF
    os.environ['PATH'] += ':'+ bin_dir
    
    cmd = " ".join(("pdp3","--config",cfg_path,rnx_path))
    
    print(cmd)
    
    os.chdir(run_dir_use)
    with Popen([cmd],
               executable="/bin/bash",
               shell=True,
               stdout=PIPE,
               stderr=PIPE) as p:
        while True:
            # Use read1() instead of read() or Popen.communicate() as both blocks until EOF
            # https://docs.python.org/3/library/io.html#io.BufferedIOBase.read1
            text = p.stdout.read1().decode("utf-8")
            print(text, end='', flush=True)
            text = p.stderr.read1().decode("utf-8")
            print(text, end='', flush=True)            
        
    return None


if __name__ == "__main__":
    cfg_template_path = "/home/psakicki/.PRIDE_PPPAR_BIN/config_template"
    rnx_path = '/home/psakicki/aaa_FOURBI/RINEXexemple/BORG00REU_R_20231700000_01D_30S_MO.rnx'
    prod_parent_dir = '/home/psakicki/aaa_FOURBI/PRODS_GNSSexemple/prods/'
    prod_ac_name = "WUM0MGXFIN" 
    tmp_dir = "/home/psakicki/aaa_FOURBI/PRODS_GNSSexemple/unzip_tmp"
    cfg_dir = "/home/psakicki/aaa_FOURBI/pride_pppar_cfg/"
    run_dir = "/home/psakicki/aaa_FOURBI/toto/"
    
    pride_pppar_runner(rnx_path,
                       cfg_template_path,
                       prod_ac_name,
                       prod_parent_dir,
                       tmp_dir,
                       cfg_dir,
                       run_dir,
                       cfg_prefix='pride_pppar_cfg_1a',
                       bin_dir=None)    