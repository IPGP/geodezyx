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

import datetime as dt
import os 
import shutil
import hatanaka
from subprocess import Popen, PIPE
import multiprocessing as mp


import subprocess

#### Import the logger
import logging
log = logging.getLogger(__name__)

def run_command(command):
    # Run the command and capture both stdout and stderr
    process = subprocess.Popen([command],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,
                               executable="/bin/bash",
                               shell=True)

    # Continuously read and print stdout and stderr
    while True:
        # Read a line from stdout
        stdout_line = process.stdout.read1().decode("utf-8") 
        if stdout_line:
            print("STDOUT: %s",stdout_line.strip())
        # Read a line from stderr
        #stderr_line = process.stderr.read1().decode("utf-8") 
        #if stderr_line:
        #    print(f"STDERR: {stderr_line.strip()}", end='', flush=True)

        # Check if the process has finished
        return_code = process.poll()
        if return_code is not None:
            print("RETURN CODE: %s",return_code)
            break

def pride_pppar_runner_mono(rnx_path,
                            cfg_template_path,
                            prod_ac_name,
                            prod_parent_dir,
                            tmp_dir,
                            cfg_dir,
                            run_dir,
                            cfg_prefix='pride_pppar_cfg_1a',
                            mode='K',
                            options_dic={},
                            bin_dir=None,
                            force=False):
    
    if not bin_dir:
        bin_dir = os.path.join(os.environ['HOME'],'.PRIDE_PPPAR_BIN')

    srt = conv.rinexname2dt(rnx_path)

    if srt < dt.datetime(2023,6,4):
        print("DAY SKIPPED")
        return None 

    doy,year = conv.dt2doy_year(srt)
    rnx_file = os.path.basename(rnx_path)
    
    site = rnx_file[:4].upper()

    ########## DEFINE DIRECTORIES
    tmp_dir_use = os.path.join(tmp_dir, year, doy) ### tmp is year/doy only, no site,
                                                   ### because the 'common' dir must be the same 
    cfg_dir_use = os.path.join(cfg_dir, year, doy, site)
    run_dir_use = os.path.join(run_dir,mode,prod_ac_name, site) ### pdp3 add year/doy by itself
    run_dir_ope = os.path.join(run_dir_use, year, doy)
    
    logs_existing = utils.find_recursive(run_dir_ope, "log*" + site.lower())
    
    if len(logs_existing) > 0 and not force:
        print("log exists for %s, skip",rnx_file)
        return None
    else:
        print("no skip",rnx_file)

    utils.create_dir(tmp_dir_use)
    utils.create_dir(cfg_dir_use)
    utils.create_dir(run_dir_use)

    ########### UNCOMPRESS RINEX
    rnx_bnm = os.path.basename(rnx_path) 
    if 'crx' in rnx_bnm or 'd.Z' in rnx_bnm or 'd.gz' in rnx_bnm:  
        rnx_path_tmp = shutil.copy(rnx_path,tmp_dir_use)
        rnx_path_use = str(hatanaka.decompress_on_disk(rnx_path_tmp,
                                                       delete=True))
    else:
        rnx_path_use = rnx_path

    ########### DOWNLOAD PRODUCTS
    if False:
        if 'MGX' in prod_ac_name:
            mgex=True
        else:
            mgex=False

        for data_center in ('ign',):   ##'whu'
            dl_prods_fct = operational.download_products_gnss
            prods=dl_prods_fct(prod_parent_dir,
                               srt,srt,
                               AC_names = (prod_ac_name,),
                               prod_types = ("sp3","clk",
                                             "bia","obx",
                                             "erp"),
                               remove_patterns=("ULA",),
                               archtype ='year/doy',
                               new_name_conv = True,
                               parallel_download=1,
                               archive_center=data_center,
                               mgex=mgex,
                               repro=0,
                               sorted_mode=False,
                               return_also_uncompressed_files=True,
                               ftp_download=False,
                               dow_manu=False)
            print(prods)
            if len(prods) >= 5:
                break
        
    ########### GENERATE CONFIG FILE
    def _find_unzip_prod(prod):
        """
        find the right products for a given day, unzip it in the temp dir,
        return the path of the 2 files
        """
    
        prod_lis =  operational.find_IGS_products_files(prod_parent_dir,
                                                        [prod],
                                                        [prod_ac_name], 
                                                        srt,
                                                        severe=False)
        if len(prod_lis) == 0:
            print("WARN: not prod found")
            prod_out = "Default"
            prod_ori = "Default"
            return prod_out, prod_ori
        elif len(prod_lis) > 1:
            print("WARN: several prod found, keep the 1st one")
            print(prod_lis)
        else:
            pass
        
        prod_ori = prod_lis[0]
        prod_out = files_rw.unzip_gz_Z(prod_ori,out_gzip_dir=tmp_dir_use)
        
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
    _change_value_in_cfg(cfg_lines, "Product directory", tmp_dir_use)
    _change_value_in_cfg(cfg_lines, "Satellite orbit", bnm(sp3_path))
    _change_value_in_cfg(cfg_lines, "Satellite clock", bnm(clk_path))
    _change_value_in_cfg(cfg_lines, "ERP", bnm(erp_path))
    _change_value_in_cfg(cfg_lines, "Quaternions", bnm(obx_path))
    _change_value_in_cfg(cfg_lines, "Code/phase bias", bnm(bia_path))
    
    date_prod_midfix = year + doy + '0000'

    cfg_name = cfg_prefix + '_' + prod_ac_name + '_' + date_prod_midfix
    cfg_path = os.path.join(cfg_dir_use,cfg_name)
    
    ### write the good config file
    with open(cfg_path,'w+') as fout:
        for l in cfg_lines:
            fout.write(l)    
    
    ##### RUN THE STUFF
    os.environ['PATH'] += ':'+ bin_dir
    
    cmd = " ".join(("pdp3","--config",cfg_path,
                    '--mode',mode,
                    rnx_path_use))
    log.info(cmd)

    
    os.chdir(run_dir_use) ## not run_dir_use, pdp3 will goes by itselft to yyyy/doy

    run_command(cmd)

    return None

if __name__ == "__main__":
    cfg_template_path = "/home/ovsgnss/.PRIDE_PPPAR_BIN/config_template"
    rnx_path = '/scratch/calcgnss/temp_stuffs/2402_tests_PF_pride/SNEG/SNEG00REU_R_20232000000_01D_30S_MO.crx.gz' 
    prod_parent_dir = '/scratch/calcgnss/prods_gnss_pride-pppar/prods'  
    #prod_ac_name = "GRG0OPSFIN"
    #prod_ac_name = "COD0MGXFIN"
    prod_ac_name = "WUM0MGXFIN" 
    tmp_dir = "/scratch/calcgnss/pride-pppar_process/tmp" 
    cfg_dir = "/scratch/calcgnss/pride-pppar_process/cfg" 
    run_dir = "/scratch/calcgnss/pride-pppar_process/run" 
    #rnx_list = utils.find_recursive('/scratch/calcgnss/temp_stuffs/2402_tests_PF_pride/','*crx*')
    rnx_list = utils.find_recursive('/vol/ovsg/acqui/GPSOVSG/rinex/2024'  ,'*psa1*d.Z*')
    rnx_list = utils.find_recursive('/home/ovsgnss/050_DATA_GNSS/baiededix/OVPF/2024','*borg*d.Z*')
    multi_process = 14
    
    kwargs_list = []

    dl_prods_fct = operational.download_products_gnss
    prods=dl_prods_fct(prod_parent_dir,
                       dt.datetime(2023,6,1),
                       dt.datetime(2023,7,31),
                       AC_names = (prod_ac_name,),
                       prod_types = ("sp3","clk",
                                     "bia","obx",
                                     "erp"),
                       remove_patterns=("ULA",),
                       archtype ='year/doy',
                       new_name_conv = True,
                       parallel_download=1,
                       archive_center='ign',
                       mgex=True, 
                       repro=0,
                       sorted_mode=False,
                       return_also_uncompressed_files=True,
                       ftp_download=False,
                       dow_manu=False)


    
    for rnx_path in rnx_list:
        kwargs = {'rnx_path' :rnx_path,
                  'cfg_template_path' : cfg_template_path,
                  'prod_ac_name' : prod_ac_name,
                  'prod_parent_dir' : prod_parent_dir,
                  'tmp_dir' : tmp_dir,
                  'cfg_dir' : cfg_dir,
                  'run_dir' : run_dir,
                  'cfg_prefix' : 'pride_pppar_cfg_1a',
                  'bin_dir' : None}
                  
        kwargs_list.append(kwargs)

    def pride_pppar_mp_wrap(kwargs_inp):
        try:
            out_runner = pride_pppar_runner_mono(**kwargs_inp)
            return out_runner 
        except Exception as e:
            log.error("%s raised, RINEX is skiped: %s",type(e).__name__,
                      kwargs_inp['rnx_path'])
            raise e
            

    if multi_process > 1:
        log.info("multiprocessing: %d cores used",multi_process)

    Pool = mp.Pool(processes=multi_process)
    results_raw = [Pool.apply_async(pride_pppar_mp_wrap, args=(x,)) for x in kwargs_list]
    results     = [e.get() for e in results_raw]
