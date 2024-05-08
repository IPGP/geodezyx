#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 18:01:49 2024

@author: psakic
"""

import datetime as dt
#### Import the logger
import logging
import multiprocessing as mp
import os
import shutil
import subprocess

import hatanaka

from geodezyx import conv
from geodezyx import files_rw
#### Import star style
from geodezyx import operational
from geodezyx import utils

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


def dl_brdc_pride_pppar(prod_parent_dir,date_list):
    brdc_lis = []    
    ######### BROADCAST
    date_list_uniq = [conv.round_dt(d,'1D','floor') for d in date_list]
    date_list_uniq = sorted(list(set(date_list_uniq)))
    
    for date in date_list_uniq:
        brdc = operational.download_gnss_rinex({'nav_rt' : ['BRDC']},
                                               prod_parent_dir,
                                               date, date,
                                               archtype = "year/doy",
                                               parallel_download=1,
                                               force=True)
        brdc_lis.append(brdc)
    
    return brdc_lis        
    

        
def dl_prods_pride_pppar(prod_parent_dir,date_list,prod_ac_name):
    
    dl_prods_fct = operational.download_gnss_products
    
    ######### ORBITS CLOCKS ETC...    
    for data_center in ('ign','whu'):   ##'whu'
        if 'MGX' in prod_ac_name:
            mgex=True
        else:
            mgex=False
            
        prods = dl_prods_fct(prod_parent_dir,
                             min(date_list),
                             max(date_list),
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
        
        if len(prods) >= 5:
            log.info("enougth products found: %s",len(prods))
            break
        
        return prods

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
                            force=False,
                            dl_prods=False,
                            default_fallback=False):
    
    if not bin_dir:
        bin_dir = os.path.join(os.environ['HOME'],'.PRIDE_PPPAR_BIN')

    srt = conv.rinexname2dt(rnx_path)

    doy,year = conv.dt2doy_year(srt)
    rnx_file = os.path.basename(rnx_path)
    hourmin_str = str(srt.hour).zfill(2) + str(srt.minute).zfill(2) 
    
    site = rnx_file[:4].upper()
    
    ########## DEFINE DIRECTORIES
    tmp_dir_use = os.path.join(tmp_dir, year, doy) ### tmp is year/doy only, no site,
                                                   ### because the 'common' dir must be the same 
    cfg_dir_use = os.path.join(cfg_dir, year, doy, site)
    run_dir_use = os.path.join(run_dir,mode,prod_ac_name, site) ### pdp3 add year/doy by itself
    run_dir_ope = os.path.join(run_dir_use, year, doy)
    run_dir_fin = run_dir_ope + "_" + hourmin_str
    
    ########### CHECK IF 
    logs_existing = utils.find_recursive(run_dir_fin, "log*" + site.lower())
    
    if len(logs_existing) > 0 and not force:
        log.info("log exists for %s in %s, skip",rnx_file,run_dir_fin)
        return None
    else:
        pass

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
    if dl_prods:
        _ = dl_prods_pride_pppar(prod_parent_dir,[srt],prod_ac_name)
        _ = dl_brdc_pride_pppar(prod_parent_dir,[srt])
        
    ########### MOVE BRDC IN TMP (so pdp3 can handle it)
    ## we must also rename the BKG brdc to fake it as the IGS one
    def _find_unzip_brdc():
        brdc_pattern = "*BRDC00WRD_S_" + year + doy + "*gz"
        prod_dir_year_doy = os.path.join(prod_parent_dir,year,doy)
        brdc_lis = utils.find_recursive(prod_dir_year_doy, brdc_pattern.strip())
        
        if len(brdc_lis) == 1: ## normal case
            brdc_ori = brdc_lis[0]
            brdc_out = files_rw.unzip_gz_Z(brdc_ori,out_gzip_dir=tmp_dir_use)
        elif len(brdc_lis) == 0:
            brdc_ori = None
            brdc_out = None
            log.warning("no brdc. found")
        elif len(brdc_lis) > 1:
            log.warning("several brdc found, keep the 1st one")
            log.warning(brdc_lis)
            brdc_ori = brdc_lis[0]
            brdc_out = files_rw.unzip_gz_Z(brdc_ori,
                                           out_gzip_dir=tmp_dir_use)
        else:
            #### this should never happend
            pass
        
        ##### FINAL STEP: rename the brdc
        if brdc_out:
            brdc_out_new = brdc_out.replace('BRDC00WRD_S_','BRDC00IGS_R_')
            brdc_out = os.rename(brdc_out, brdc_out_new)
        
        
        return brdc_out, brdc_ori
        
        
    ########### GENERATE CONFIG FILE
    def _find_unzip_prod(prod):
        """
        find the right products for a given day, unzip it in the temp dir,
        return the path of the 2 files
        """
        if "ULT" in prod_ac_name: 
            add_hourly_file = True 
        else:
            add_hourly_file = False
        
        find_prod_epoch_ini = srt

        smart_ultra=True
        if "ULT" in prod_ac_name and smart_ultra:
            delta_epoch_max = 23
        else:
            delta_epoch_max = 0

        
        find_prods =operational.find_IGS_products_files 
        
        #### this loop is to find former ULTRA if the latest is missing
        # it works for ULTRA only, for other latencies, delta_epoch_max = 0
        for i_d_epo in range(delta_epoch_max):

            find_prod_epoch = find_prod_epoch_ini - dt.timedelta(seconds=3600*i_d_epo) 

            if smart_ultra and i_d_epo > 0:
                log.warn("no prods found for epoch %s, but we try %i hour before (%s)",
                        find_prod_epoch_ini, i_d_epo, find_prod_epoch)

            prod_lis =  find_prods(prod_parent_dir,
                                   [prod],
                                   [prod_ac_name], 
                                   find_prod_epoch,
                                   severe=False,
                                   regex_old_naming=False,
                                   regex_igs_tfcc_naming=False,
                                   add_hourly_file=add_hourly_file)
            if len(prod_lis) > 0:
                break

        #### differents behaviors depending on the number of files found
        if len(prod_lis) == 1: ## normal case
            prod_ori = prod_lis[0]
            prod_out = files_rw.unzip_gz_Z(prod_ori,out_gzip_dir=tmp_dir_use)
        elif len(prod_lis) == 0 and default_fallback:
            log.warning("no prod. %s found, fallback to Default in cfg file",prod)
            prod_out = "Default"
            prod_ori = "Default"
        elif len(prod_lis) == 0 and not default_fallback:
            log.warning("no prod. %s found, no fallback to Default set, aborting...",prod)
            prod_out = None  
            prod_ori = None 
        elif len(prod_lis) > 1:
            log.warning("several prod found, keep the 1st one")
            log.warning(prod_lis)
            prod_ori = prod_lis[0]
            prod_out = files_rw.unzip_gz_Z(prod_ori,out_gzip_dir=tmp_dir_use)
        else:
            #### this should never happend
            pass
        
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
    brdc_path,_ = _find_unzip_brdc()

    if not any((sp3_path,bia_path,clk_path,obx_path,erp_path)) and not default_fallback:
        log.error("a prod. at least is missing and no fallback to Default is set, abort")
        return  None
    
    ## alias for basename
    bnm = os.path.basename
    
    ### change the values in the config file template
    _change_value_in_cfg(cfg_lines, "Product directory", tmp_dir_use)
    _change_value_in_cfg(cfg_lines, "Satellite orbit", bnm(sp3_path))
    _change_value_in_cfg(cfg_lines, "Satellite clock", bnm(clk_path))
    _change_value_in_cfg(cfg_lines, "ERP", bnm(erp_path))
    _change_value_in_cfg(cfg_lines, "Quaternions", bnm(obx_path))
    _change_value_in_cfg(cfg_lines, "Code/phase bias", bnm(bia_path))
    
    date_prod_midfix = year + doy + hourmin_str  

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
    
    if force and os.path.isdir(run_dir_fin):
        shutil.rmtree(run_dir_fin)

    os.rename(run_dir_ope, run_dir_fin)
    
    return None


def pride_pppar_mp_wrap(kwargs_inp):
    try:
        out_runner = operational.pride_pppar_runner_mono(**kwargs_inp)
        return out_runner 
    except Exception as e:
        log.error("%s raised, RINEX is skiped: %s",
                  type(e).__name__,
                  kwargs_inp['rnx_path'])
        raise e

def pride_pppar_runner(rnx_path_list,
                       cfg_template_path,
                       prod_ac_name,
                       prod_parent_dir,
                       tmp_dir,
                       cfg_dir,
                       run_dir,
                       multi_process = 1,
                       cfg_prefix='pride_pppar_cfg_1a',
                       mode='K',
                       options_dic={},
                       bin_dir=None,
                       force=False,
                       dl_prods=False):
    
    date_list = [conv.rinexname2dt(rnx) - dt.timedelta(seconds=0) for rnx in rnx_path_list] 
    
    _ = dl_prods_pride_pppar(prod_parent_dir,date_list,prod_ac_name)
    _ = dl_brdc_pride_pppar(prod_parent_dir, date_list)
        
    kwargs_list = []
    for rnx_path in rnx_path_list:
        kwargs = {'rnx_path' :rnx_path,
                  'cfg_template_path' : cfg_template_path,
                  'prod_ac_name' : prod_ac_name,
                  'prod_parent_dir' : prod_parent_dir,
                  'tmp_dir' : tmp_dir,
                  'cfg_dir' : cfg_dir,
                  'run_dir' : run_dir,
                  'cfg_prefix' : cfg_prefix,
                  'bin_dir' : bin_dir,
                  'mode' : mode,
                  'options_dic' : options_dic,
                  'force': force,
                  'dl_prods' : dl_prods}
                  
        kwargs_list.append(kwargs)
            
    if multi_process > 1:
        log.info("multiprocessing: %d cores used",multi_process)
    
    Pool = mp.Pool(processes=multi_process)
    results_raw = [Pool.apply_async(pride_pppar_mp_wrap, args=(x,)) for x in kwargs_list]
    results     = [e.get() for e in results_raw]

