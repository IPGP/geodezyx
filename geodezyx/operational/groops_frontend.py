#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 10:54:49 2023

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


#### Import the logger
import logging
log = logging.getLogger(__name__)




def log_subprocess_output(pipe,logger=None,file=None,file2=None):
    """
    Intern fuction to write the stdout/err in the console logger/logfile
    """
    for line in iter(pipe.readline, b''): # b'\n'-separated lines
        line_clean = line.decode().strip()
        if file:
            file.write(line_clean + "\n") 
        if file2:
            file2.write(line_clean + "\n") 
        if logger:
            logger(line_clean)
    return 
    
    
    
def subprocess_frontend2(cmd_in,
                         save_log=True,
                         log_dir=None,
                         log_name_out="out.log",
                         log_name_err="err.log",
                         logname_timestamp=True,
                         err_also_in_outfile=True):
    """
    A generic frontend to run an extern command through subprocess
    and write the stdout and stderr outputs in log files.
    
    Parameters
    ----------
    cmd_in : str
        The subprocess command.
    save_log : str, optional
        export as log the stdout/stderr in files. The default is True.
    log_dir : str, optional
        directory where the logs will be stored. The default is None.
    log_name_out : str, optional
        filename of the stdout log. The default is "out.log".
    log_name_err : str, optional
        filename of the stderr log. The default is "err.log".
    logname_timestamp : str, optional
        add a timestamp as prefix. The default is True.
    err_also_in_outfile : bool, optional
        write also the stderr in the stdout log. The default is True.

    Returns
    -------
    exitcode : int
        the exit code of the command.
        
    Notes
    -----
    Inspired by:
    https://stackoverflow.com/questions/21953835/run-subprocess-and-print-output-to-logging
    And for the threading:
    https://stackoverflow.com/questions/6809590/merging-a-python-scripts-subprocess-stdout-and-stderr-while-keeping-them-disti

    """
    
    from subprocess import Popen, PIPE, STDOUT
    now = utils.get_timestamp()
    
    #### manage the paths of the output logs
    if save_log:
        if not log_dir:
            log_dir = os.getcwd()
    
        if logname_timestamp:
            prefix = now + "_"
        else:
            prefix = ""
    
        err_file = open(log_dir + "/" + prefix + log_name_err, "w+")
        out_file = open(log_dir + "/" + prefix + log_name_out, "w+")
        
        if err_also_in_outfile:
            err_in_outfile = out_file
        else:
            err_in_outfile = None
            
    else:
        out_file = None
        err_file = None
        err_in_outfile = None
                                
    ## arguments must be splitted 
    if type(cmd_in) is str:
        cmd_split = cmd_in.split()
    else:
        cmd_split = cmd_in
    
    
    ####### Run the command here  #####################
    process = Popen(cmd_split,
                    stdout=PIPE,
                    stderr=PIPE)
    ##################################################
    
    
    ### Get out/err simultaneously
    stdout_thread = Thread(target=log_subprocess_output,
                           args=(process.stdout,log.info,out_file))
    stderr_thread = Thread(target=log_subprocess_output,
                           args=(process.stderr,log.error,err_file,
                                 err_in_outfile))
    
    stderr_thread.start()
    stdout_thread.start()

    #while (process.poll() is None) and stdout_thread.is_alive() and stderr_thread.is_alive():
    while stdout_thread.is_alive() or stderr_thread.is_alive():
        pass ### do nothing while the threads are runing

    exitcode = process.wait() # 0 means success

    if save_log:
        out_file.close()
        err_file.close()       
    
    return exitcode


def subprocess_frontend3(cmd_in,
                         save_log=True,
                         log_dir=None,
                         log_name_out="out.log",
                         logname_timestamp=True):
    """
    A generic frontend to run an extern command through subprocess
    and write the stdout and stderr outputs in log files.
    
    Parameters
    ----------
    cmd_in : str
        DESCRIPTION.
    save_log : str, optional
        export as log the stdout/stderr in files. The default is True.
    log_dir : str, optional
        directory where the logs will be stored. The default is None.
    log_name_out : str, optional
        filename of the stdout log. The default is "out.log".
    logname_timestamp : str, optional
        add a timestamp as prefix. The default is True.
        
    Returns
    -------
    exitcode : int
        the exit code of the command.
        
    Notes
    -----
    Inspired by:
    https://stackoverflow.com/questions/21953835/run-subprocess-and-print-output-to-logging
    And for the threading:
    https://stackoverflow.com/questions/6809590/merging-a-python-scripts-subprocess-stdout-and-stderr-while-keeping-them-disti

    """
    
    from subprocess import Popen, PIPE, STDOUT
    now = utils.get_timestamp()
    
    #### manage the paths of the output logs
    if save_log:
        if not log_dir:
            log_dir = os.getcwd()
    
        if logname_timestamp:
            prefix = now + "_"
        else:
            prefix = ""
    
        out_file = open(log_dir + "/" + prefix + log_name_out, "w+")
        
    else:
        out_file = None
                                
    ## arguments must be splitted 
    if type(cmd_in) is str:
        cmd_split = cmd_in.split()
    else:
        cmd_split = cmd_in
    
    
    ####### Run the command here  #####################
    process = Popen(cmd_split,
                    stdout=PIPE,
                    stderr=STDOUT)
    # Get out/err simultaneously with stderr=STDOUT
    ##################################################
    
    ### Get out/err simultaneously
    log_subprocess_output(process.stdout,log.info,out_file)

    exitcode = process.wait() # 0 means success
    
    if save_log:
        out_file.close()
    
    return exitcode

def groops_basic_runner(xml_config_path="",
                        global_var_dict=dict(),
                        xml_var_dict=dict(),
                        quiet=False,
                        log_dir="/home/ovsgnss/020_CALC/groops_process/031_groops_frontend_logs",
                        groops_bin_path='/opt/softs_gnss/groops/bin/groops'):


    global_args_str = ""
    for key, val in global_var_dict.items():
        global_args_str = global_args_str + " ".join((" --global",key+"="+str(val)))
        
        
    command = " ".join((groops_bin_path,global_args_str,xml_config_path))
    log.info("groops command:")
    log.info("%s",command)
    
    if not quiet:
        subprocess_frontend2(command,
                             save_log=True,
                             log_dir=log_dir,
                             log_name_out="out.log",
                             logname_timestamp=True)

def vmf_tropo_downloader(output_dir,
                         startdate = dt.datetime(2019,1,1),
                         enddate = dt.datetime(2019,1,1),
                         model = "VMF3",
                         version = "OP"):
    
    Dates_range = conv.dt_range(startdate,enddate)
    
    url_root = "https://vmf.geo.tuwien.ac.at"
    dir_root = output_dir
    subdir   = "trop_products/GRID/1x1/"
    
    Files_hour_list = []
    for currdate in Dates_range:
        subdir_date = os.path.join(subdir,
                                   model,
                                   model + "_" +  version,
                                   str(currdate.year))
        
        
        url_date = os.path.join(url_root,subdir_date)
        dir_date = os.path.join(dir_root,subdir_date)
        utils.create_dir(dir_date)

        date_str = conv.dt2str(currdate,"%Y%m%d") 
        
        for hourext in ("H00","H06","H12","H18"):      
            filename_hour = model + "_" + date_str + "." + hourext
            url_hour = os.path.join(url_date,filename_hour)
            outfile_hour = os.path.join(dir_date,filename_hour)
            
            log.info("download %s to %s",url_hour,outfile_hour)
            
            Files_hour_list.append(outfile_hour)
            if not utils.empty_file_check(outfile_hour):
                log.info("%s exists, download skipped",outfile_hour)
                continue
            else:
                DownObjt = urllib.request.urlretrieve(url_hour, 
                                                      outfile_hour)   
    return Files_hour_list
    
#%%

###############################################################################
###############################################################################

def groops_ppp_full_runner(rinex_path,
                           project_name,
                           igs_ac_10char,
                           vmf_tropo_root_dir = "/scratch/calcgnss/prods_tropo_vmf3",
                           prods_gnss_root_dir = "/scratch/calcgnss/prods_gnss/",
                           config_files_root_dir = "/opt/softs_gnss/groops/config/030_prototype3/"):
    
    ###############################################################################
    ######## Set python fct variables
    
    debug_sleep_time = 2
    quiet=False


    prods_gnss_dir = os.path.join(prods_gnss_root_dir,igs_ac_10char)
    
    
    if igs_ac_10char[4:7] == "MGX":
        mgex=True
    else:
        mgex=False
    
    if (igs_ac_10char[4] == "R") and (igs_ac_10char[6] == "3"):
        repro=3
    else:
        repro=0
    
    project_name_use = project_name + "_" + igs_ac_10char 
    
    site4char = os.path.basename(rinex_path)[:4]
    date_rnx = conv.rinexname2dt(rinex_path)
    date_rnx_mjd = int(conv.dt2MJD(date_rnx))
    date_rin_ymd = conv.dt2str(date_rnx,"%Y-%m-%d")
    
    ###############################################################################
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ######## VMF TROPO
    
    ###############################################################################
    ######## Download Tropo
    log.info("\n****** VMF troposphere grids download *****************************")
    
    VMF_files  = vmf_tropo_downloader(vmf_tropo_root_dir,
                                      startdate = date_rnx  - dt.timedelta(days=0),
                                      enddate   = date_rnx)
    
    ###############################################################################
    ######## Convert Tropo
    log.info("\n****** VMF troposphere grids conversion ***************************")
    
    xml_config_path=config_files_root_dir + "010_groopsConvert_tropoVmf3_v03a.xml"
    
    global_var_dict = dict()
    global_var_dict["timeStart"] = date_rnx_mjd
    global_var_dict["timeEnd"]   = date_rnx_mjd
    global_var_dict["groopsInpVmfGridFile1"] = VMF_files[-4]
    global_var_dict["groopsInpVmfGridFile2"] = VMF_files[-3]
    global_var_dict["groopsInpVmfGridFile3"] = VMF_files[-2]
    global_var_dict["groopsInpVmfGridFile4"] = VMF_files[-1]
    
    groops_basic_runner(xml_config_path,
                        global_var_dict,
                        quiet=quiet)
    
    ###############################################################################
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ######## PRODUCTS
    
    ###############################################################################
    ######## Check if converted products exists
    log.info("****** Converted GNSS products existance check ********************")
    
    paaath = "/home/ovsgnss/020_CALC/groops_process/021_conv_igs_prods"
    conv_prod_dir = os.path.join(paaath,igs_ac_10char,date_rin_ymd)
    
    All_conv_prods = utils.find_recursive(conv_prod_dir,"*dat")
    Clk_conv_prods = utils.find_recursive(conv_prod_dir,"clock*dat")
    Orb_conv_prods = utils.find_recursive(conv_prod_dir,"orbit*dat")
        
    if len(Orb_conv_prods) > 5:
        download_prods = False
        log.info("Products download/conversion skipped, %s converted files found in %s",
                 len(All_conv_prods),conv_prod_dir)
    else:
        download_prods = True
        log.info("Products download/conversion will start, no files found in %s",conv_prod_dir)
    
    ###############################################################################
    ######## Download Products
    
    if download_prods:
        log.info("****** GNSS products download *********************************")
    
        download_fct = operational.multi_downloader_orbs_clks_2
        
        Prod_types  = ["sp3","clk","bia","obx"]
        AC_names    = [igs_ac_10char]
        Prods_out_list = download_fct(prods_gnss_dir,
                                      date_rnx,
                                      date_rnx,
                                      AC_names,
                                      Prod_types,
                                      archtype ='week',
                                      archive_center="cddis",
                                      parallel_download=1,
                                      repro=repro,
                                      mgex=mgex)
        
        if not Prods_out_list:
            log.warning("No downloaded products found locally, abort")  
            raise Exception
        
        ###### Check the downloaded products (IGS AC ID)
        Igs_ac_10char_list_out = [os.path.basename(e)[:10] for e in Prods_out_list]
        Igs_ac_10char_list_out = list(set(Igs_ac_10char_list_out))
            
        if len(Igs_ac_10char_list_out) > 1:
            log.warning("Several downloaded product types found: %s",
                        Igs_ac_10char_list_out)
               
        if  Prods_out_list and (igs_ac_10char != Igs_ac_10char_list_out[0]):
            log.warning("Product type downloaded != requested: %s",
                        Igs_ac_10char_list_out[0],
                        igs_ac_10char)
    
    ###############################################################################
    ######## Convert Products: clocks
    if download_prods:
        log.info("****** GNSS products conversion: clocks *******************")
    
        xml_config_path=config_files_root_dir + "021_groopsConvert_IgsProds_clock_v03a.xml"
        
        global_var_dict = dict()
        global_var_dict["timeStart"] = date_rnx_mjd
        global_var_dict["timeEnd"] = date_rnx_mjd + 1
        global_var_dict["igsAC10Char"] = igs_ac_10char
        global_var_dict["inpIgsProdsDir"] = prods_gnss_dir
        groops_basic_runner(xml_config_path,
                            global_var_dict,
                            quiet=quiet)

    
        time.sleep(debug_sleep_time)    
    
    ###############################################################################
    ######## Convert Products: bias
    if download_prods:
        log.info("****** GNSS products conversion: bias *******************")
    
        xml_config_path=config_files_root_dir + "022_groopsConvert_IgsProds_bias_v03a.xml"
        
        groops_basic_runner(xml_config_path,
                            global_var_dict,
                            quiet=quiet)

        time.sleep(debug_sleep_time)    
    
    ###############################################################################
    ######## Convert Products: orbits
    if download_prods:
        log.info("****** GNSS products conversion: orbits *******************")
    
        xml_config_path=config_files_root_dir + "023_groopsConvert_IgsProds_orbit_v03a.xml"
        
        groops_basic_runner(xml_config_path,
                            global_var_dict,
                            quiet=quiet)

        time.sleep(debug_sleep_time)    

    ###############################################################################
    ######## Convert Products: alternative attitude
    if download_prods:
        log.info("****** GNSS products conversion: alternative attitude *****")
    
        xml_config_path=config_files_root_dir + "024_groopsConvert_IgsProds_alt_attitude_v03a.xml"
        
        groops_basic_runner(xml_config_path,
                            global_var_dict,
                            quiet=quiet)

        time.sleep(debug_sleep_time)    



    ###############################################################################
    ######## Convert Products: attitude
    if download_prods:
        log.info("****** GNSS products conversion: attitude *****")
    
        xml_config_path=config_files_root_dir + "025_groopsConvert_IgsProds_attitude_v03a.xml"
        
        groops_basic_runner(xml_config_path,
                            global_var_dict,
                            quiet=quiet)

        time.sleep(debug_sleep_time)    


    
    ###############################################################################
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ######## RINEX & STATION INFO
    
    ###############################################################################
    ######## convert sitelog > station info
    
    # log.info("\n****** sitelog > station info conversion **************************")
    
    # xmltree = et.parse(xml_config_generic)
    # # get input files
    # xmltree.find('.//rinexo').text
    
    # xml_config_path=config_files_root_dir + "022_groopsConvert_GnssObs_v02b.xml"
    
    # global_var_dict = dict()
    # global_var_dict["timeStart"] = date_rnx_mjd
    # global_var_dict["timeEnd"]   = date_rnx_mjd + 1
    # global_var_dict["site4Char"] = site4char
    # global_var_dict["inpGnssObsRnxFile"] = rinex_path
    
    # groops_basic_runner(xml_config_path,
    #                     global_var_dict)
    
    ###############################################################################
    ######## Convert RINEX
    
    log.info("****** RINEX observation conversion *******************************")
    
    xml_config_path=config_files_root_dir + "031_groopsConvert_GnssObs_v03a.xml"
    
    global_var_dict = dict()
    global_var_dict["timeStart"] = date_rnx_mjd
    global_var_dict["timeEnd"]   = date_rnx_mjd + 1
    global_var_dict["site4Char"] = site4char
    global_var_dict["inpGnssObsRnxFile"] = rinex_path
    
    groops_basic_runner(xml_config_path,
                        global_var_dict,
                        quiet=quiet)

    
    ###############################################################################
    ######## Edit station list
    
    log.info("****** Station list edition ***************************************")
    
    station_list_path='/opt/softs_gnss/groops/stationlists/station_list_OPERA_01a.txt'
    F=open(station_list_path,"w+")
    F.write(site4char)
    F.close()
    
    
    ###############################################################################
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ######## Run Processing
    
    log.info("****** Processing run *********************************************")
    
    xml_config_path=config_files_root_dir + "040_groopsGnssProcessing_v03a.xml"
    
    global_var_dict = dict()
    global_var_dict["igsAC10Char"] = igs_ac_10char
    global_var_dict["projectName"] = project_name_use
    global_var_dict["timeStart"] = date_rnx_mjd
    global_var_dict["timeEnd"]   = date_rnx_mjd + 1
    
    groops_basic_runner(xml_config_path,
                        global_var_dict,
                        quiet=quiet)

    
    
    return

### ===========================================================================
### ===========================================================================

rinex_path = "/vol/ovpf/miroir_ovpf/DonneesAcquisition/geodesie/GPSData/2023/008/fjag0080.23d.Z"
rinex_root_path = "/vol/ovpf/miroir_ovpf/DonneesAcquisition/geodesie/GPSData/"

igs_ac_10char = "IGS0OPSFIN"
igs_ac_10char = "IGS2R03FIN"
igs_ac_10char = "COD0R03FIN"
igs_ac_10char = "GRG6RE3FIN"
igs_ac_10char = "TUG0R03FIN"
igs_ac_10char = "GFZ0MGXRAP"
igs_ac_10char = "COD0MGXFIN"
igs_ac_10char = "GRG0MGXFIN"
igs_ac_10char = "GRG0OPSFIN"

Rinexs = operational.rinex_finder(rinex_root_path,
                                  specific_sites=["fjag"], #,"borg","DERG","FJAG"],
                                  start_epoch=dt.datetime(2023,1,1),
                                  end_epoch=dt.datetime(2023,1,15))

print(Rinexs)
project_name="test_calc_PF_03"

for rnx in Rinexs:
    groops_ppp_full_runner(rnx,
                           project_name,
                           igs_ac_10char,
                           vmf_tropo_root_dir = "/scratch/calcgnss/prods_tropo_vmf3",
                           prods_gnss_root_dir = "/scratch/calcgnss/prods_gnss/",
                           config_files_root_dir = "/opt/softs_gnss/groops/config/030_prototype3/")