#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 15:48:26 2022

@author: psakicki
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
import os 
import re
import subprocess
import time

#manage XML
from xml.etree import ElementTree as et
from xml.etree.ElementTree import Element

#### geodeZYX modules
from geodezyx import conv
from geodezyx import operational
from geodezyx import utils

#### Import the logger
import logging
log = logging.getLogger(__name__)

##########  END IMPORT  ##########

def anubis_runner(rnx_inp,
                  out_dir_main,
                  xml_config_generic,
                  period=None,
                  interval=None,
                  quiet=False,
                  download_nav=True,
                  download_sp3=True,
                  anubis_path="/home/psakicki/SOFTWARE/ANUBIS/anubis-3.3.3613-lin-shared-64b"):
    """
    Run an Anubis quality check
    Designed for Anubis 3.3

    Parameters
    ----------
    rnx_inp : str or list
        Input RINEXs. can be a RINEX path list,
        or the path of the parent archive directory
    out_dir_main : str
        output main directory.
    xml_config_generic : str
        path of the generic XML configuration file.
        will be stored in <out_dir_main>/inp
    period : None or int
        nominal file period in the RINEX (in hours)
        if None is given, guess based on the RINEX name
        The default is None.    
    interval : None or int
        nominal data interval in the RINEX (in sec)
        if None is given, guess based on the RINEX name
        The default is None.
    quiet : bool, optional
        if True, do not run the Anubis QC.
        Anubis QC rsults are stored in stored in <out_dir_main>/out
        The default is False.
    download_nav : bool, optional
        Download automatically the Broadcast navigation files. 
        will be stored in <out_dir_main>/nav
        The default is True.
    download_sp3 : bool, optional
        Download automatically the SP3 orbit files. CODE's MGEX are used per default.
        will be stored in <out_dir_main>/nav
        The default is True.
    anubis_path : str, optional
        path of the Anubis executable. 
        The default is "/home/psakicki/SOFTWARE/ANUBIS/anubis-3.3.3613-lin-shared-64b".

    Returns
    -------
    XML_INP_LIST : list
        list of the generated XML config files.

    """
    
    
    if utils.is_iterable(rnx_inp):
        RNXLIST = rnx_inp 
    else:
        RNXLIST = operational.rinex_finder(rnx_inp) 
    
    XML_INP_LIST = []
        
    for rnx_path in RNXLIST:
        log.info("current RINEX: %s",rnx_path)
        
        rnx_name = os.path.basename(rnx_path)
        
        ###### MANAGE PERIOD AND INTERVAL
        ### check if a valid period/interval is avaiable
        if ((period is None) or (interval is None)) and (not conv.rinex_regex_search_tester(rnx_name,False,True)):
            log.error("the input RINEX %s does not have a long name, unable to detect period/interval, set them manually",rnx_name)
            raise Exception

        if period is None:
            period_ok = conv.period_from_rinex_name(rnx_name)
            print("PERIOD",period_ok)
        else:
            period_ok = period
            
        if interval is None:
            interval_ok = conv.interval_from_rinex_name(rnx_name)
            print("INTERVAL",interval_ok)
        else:
            interval_ok = interval
            
        log.info("interval:%s,period:%s",interval_ok,period_ok)
        
        ### manage dates
        date_start = conv.rinexname2dt(rnx_path)
        date_end   = date_start + dt.timedelta(seconds=period_ok) - dt.timedelta(seconds=1)
        
        date_start_str = conv.dt2str(date_start,"%Y-%m-%d %H:%M:%S")
        date_doy_str   = conv.dt2str(date_start,"%Y%j%H%M")
        date_doy_short_str = conv.dt2str(date_start,"%Y%j")
        date_end_str   = conv.dt2str(date_end,"%Y-%m-%d %H:%M:%S")
        date_name_str  = conv.dt2str(date_start,"%Y%m%d%H%M")
        
        ### site name
        if re.search(conv.rinex_regex_long_name(),rnx_path):
            site = os.path.basename(rnx_path)[:9].upper()
        else:
            site = os.path.basename(rnx_path)[:4].upper()
    
        site_date = site + "_" + date_name_str + "_" + date_doy_str
    
        ### create directories
        utils.create_dir(out_dir_main)
        nav_dir = utils.create_dir(out_dir_main + "/nav")
        inp_dir = utils.create_dir(out_dir_main + "/inp")
        out_dir = utils.create_dir(out_dir_main + "/out")
        
        xml_path_ope = inp_dir + "/" + site_date + "_inp.xml"
    
        out_xml = os.path.join(out_dir, site_date + ".xml") 
        out_xtr = os.path.join(out_dir, site_date + ".xtr") 
        out_log = os.path.join(out_dir, site_date + ".log") 

        ######## MANAGE NAV FILE #########                
        ### manage the brdc-file download
        nav_path = ""
        
        potential_nav_file = "BRDC00GOP_R_" + date_doy_short_str + "0000_01D_MN.rnx.gz"
        potential_nav_path = os.path.join(nav_dir,potential_nav_file)
        
        nav_file_exists = os.path.isfile(potential_nav_path)
        
        #### the nav file already exists
        print(potential_nav_path)
        if nav_file_exists:
            log.info("%s already exists, download skipped ;)",potential_nav_file)
            nav_path = potential_nav_path
        
        #### the nav file does not exsits but we want to download it
        elif not nav_file_exists and download_nav:
            statdico = dict()
            statdico['brdc'] = ['BRDC']
            brdc_list = operational.multi_downloader_rinex(statdico,
                                                           nav_dir,
                                                           date_start,date_end,
                                                           "/",
                                                           parallel_download=1,
                                                           sorted_mode=0,
                                                           final_archive_for_sup_check=nav_dir)    
            if brdc_list:
                nav_path = brdc_list[0]
            else:
                nav_path = ""
        
        #### we do not consider the nav file
        else:
            nav_path = ""
            
        
        ######## MANAGE SP3 FILE #########                
        ### manage the SP3-file download
        sp3_path = ""
        if download_sp3:
            sp3_list = operational.multi_downloader_orbs_clks_2(nav_dir,
                                                                date_start,date_end,
                                                                AC_names = ["COD"],
                                                                prod_types = ["SP3"],
                                                                remove_patterns=("ULA",),
                                                                archtype ='/',
                                                                new_name_conv = True,
                                                                parallel_download=1,
                                                                archive_center='ign',
                                                                mgex=True,
                                                                repro=0,
                                                                sorted_mode=False,
                                                                return_also_uncompressed_files=True,
                                                                ftp_download=False,
                                                                dow_manu=False) 
            
            if sp3_list:
                sp3_path = sp3_list[0]    
            else:
                sp3_path = ""
            
        ### start to modifiy the XML
        xmltree = et.parse(xml_config_generic)
    
        if type(nav_path) is list: #hotfix 230406 ... must be investigated!
            nav_path = nav_path[0]

        # set input files
        xmltree.find('.//rinexo').text = rnx_path
        xmltree.find('.//rinexn').text = nav_path
        xmltree.find('.//sp3').text    = sp3_path
    
        # set output files
        xmltree.find('.//xtr').text = out_xtr
        xmltree.find('.//xml').text = out_xml
        xmltree.find('.//log').text = out_log
    
        # set dates
        xmltree.find('.//beg').text = date_start_str
        xmltree.find('.//end').text = date_end_str
        xmltree.find('.//int').text = str(interval_ok)
        xmltree.find('.//rec').text = site
    
        # set site infos
        elem_site = xmltree.find('site')
        elem_site.set('id', site)
        elem_site.set('name', site)
        elem_site.set('domes', "00000X000")
    
        # write output operational xml 
        xmltree.write(xml_path_ope)
        XML_INP_LIST.append(xml_path_ope)
        
        ### run Anubis
        if not quiet:
            #os.chdir(os.path.dirname(xml_path_ope))
            command = anubis_path + " -x " + xml_path_ope
            
            log.info("command: %s",command)
            
            process1 = subprocess.run(command,
                                      capture_output=True,
                                      text=True,
                                      shell=True)
            
            time.sleep(2)
            
    return XML_INP_LIST
