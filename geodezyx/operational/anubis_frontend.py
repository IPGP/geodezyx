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
                  silent=False,
                  download_nav=True,
                  download_sp3=True,
                  anubis_path="/home/psakicki/SOFTWARE/ANUBIS/anubis-3.3.3613-lin-shared-64b"):
    """
    Run an Anubis quality check
    Designed for Anubis 3.3

    Parameters
    ----------
    rnx_inp : str or list
        Input RINEXs. can be a path list, or the path of the parent archive directory
    out_dir_main : str
        output main directory.
    xml_config_generic : str
        path of the generic XML configuration file.
        will be stored in <out_dir_main>/inp
    silent : bool, optional
        if True, do not run the Anubis QC.
        Anubis QC rsults are stored in stored in <out_dir_main>/out
        The default is False.
    download_nav : bool, optional
        Download automatically the Broadcast navigation files. 
        will be stored in <out_dir_main>/nav
        The default is True.
    download_sp3 : TYPE, optional
        Download automatically the SP3 orbit files. MGEX CODE are used per default.
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

        nav_path = rnx_path[:-1] +  "n"
    
        date_start = conv.rinexname2dt(rnx_path)
        date_end = date_start + dt.timedelta(days=1) - dt.timedelta(seconds=1)
        
        date_start_str = conv.dt2str(date_start,"%Y-%m-%d %H:%M:%S")
        date_doy_str   = conv.dt2str(date_start,"%Y%j")
        date_end_str   = conv.dt2str(date_end,"%Y-%m-%d %H:%M:%S")
        date_name_str  = conv.dt2str(date_start,"%Y%m%d")
        
        ### site name
        if re.search(conv.rinex_regex_new_name(),rnx_path):
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
        
        ### manage the brdc-file download
        nav_path = ""
        if download_nav:
            statdico = dict()
            statdico['brdc'] = ['BRDC']
            brdc_list = operational.multi_downloader_rinex(statdico,
                                                           nav_dir,
                                                           date_start,date_end,
                                                           "/",
                                                           parallel_download=1,
                                                           sorted_mode=0)    
            if brdc_list:
                nav_path = brdc_list[0]
            else:
                nav_path = ""
        
        
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
        if not silent:
            #os.chdir(os.path.dirname(xml_path_ope))
            command = anubis_path + " -x " + xml_path_ope
            
            log.info("command: %s",command)
            
            process1 = subprocess.run(command,
                                      capture_output=True,
                                      text=True,
                                      shell=True)
            
    return XML_INP_LIST