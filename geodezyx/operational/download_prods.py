#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic
This sub-module of geodezyx.operational contains functions to download
gnss data and products from distant IGS servers. 
it can be imported directly with:
from geodezyx import operational
The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License
Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
from ftplib import FTP, FTP_TLS
# import glob
import itertools
import multiprocessing as mp
import numpy as np
import os
# import pandas as pd 
import re
import shutil
# import urllib
# import ftplib

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

import geodezyx.operational.download_utils as dlutils

#### Import star style
# from geodezyx import *                   # Import the GeodeZYX modules
# from geodezyx.externlib import *         # Import the external modules
# from geodezyx.megalib.megalib import *   # Import the legacy modules names


#### Import the logger
import logging
log = logging.getLogger(__name__)

##########  END IMPORT  ##########



#  _____               _            _         _____                      _                 _           
# |  __ \             | |          | |       |  __ \                    | |               | |          
# | |__) | __ ___   __| |_   _  ___| |_ ___  | |  | | _____      ___ __ | | ___   __ _  __| | ___ _ __ 
# |  ___/ '__/ _ \ / _` | | | |/ __| __/ __| | |  | |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |/ _ \ '__|
# | |   | | | (_) | (_| | |_| | (__| |_\__ \ | |__| | (_) \ V  V /| | | | | (_) | (_| | (_| |  __/ |   
# |_|   |_|  \___/ \__,_|\__,_|\___|\__|___/ |_____/ \___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|\___|_|   
                                                                                                      
                                                                                                      

############################################################################
######## PRODUCTS DOWNLOADER
############################################################################

def download_products_gnss(archive_dir,
                           startdate,enddate,
                           AC_names = ("wum","cod"),
                           prod_types = ("sp3","clk"),
                           remove_patterns=("ULA",),
                           archtype ='week',
                           new_name_conv = True,
                           parallel_download=4,
                           archive_center='ign',
                           mgex=False,
                           repro=0,
                           sorted_mode=False,
                           return_also_uncompressed_files=True,
                           ftp_download=False,
                           dow_manu=False):
    """
    dow_manu = False, no dow manu, consider the converted dow from the time span, regular case
    dow_manu = None, no dow in the REGEX, the crawler will search only for the week
    dow_manu = 0 or 7: the dow in question    
    
    
    to control the lattency with the new naming convention, simply add it completly in the AC name
    e.g. IGS0OPSRAP
    
    Note
    ----
    The new naming convention has been fully adopted since GPS Week 2238-0

    
    """
    
    if mgex:
        mgex_str = "mgex/"
    else:
        mgex_str = ""
        
    if repro:
        repro_str = "repro" + str(repro) + "/"
    else:
        repro_str = ""
        
    if not utils.is_iterable(remove_patterns):
        remove_patterns = [remove_patterns]
    
    secure_ftp = False
    
    log.info("data center used : %s",archive_center)
    log.info("mgex/repro : %s/%s",mgex,repro)
    
    if archive_center == "cddis":
        arch_center_main    = 'gdc.cddis.eosdis.nasa.gov'
        arch_center_basedir = '/pub/gps/products/' + mgex_str
        ftp_download = True
        secure_ftp = True
        parallel_download = 1
        log.info('cddis as data center, FTP and no parallel download forced')
        
    elif archive_center == "cddis_glonass":
        arch_center_main    = 'cddis.gsfc.nasa.gov'
        arch_center_basedir = '/pub/glonass/products/' + mgex_str
        
    elif archive_center == "ign":
        arch_center_main    = 'igs.ign.fr'
        arch_center_basedir = '/pub/igs/products/' + mgex_str  
        
    elif archive_center == "ign_iono":
        arch_center_main    = 'igs-rf.ign.fr'
        arch_center_basedir = '/pub/'  

    elif archive_center == "ensg":
        arch_center_main    = 'igs.ensg.ign.fr'
        arch_center_basedir = '/pub/igs/products/' + mgex_str
        
    elif archive_center == "whu":
        arch_center_main    = "igs.gnsswhu.cn"
        arch_center_basedir = "/pub/gps/products/" + mgex_str
        
    elif archive_center == "ign_rf":
        arch_center_main    = 'igs-rf.ign.fr'
        arch_center_basedir = '/pub/' + mgex_str  

    elif archive_center == "ensg_rf":
        arch_center_main    = 'igs-rf.ensg.ign.fr'
        arch_center_basedir = '/pub/' + mgex_str
        
    elif archive_center == "acc_xpr_mgex_cmb":
        ##### DO NOT WORK !!!!!
        ## and this archive is not mainteed anyway (2023-01)
        arch_center_main    = 'http://igsacc.s3-eu-central-1.amazonaws.com/products/mgex/final/2069/igm20694.sp3.Z'
        ftp_download = False
        log.info('ACC experimental mgex combi. as data center, HTTP download forced')


    Dates_list = conv.dt_range(startdate,enddate)

    wwww_dir_previous = None
    pool = mp.Pool(processes=parallel_download) 

    ## internal fct to create the FTP objects
    
    class MyFTP_TLS(FTP_TLS):
        """Explicit FTPS, with shared TLS session"""
        ### This new class is to avoid the error 
        ### ssl.SSLEOFError: EOF occurred in violation of protocol (_ssl.c:2396)
        ### source:
        ### https://stackoverflow.com/questions/14659154/ftps-with-python-ftplib-session-reuse-required
        def ntransfercmd(self, cmd, rest=None):
            conn, size = FTP.ntransfercmd(self, cmd, rest)
            if self._prot_p:
                conn = self.context.wrap_socket(conn,
                                                server_hostname=self.host,
                                                session=self.sock.session)  # this is the fix
            return conn, size
        
    def ftp_objt_create(secure_ftp_inp,chdir=""):
        
        # define the right constructor
        if secure_ftp_inp:
            ftp_constuctor = MyFTP_TLS
            #ftp=ftp_constuctor()
            #ftp.set_debuglevel(2)
            #ftp.connect(arch_center_main)
            #ftp.login('anonymous','')
            #ftp.prot_p()
        else:     
            ftp_constuctor = FTP
            #ftp = ftp_constuctor(arch_center_main)
            #ftp.login()
            
        ## create a list of FTP object for multiple downloads
        Ftp_obj_list_out = [ftp_constuctor(arch_center_main) for i in range(parallel_download)]
        if secure_ftp:
            [f.login('anonymous','') for f in Ftp_obj_list_out]    
            [f.prot_p() for f in Ftp_obj_list_out]    
        else:
            [f.login() for f in Ftp_obj_list_out]    
            
        # define the main obj for crawling
        ftp_main = Ftp_obj_list_out[0]
        
        # change the directory of the main ftp obj if we ask for it
        if chdir:
            log.info("Move to: %s",chdir)
            ftp_main.cwd(chdir)
        
        return ftp_main, Ftp_obj_list_out
    
    ###################################################################
    ########### Remote file search      

    Potential_localfiles_list_all = []
    
    ### check if the pattern of the wished products are in the listed daily files
    for ipatt_tup, patt_tup in enumerate(list(itertools.product(Dates_list,
                                                                AC_names,
                                                                prod_types))):
        dt_cur , ac_cur , prod_cur  = patt_tup
        wwww , dow = conv.dt2gpstime(dt_cur)
        
        #### Manage the cases of manual DOW
        if type(dow_manu) is int:
            dow = dow_manu
        elif dow_manu is None:
            dow = ""
        elif dow_manu is False:
            pass        
        else:
            dow = str(dow_manu)
               
        log.info("*** Search prods. for %s-%s, AC/prod: %s/%s",wwww,dow,ac_cur,prod_cur)
        wwww_dir = os.path.join(arch_center_basedir,str(wwww),repro_str)
        
        n_ftp_ask = 500 ## Max interrogation of the FTP server to avoid 
                        ## potential errors
                        ## An new FTP instance is created if above it 
        
        if np.mod(ipatt_tup,n_ftp_ask) == 0:
            log.info("Create a new FTP instance")
            ftp, Ftp_obj_list = ftp_objt_create(secure_ftp)
            
        if wwww_dir_previous != wwww_dir or np.mod(ipatt_tup,n_ftp_ask) == 0:
            log.info("Move to: %s",wwww_dir)
            try:
                ftp.cwd(wwww_dir)
            except:
                log.warning("%s do not exists, skiping...",wwww_dir)
            Files_listed_in_FTP = ftp.nlst()
            wwww_dir_previous = wwww_dir
            if len(Files_listed_in_FTP) == 0:
                log.warning("no files found in directory %s",wwww_dir)
            
            
        Files_remote_date_list = []

        pattern_old_nam = ac_cur.lower()+".*"+str(wwww)+str(dow)+".*"+prod_cur.lower()+"\..*"
        Files = [f for f in Files_listed_in_FTP if re.search(pattern_old_nam,f)]
        
        pattern_new_nam = ""
        
        Files_new_nam = []
        if new_name_conv: ### search for new name convention
            
            if dow is None:
                log.error("dow == None and search for new name convention, Error ...")
                raise  Exception()
                
            ac_newnam   = ac_cur.upper()

            doy_newnam  = "".join(reversed(conv.dt2doy_year(conv.gpstime2dt(wwww,dow))))
            prod_newnam = prod_cur.upper()
            
            pattern_new_nam = utils.join_improved(".*",ac_newnam,doy_newnam,prod_newnam)
            pattern_new_nam = ".*" + pattern_new_nam + "\..*"

            Files_new_nam   = [f for f in Files_listed_in_FTP if re.search(pattern_new_nam,f)]
        
        
        log.info("Regex : %s %s",pattern_old_nam,pattern_new_nam)        
        Files = Files + Files_new_nam
        
        if len(Files) == 0:
            log.warning("no product found for" + " %s"*len(patt_tup),
                        *patt_tup)
        
        Files_remote_date_list = Files_remote_date_list + Files
            
        ### exclude some pattern
        for negpatt in remove_patterns:
            Files_remote_date_list = [e for e in Files_remote_date_list if not re.search(negpatt,e)]
            
            
        ###################################################################
        ########### Download
            

        archive_dir_specif = dlutils.effective_save_dir_orbit(archive_dir,
                                                              ac_cur,
                                                              dt_cur,
                                                              archtype)
        
        if len(Files_remote_date_list) > 0:
            utils.create_dir(archive_dir_specif)
        
        ### Generation of the Download fct inputs
        Files_remote_date_chunck = utils.chunkIt(Files_remote_date_list,
                                                 parallel_download)
        Downld_tuples_list = []
        Potential_localfiles_list = []

        if ftp_download: ### FTP Download
            for ftpobj , Chunk in zip(Ftp_obj_list,Files_remote_date_chunck):
                for filchunk in Chunk:
                    Potential_localfiles_list.append(os.path.join(archive_dir_specif,filchunk))
                    if parallel_download == 1:
                        Downld_tuples_list.append((ftpobj,filchunk,archive_dir_specif))
                    else:
                        Downld_tuples_list.append((arch_center_main,wwww_dir,
                                                   filchunk,archive_dir_specif))
        else: ### HTTP download
            Downld_tuples_list = itertools.product(["/".join(('ftp://' + arch_center_main,wwww_dir,f)) for f in Files_remote_date_list],[archive_dir_specif])
            [Potential_localfiles_list.append(os.path.join(archive_dir_specif,f)) for f in Files_remote_date_list]
        
        Potential_localfiles_list_all = Potential_localfiles_list_all +  Potential_localfiles_list 

        ### Actual Download
        if ftp_download and parallel_download == 1:
            for tup in Downld_tuples_list:
                dlutils.FTP_downloader(*tup)
        elif ftp_download and parallel_download > 1:
            _ = pool.map_async(dlutils.FTP_downloader_wo_objects,
                               Downld_tuples_list)
        elif not ftp_download and parallel_download == 1:
            for tup in Downld_tuples_list:
                dlutils.downloader_wrap(tup)
        elif not ftp_download and parallel_download > 1:
            _ = pool.map(dlutils.downloader_wrap,Downld_tuples_list)
    

    ###################################################################
    ########### Final Independent files existence check

    Localfiles_lis = []
    if not return_also_uncompressed_files:
        Pot_locfiles_list_use = Potential_localfiles_list_all
    else:
        Pot_locfiles_list_use = []
        for localfile in Potential_localfiles_list_all:
            Pot_compress_name_list = [localfile]
            Pot_compress_name_list.append(localfile.replace(".gz",""))
            Pot_compress_name_list.append(localfile.replace(".Z",""))
            Pot_compress_name_list = list(set(Pot_compress_name_list))
            
            Pot_locfiles_list_use = Pot_locfiles_list_use + Pot_compress_name_list
            
    for pot_localfile in Pot_locfiles_list_use:
        if os.path.isfile(pot_localfile):
            Localfiles_lis.append(pot_localfile)
    
    pool.close()
    return Localfiles_lis


def multi_downloader_orbs_clks_2():
    log.warn('multi_downloader_orbs_clks_2 is a legacy alias for the newly renamed function download_products_gnss')
    #return download_products_gnss(**kwargs)


def orbclk_long2short_name(longname_filepath_in,
                           rm_longname_file=False,
                           center_id_last_letter=None,
                           center_manual_short_name=None,
                           force=False,
                           dryrun=False,
                           output_dirname=None):
    """
    Rename a long naming new convention IGS product file to the short old
    convention
    Naming will be done automaticcaly based on the 3 first charaters of the
    long AC id
    e.g. CODE => cod, GRGS => grg, NOAA => noa ...

    Parameters
    ----------
    longname_filepath_in : str
        Full path of the long name product file

    rm_longname_file : bool
        Remove the original long name product file

    center_id_last_letter : str
        replace the last letter of the short AC id by another letter
        (see note below)
        
    center_manual_short_name : str
        replace completely the long name with this one
        overrides center_id_last_letter
        
    force : bool
        if False, skip if the file already exsists

    dryrun : bool
        if True, don't rename effectively, just output the new name
        
    output_dirname : str
        directory where the output shortname will be created
        if None, will be created in the same folder as the input longname

    Returns
    -------
    shortname_filepath : str
        Path of the  short old-named product file

    Note
    ----
    if you rename MGEX orbits, we advise to set
    center_id_last_letter="m"
    the AC code name will be changed to keep a MGEX convention
    (but any other caracter can be used too)

    e.g. for Bern's products, the long id is CODE

    if center_id_last_letter=None, it will become cod,
    if center_id_last_letter=m, it will become com

    """

    log.info("will rename " + longname_filepath_in)

    longname_basename = os.path.basename(longname_filepath_in)
    longname_dirname  = os.path.dirname(longname_filepath_in)
    
    if not output_dirname:
        output_dirname = longname_dirname

    center = longname_basename[:3]

    if center_manual_short_name:
        center = center_manual_short_name
    elif center_id_last_letter:
        center_as_list = list(center)
        center_as_list[-1] = center_id_last_letter
        center = "".join(center_as_list)

    yyyy   = int(longname_basename.split("_")[1][:4])
    doy    = int(longname_basename.split("_")[1][4:7])

    day_dt = conv.doy2dt(yyyy,doy)

    wwww , dow = conv.dt2gpstime(day_dt)

    shortname_prefix = center.lower() + str(wwww).zfill(4) + str(dow)

    ### Type handeling
    if   "SP3" in longname_basename:
        shortname = shortname_prefix + ".sp3"
    elif "CLK" in longname_basename:
        shortname = shortname_prefix + ".clk"
    elif "ERP" in longname_basename:
        shortname = shortname_prefix + ".erp"
    elif "BIA" in longname_basename:
        shortname = shortname_prefix + ".bia"
    elif "SNX" in longname_basename:
        shortname = shortname_prefix + ".snx"
    else:
        log.error("filetype not found for " + longname_basename)
    
    ### Compression handeling
    if longname_basename[-3:] == ".gz":
        shortname = shortname + ".gz"
    elif longname_basename[-2:] == ".Z":
        shortname = shortname + ".Z"
        
    shortname_filepath = os.path.join(output_dirname , shortname)
    
    if not force and os.path.isfile(shortname_filepath):
        log.info("skip " + longname_filepath_in)
        log.info(shortname_filepath + " already exists")        
        return shortname_filepath
    
    if not dryrun:
        log.info("renaming " + longname_filepath_in + " => " + shortname_filepath)
        shutil.copy2(longname_filepath_in , shortname_filepath)

    if rm_longname_file and not dryrun:
        log.info("remove " + longname_filepath_in)
        os.remove(longname_filepath_in)

    return shortname_filepath


 #  ______                _   _                _____                                         _ 
 # |  ____|              | | (_)              / ____|                                       | |
 # | |__ _   _ _ __   ___| |_ _  ___  _ __   | |  __ _ __ __ ___   _____ _   _  __ _ _ __ __| |
 # |  __| | | | '_ \ / __| __| |/ _ \| '_ \  | | |_ | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
 # | |  | |_| | | | | (__| |_| | (_) | | | | | |__| | | | (_| |\ V /  __/ |_| | (_| | | | (_| |
 # |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \_____|_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
 #                                                                        __/ |                
 #                                                                       |___/ 
 
 
def multi_downloader_orbs_clks(archive_dir,startdate,enddate,calc_center='igs',
                               sp3clk='sp3',archtype ='year/doy',parallel_download=4,
                               archive_center='ign',repro=0,sorted_mode=False,
                               force_weekly_file=False, return_also_uncompressed_files=True):

    """
    Download IGS products. Can manage MGEX products too
    (see archive_center argument)

    Parameters
    ----------
    archive_dir : str
        Parent archive directory where files will be stored

    startdate & enddate : datetime
        Start and End of the wished period

    calc_center : str or list of str
        calc_center can be a string or a list, describing the calc center
        e.g. 'igs','grg','cod','jpl' ...

    sp3clk : str
        Product type, can handle :

            'clk'

            'clk_30s'

            'sp3'
            
            'snx'
            
            'sum'

            'erp'

            'bia'

    archive_center : str
        server of download, "regular" IGS or MGEX, can handle :

            'cddis'

            'cddis_mgex'

            'cddis_mgex_longname'

            'ign'

            'ign_mgex'

            'ign_mgex_longname'
            
            'gfz_local'

    archtype: str
        string describing how the archive directory is structured, e.g :

            stat

            stat/year

            stat/year/doy

            year/doy

            year/stat

            week/dow/stat

            ... etc ...

    repro : int
        number of the IGS reprocessing (0 = routine processing)

    sorted_mode : bool
        if False:
            using the map multiprocess fct so the download order will
            be scrambled
        if True:
            using the apply multiprocess fct so the download order will be
            in the chrono. order
        The scrambled (False) is better, bc. it doesn't create zombies processes

    Returns
    -------
    localfiles_lis : list of str
        list of downloaded products paths
    """
    
    log.error("multi_downloader_orbs_clks IS DISCONTINUED, use multi_downloader_orbs_clks_2")
    return None

#     if type(calc_center) is str:
#         calc_center =  [ ''.join([calc_center]) ] # POURQUOI CETTE LIGNE ?? (150717)
# #        calc_center =  [calc_center]

#     pool = mp.Pool(processes=parallel_download)
#     urllist = []
#     savedirlist = []

#     if sp3clk == 'clk':
#         typ = 'Clocks'
#     elif sp3clk == 'clk_30s':
#         typ = 'Clocks (30s)'
#     elif sp3clk == 'sp3':
#         typ = 'Orbits'
#     elif sp3clk == 'erp':
#         typ = 'ERP'
#     elif sp3clk == 'snx':
#         typ = 'SINEXs'
#     elif sp3clk == 'sum':
#         typ = 'SUM files'
#     elif sp3clk == 'bia':
#         typ = 'ISBs'
#     else:
#         typ = '????'

#     log.info("generating the list of potential " + typ + " ...")
    
#     for cc in calc_center:
#             curdate = startdate
#             while curdate <= enddate:
#                 if re.search("igs([0-9][0-9]|yy|YY)P",cc):
#                     cc = "igs" + str(curdate.year)[2:] + "P"
#                     log.info("INFO : IGS reference frame snx/ssc, correcting the year : " + cc)

#                 url = ''
#                 if archive_center == 'cddis':
#                     url = orbclk_cddis_server(curdate,cc,repro=repro,sp3clk=sp3clk,
#                                               force_weekly_file=force_weekly_file)
#                 elif archive_center == 'cddis_mgex':
#                     url = orbclk_cddis_server(curdate,cc,repro=repro,sp3clk=sp3clk,
#                                               mgex=True,force_weekly_file=force_weekly_file)
#                 elif archive_center == 'cddis_mgex_longname':
#                     url = orbclk_cddis_server(curdate,cc,repro=repro,
#                                               sp3clk=sp3clk,mgex=True,longname=True,
#                                               force_weekly_file=force_weekly_file)
#                 elif archive_center == 'ign':
#                     url = orbclk_ign_server(curdate,cc,repro=repro,sp3clk=sp3clk,
#                                               mgex=False,force_weekly_file=force_weekly_file)
#                 elif archive_center == 'ign_mgex':
#                     url = orbclk_ign_server(curdate,cc,repro=repro,sp3clk=sp3clk,
#                                               mgex=True,force_weekly_file=force_weekly_file)
#                 elif archive_center == 'ign_mgex_longname':
#                     url = orbclk_ign_server(curdate,cc,repro=repro,
#                                               sp3clk=sp3clk,mgex=True,longname=True,
#                                               force_weekly_file=force_weekly_file)
#                 elif archive_center == 'gfz_local':
#                     url = orbclk_gfz_local_server(curdate,cc,repro=repro,
#                                               sp3clk=sp3clk)

#                 else:
#                     log.error('ERR : Wrong archive_center name !!! :' + archive_center)
#                 urllist.append(url)
#                 savedir = dlutils.effective_save_dir_orbit(archive_dir,
#                                                            cc,
#                                                            curdate,
#                                                            archtype)
#                 savedirlist.append(savedir)
#                 curdate = curdate + dt.timedelta(days=1)

#     savedirlist = [x for (y,x) in sorted(zip(urllist,savedirlist))]
#     urllist = sorted(urllist)

#     log.info(" ... done")
#     log.info(str(len(urllist)) + " potential " + typ)

#     if not sorted_mode:
#         _ = pool.map(dlutils.downloader_wrap,list(zip(urllist,savedirlist)))
#     else:
#          results = [pool.apply_async(dlutils.downloader , args=(u,sd)) for u,sd in zip(urllist,savedirlist)]

#     localfiles_lis = []
    
#     if not return_also_uncompressed_files:
#         for url , savedir in zip(urllist,savedirlist):
#             localfile = os.path.join(savedir,os.path.basename(url))
#             if os.path.isfile(localfile):
#                 localfiles_lis.append(localfile)
#     else:
    
#         for url , savedir in zip(urllist,savedirlist):

#             localfile = os.path.join(savedir,os.path.basename(url))        
            
#             Pot_compress_files_list = [localfile]
#             Pot_compress_files_list.append(localfile.replace(".gz",""))
#             Pot_compress_files_list.append(localfile.replace(".Z",""))
#             Pot_compress_files_list = list(set(Pot_compress_files_list))
            
#             for potential_exisiting_file in Pot_compress_files_list:
#                 if os.path.isfile(potential_exisiting_file):
#                     localfiles_lis.append(potential_exisiting_file)
    
#     pool.close()
#     return localfiles_lis



# def force_weekly_file_fct(force_weekly_file,sp3clk,day_in):
#     if force_weekly_file == False:
#         day = day_in

#     elif type(force_weekly_file) is str or type(force_weekly_file) is int:
#         log.info("INFO : The weekly file will be downloaded (DoW = %s)",force_weekly_file)
#         log.info("       Check force_weekly_file option if you don't want it")
#         day = force_weekly_file

#     elif sp3clk in ("erp","sum") and force_weekly_file == True:
#         log.info("The weekly file (DoW = 7) will be downloaded for " + sp3clk.upper())
#         day = 7
        
#     return day


# def orbclk_cddis_server(date,center='igs', sp3clk = 'sp3', repro=0, mgex=False,
#                         longname = False, force_weekly_file=False):
#     """
#     longname is experimental and only for MGEX yet !! (180426)
#     """
#     urlserver='ftp://cddis.gsfc.nasa.gov/pub/gps/products/'
#     if mgex:
#         urlserver = urlserver + 'mgex/'
#     if repro == 0:
#         rep_fldr = ''
#     else:
#         rep_fldr = 'repro' + str(repro)

#     if repro != 0:
#         center     = list(center)
#         center[-1] = str(repro)
#         center = ''.join(center)
#     if repro == 3:
#         longname = True
        
#     if center in ("cod","cof","co2","cf2") and sp3clk == "sp3":
#         log.info("CODE orbit extension changed to eph")
#         sp3clk = "eph"

#     # date definition
#     week, day = conv.dt2gpstime(date)

#     ## force_weekly_file handeling
#     day = force_weekly_file_fct(force_weekly_file_fct,sp3clk,day)

#     if not longname: # e.g. gbm19903.sp3.Z
#         if not 'igu' in center:
#             orbname = center + str(week).zfill(4)  + str(day) +'.'+ sp3clk +'.Z'
#             url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)
#         else:
#             igusuffix = center[-2:]
#             center    = center[:-2]
#             orbname = center + str(week).zfill(4)  + str(day)  + '_' + igusuffix +'.'+ sp3clk + '.Z'
#             url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)
#     else: # e.g. COD0MGXFIN_20180580000_01D_05M_ORB.SP3.gz
#         if len(center) == 3:
#             center = center.upper() + "0"
#         else:
#             center = center.upper()

#         datelong = date.strftime("%Y%j")

#         if "SHA" in center:
#             if sp3clk == "sp3":
#                 sp3clk_long = "_01D_15M_ORB.SP3"
#             elif sp3clk == "clk":
#                 sp3clk_long = "_01D_05M_CLK.CLK"
#             elif sp3clk == "erp":
#                 sp3clk_long = "_03D_12H_ERP.ERP"
#             elif sp3clk == "bia":
#                 sp3clk_long = "_01D_01D_OSB.BIA"
#             elif sp3clk == "snx":
#                 sp3clk_long = "_01D_000_SOL.SNX"

#             orbname = center + "MGXRAP_" + datelong + "0000"  + sp3clk_long + ".gz"
#         else:
#             if sp3clk == "sp3":
#                 sp3clk_long = "_01D_15M_ORB.SP3"
#             elif sp3clk == "clk":
#                 sp3clk_long = "_01D_30S_CLK.CLK"
#             elif sp3clk == "erp":
#                 sp3clk_long = "_03D_12H_ERP.ERP"
#             elif sp3clk == "bia":
#                 sp3clk_long = "_01D_01D_OSB.BIA"
#             elif sp3clk == "snx":
#                 sp3clk_long = "_01D_000_SOL.SNX"

#             orbname = center + "MGXFIN_" + datelong + "0000"  + sp3clk_long + ".gz"

#         url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)

#     return url



# def orbclk_ign_server(date,center='igs', sp3clk = 'sp3', repro=0, mgex=False,
#                         longname = False,force_weekly_file=False):
#     """
#     longname is experimental and only for MGEX yet !! (180426)
#     Must be merged with orbclk_cddis_server to make a big equivalent fct (180523)
#     """
#     urlserver='ftp://igs.ign.fr/pub/igs/products/'
#     if mgex:
#         urlserver = urlserver + 'mgex/'
#     if repro == 0:
#         rep_fldr = ''
#     else:
#         rep_fldr = 'repro' + str(repro)
#     if repro != 0:
#         center     = list(center)
#         center[-1] = str(repro)
#         center = ''.join(center)
#     week, day = conv.dt2gpstime(date)

#     ## force_weekly_file handeling
#     day = force_weekly_file_fct(force_weekly_file,sp3clk,day)

#     if not longname: # e.g. gbm19903.sp3.Z
#         if not 'igu' in center:
#             orbname = center + str(week).zfill(4)  + str(day) +'.'+ sp3clk +'.Z'
#             url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)
#         else:
#             igusuffix = center[-2:]
#             center    = center[:-2]
#             orbname = center + str(week).zfill(4)  + str(day)  + '_' + igusuffix +'.'+ sp3clk + '.Z'
#             url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)
#     else: # e.g. COD0MGXFIN_20180580000_01D_05M_ORB.SP3.gz
#         if len(center) == 3:
#             center = center.upper() + "0"
#         else:
#             center = center.upper()

#         datelong = date.strftime("%Y%j")

#         if "SHA" in center:
#             if sp3clk == "sp3":
#                 sp3clk_long = "_01D_15M_ORB.SP3"
#             elif sp3clk == "clk":
#                 sp3clk_long = "_01D_05M_CLK.CLK"
#             elif sp3clk == "erp":
#                 sp3clk_long = "_03D_12H_ERP.ERP"
#             elif sp3clk == "bia":
#                 sp3clk_long = "_01D_01D_OSB.BIA"

#             orbname = center + "MGXRAP_" + datelong + "0000"  + sp3clk_long + ".gz"
#         else:
#             if sp3clk == "sp3":
#                 sp3clk_long = "_01D_15M_ORB.SP3"
#             elif sp3clk == "clk":
#                 sp3clk_long = "_01D_30S_CLK.CLK"
#             elif sp3clk == "erp":
#                 sp3clk_long = "_03D_12H_ERP.ERP"
#             elif sp3clk == "bia":
#                 sp3clk_long = "_01D_01D_OSB.BIA"

#             orbname = center + "MGXFIN_" + datelong + "0000"  + sp3clk_long + ".gz"

#         url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)

#     return url

# def orbclk_igscb_server(date,center='gfz', sp3clk = 'sp3',repro=0):
#     urlserver='ftp://igscb.jpl.nasa.gov/pub/product/'
#     if repro == 0:
#         rep_fldr = ''
#     elif repro == 3:
#         rep_fldr = 'repro' + str(repro)
#     if repro != 0:
#         center     = list(center)
#         center[-1] = str(repro)
#         center = ''.join(center)
#     week, day = conv.dt2gpstime(date)
#     if not 'igu' in center:
#         orbname = center + str(week).zfill(4) + str(day) +'.' + sp3clk + '.Z'
#         url = os.path.join(urlserver,str(week).zfill(4) ,rep_fldr,orbname)
#     else:
#         igusuffix = center[-2:]
#         center    = center[:-2]
#         orbname = center + str(week).zfill(4) + str(day)  + '_' + igusuffix +'.' + sp3clk + '.Z'
#         url = os.path.join(urlserver,str(week).zfill(4) ,rep_fldr,orbname)
#     return url

# def orbclk_gfz_local_server(date,center='gfz', sp3clk='sp3',repro=0):
#     if repro == 0:
#         urlserver = '/dsk/igs_archive/IGS/SAVE_PROD_1d/'
#     elif repro == 3:
#         urlserver = '/dsk/repro3/ARCHIVE/IGS/SAVE_PROD_1d/'
#     else:
#         log.error("check the repro !!!")
#         raise Exception

#     week, day = conv.dt2gpstime(date)
    
#     orbname = center + str(week).zfill(4) + str(day) +'.' + sp3clk + '.Z'
#     url = os.path.join(urlserver,str(week).zfill(4),orbname)
    
#     return url




