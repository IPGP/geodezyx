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
# import datetime as dt
from ftplib import FTP
# import glob
# import itertools
# import multiprocessing as mp
import os
import pandas as pd 
# import re
import shutil
import urllib
import ftplib

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

#### Import star style
# from geodezyx import *                   # Import the GeodeZYX modules
# from geodezyx.externlib import *         # Import the external modules
# from geodezyx.megalib.megalib import *   # Import the legacy modules names


#### Import the logger
import logging
log = logging.getLogger(__name__)

##########  END IMPORT  ##########


def start_end_date_easy(start_year,start_doy,end_year,end_doy):
    """
    generates start/end datetimes from a start/end year/day of year

    Parameters
    ----------
    start_year : int
        start year.
    start_doy : int
        start day of year.
    end_year : int
        end year.
    end_doy : int
        end day of year.

    Returns
    -------
    start : datetime
        converted start datetime.
    end : datetime
        converted end datetime.

    """
    start = conv.doy2dt(start_year,start_doy)
    end   = conv.doy2dt(end_year,end_doy)
    return start , end

def effective_save_dir_orbit(parent_archive_dir,
                             calc_center,date,
                             archtype ='year/doy/'):
    """
    INTERNAL_FUNCTION

    archtype =
        stat
        stat/year
        stat/year/doy
        year/doy
        year/stat
        week/dow
        wkwwww : use a GFZ's CF-ORB wk<wwww> naming
        OR only '/' for a dirty saving in the parent folder
        ... etc ...
    """
    if archtype == '/':
        return parent_archive_dir

    out_save_dir = parent_archive_dir
    fff = archtype.split('/')
    #year = str(date.year)
    #doy = conv.dt2doy(date)
    week, dow = conv.dt2gpstime(date)

    for f in fff:
        if "wkwwww" in f:
            f_evaluated = "wk" + str(week).zfill(4)
        else:
            f_evaluated = eval(f)
        out_save_dir = os.path.join(out_save_dir,str(f_evaluated))
    return out_save_dir



#  _    _ _______ _______ _____    _____                      _                 _ 
# | |  | |__   __|__   __|  __ \  |  __ \                    | |               | |
# | |__| |  | |     | |  | |__) | | |  | | _____      ___ __ | | ___   __ _  __| |
# |  __  |  | |     | |  |  ___/  | |  | |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |
# | |  | |  | |     | |  | |      | |__| | (_) \ V  V /| | | | | (_) | (_| | (_| |
# |_|  |_|  |_|     |_|  |_|      |_____/ \___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|
                                                                                
                                                                                


#### HTTP classic Download

def downloader(url,savedir,force = False,
               check_if_file_already_exists_uncompressed=True):
    """
    general function to download a file

    INTERNAL_FUNCTION
    """

    if type(url) is tuple:
        need_auth = True
        username = url[1]
        password = url[2]
        url = url[0]
    else:
        need_auth = False

    url_print = str(url)

    rnxname = os.path.basename(url)

    Pot_compress_files_list = [os.path.join(savedir , rnxname)]

    if check_if_file_already_exists_uncompressed:
        Pot_compress_files_list.append(os.path.join(savedir , rnxname.replace(".gz","")))
        Pot_compress_files_list.append(os.path.join(savedir , rnxname.replace(".Z","")))
        Pot_compress_files_list = list(set(Pot_compress_files_list))

    for f in Pot_compress_files_list:
        if os.path.isfile(f) and (not force):
            log.info(os.path.basename(f) + " already exists locally ;)")
            return None
        
    ##### LOCAL FILE (particular case for GFZ)
    if os.path.isfile(url):
        log.info("INFO : downloader : the is a local file, a simple copy will be used")
        log.info("       URL : %s",url)
        shutil.copy(url,savedir)
    
    ##### REMOTE FILE (General case)
    elif ("ftp" in url) or ("http" in url) :
        # managing a authentification
        if need_auth:
            if 'ftp://' in url: # FTP with Auth
                url = url.replace('ftp://','ftp://' + username + ':' + password + '@')
                opener = urllib.request.build_opener()
            else: # HTTP with Auth
                password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
                top_level_url = url
                password_mgr.add_password(None, top_level_url, username, password)
                handler = urllib.request.HTTPBasicAuthHandler(password_mgr)
                # create "opener" (OpenerDirector instance)
                opener = urllib.request.build_opener(handler)
        else: # FTP or HTTP without Auth
            opener = urllib.request.build_opener()
    
        # use the opener to fetch a URL
        try:
            f = opener.open(url)
        except (urllib.error.HTTPError , urllib.error.URLError):
            log.warning("%s not found :(",rnxname)
            log.warning(url_print)
            return ""
        log.info("%s found, downloading :)",rnxname)
        data = f.read()
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        outpath = os.path.join(savedir , rnxname)
        with open(outpath, "wb") as code:
            code.write(data)
        return_str = outpath
    else:
        log.error("something goes wrong with the URL")
        log.error(url)
        return_str = ""

    return return_str


def downloader_wrap(intup):
    downloader(*intup)
    return None

 
#  ______ _______ _____    _____                      _                 _ 
# |  ____|__   __|  __ \  |  __ \                    | |               | |
# | |__     | |  | |__) | | |  | | _____      ___ __ | | ___   __ _  __| |
# |  __|    | |  |  ___/  | |  | |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |
# | |       | |  | |      | |__| | (_) \ V  V /| | | | | (_) | (_| | (_| |
# |_|       |_|  |_|      |_____/ \___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|
                                                                        
 
#### FTP DOWNLOAD                                                                                      
def _ftp_dir_list_files(ftp_obj_in):
    files = []
    try:
        files = ftp_obj_in.nlst()
    except ftplib.error_perm as resp:
        if str(resp) == "550 No files found":
            log.warning("No files in this directory" + ftp_obj_in.pwd())
        else:
            raise
    return files


def ftp_files_crawler(urllist,savedirlist):
    """
    filter urllist,savedirlist generated with multi_downloader_rinex with an
    optimized FTP crawl

    """
    ### create a DataFrame based on the urllist and savedirlist lists
    DF = pd.concat((pd.DataFrame(urllist),pd.DataFrame(savedirlist)),axis=1)
    DF_orig = DF.copy()
    
    ### rename the columns
    if DF.shape[1] == 4:
        loginftp = True 
        DF.columns = ('url','user','pass','savedir')
    else:
        loginftp = False
        DF.columns = ('url','savedir')
        DF['user'] = "anonymous"
        DF['pass'] = ""
        
    ### Do the correct split for the URLs
    DF = DF.sort_values('url')
    DF["url"]      = DF['url'].str.replace("ftp://","")
    DF["dirname"]  = DF["url"].apply(os.path.dirname)
    DF["basename"] = DF["url"].apply(os.path.basename)
    DF["root"]     = [e.split("/")[0] for e in DF["dirname"].values]
    DF["dir"]      = [e1.replace(e2,"")[1:] for (e1,e2) in zip(DF["dirname"],
                                                                   DF["root"])]
    DF["bool"] = False
    
    #### Initialisation of the 1st variables for the loop
    prev_row_ftpobj = DF.iloc[0]
    prev_row_cwd    = DF.iloc[0]
    FTP_files_list = []
    count_loop = 0  # restablish the connexion after 50 loops (avoid freezing)
    #### Initialisation of the FTP object
    FTPobj = ftplib.FTP(prev_row_ftpobj.root,
                        prev_row_ftpobj.user, 
                        prev_row_ftpobj['pass'])
    
    for irow,row in DF.iterrows():
        count_loop += 1
        
        ####### we recreate a new FTP object if the root URL is not the same
        if row.root != prev_row_ftpobj.root or count_loop > 20:
    
            FTPobj = ftplib.FTP(row.root,
                                row.user, 
                                row['pass'])
            
            prev_row_ftpobj = row
            count_loop = 0
            
        ####### we recreate a new file list if the date path is not the same        
        if (prev_row_cwd.dir != row.dir) or irow == 0:
            log.info("chdir " + row.dirname)
            FTPobj.cwd("/")
    
            try: #### we try to change for the right folder
                FTPobj.cwd(row.dir)
            except: #### If not possible, then no file in the list
                FTP_files_list = []
              
            FTP_files_list = _ftp_dir_list_files(FTPobj)
            prev_row_cwd = row 
               
        ####### we check if the files is avaiable
        if row.basename in FTP_files_list:
            DF.loc[irow,'bool'] = True
            log.info(row.basename + " found on server :)")
        else:
            DF.loc[irow,'bool'] = False
            log.warning(row.basename  + " not found on server :(")
    
    
    DFgood = DF[DF['bool']].copy()
    
    DFgood['url'] = 'ftp://' + DFgood['url']  
    
    ### generate the outputs
    if loginftp:
        urllist_out = list(zip(DFgood.url,DFgood.user,DFgood['pass']))
    else:
        urllist_out = list(DFgood.url)    
    
    savedirlist_out = list(DFgood.savedir)
    
    return urllist_out, savedirlist_out

def FTP_downloader(ftp_obj,filename,localdir):   
    localpath = os.path.join(localdir,filename)
    
    if not utils.empty_file_check(localpath):
        log.info(filename + " already exists ;)")
        bool_dl = True
    else:
        try:
            localfile = open(localpath, 'wb')
            ftp_obj.retrbinary('RETR ' + filename, localfile.write, 1024)
            localfile.close()
            bool_dl = True
            log.info(filename + " downloaded :)")
    
        except:
            log.info(localpath + " download failed :(")
            bool_dl = False
    
    return localpath , bool_dl

def FTP_downloader_wo_objects(tupin):
    arch_center_main,wwww_dir,filename,localdir = tupin
    ftp_obj_wk = FTP(arch_center_main)
    ftp_obj_wk.login()
    ftp_obj_wk.cwd(wwww_dir)
    localpath , bool_dl = FTP_downloader(ftp_obj_wk,filename,localdir)
    ftp_obj_wk.close()
    return localpath , bool_dl
    
