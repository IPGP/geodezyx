#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.operational contains functions to download
gnss data and products from distant IGS servers. 

it can be imported directly with:
from geodezyx import operational

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU GPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (GFZ, pierre.sakic@gfz-postdam.de)
GitHub repository :
https://github.com/GeodeZYX/GeodeZYX-Toolbox_v4
"""

#_____________ _   _ ________   __  _____                      _                 _
#|  __ \|_   _| \ | |  ____\ \ / / |  __ \                    | |               | |
#| |__) | | | |  \| | |__   \ V /  | |  | | _____      ___ __ | | ___   __ _  __| | ___ _ __
#|  _  /  | | | . ` |  __|   > <   | |  | |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |/ _ \ '__|
#| | \ \ _| |_| |\  | |____ / . \  | |__| | (_) \ V  V /| | | | | (_) | (_| | (_| |  __/ |
#|_|  \_\_____|_| \_|______/_/ \_\ |_____/ \___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|\___|_|






########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
from ftplib import FTP
import glob
import itertools
import multiprocessing as mp
import os 
import re
import shutil
import urllib

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

##########  END IMPORT  ##########



############################################################################
######## RINEX DOWNLOADER
############################################################################


def igs_garner_server(stat,date):
    # plante si trop de requete
    urlserver = "ftp://garner.ucsd.edu/pub/rinex/"
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join( urlserver , str(date.year) , conv.dt2doy(date) , rnxname )
    return url

def igs_cddis_server(stat,date):
    # a privilegier
    urlserver = "ftp://cddis.gsfc.nasa.gov/gps/data/daily/"
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , conv.dt2doy(date) , date.strftime('%y') + 'd' , rnxname)
    return url

def igs_cddis_nav_server(stat,date):
    # a privilegier
    urlserver = "ftp://cddis.gsfc.nasa.gov/gps/data/daily/"
    rnxname = conv.statname_dt2rinexname(stat.lower(),date,'n.Z')
    url = os.path.join(urlserver , str(date.year) , conv.dt2doy(date) , date.strftime('%y') + 'n' , rnxname)
    return url

def rgp_ign_smn_server(stat,date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , conv.dt2doy(date) , 'data_30' , rnxname)
    return url

def rgp_ign_smn_1Hz_server(stat,date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"

    urls = []

    for h in range(24):
        date_session = date
        date_session = date_session.replace(hour=h)

        print(date_session , 'session' , h)
        rnxname = conv.statname_dt2rinexname(stat.lower(),date_session ,
                                             session_a_instead_of_daily_session = 1)
        url = os.path.join(urlserver , str(date.year) , conv.dt2doy(date) ,
                           'data_1' , rnxname)

        urls.append(url)

    return urls

def unavco_server(stat,date):
    urlserver='ftp://data-out.unavco.org/pub/rinex'
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , 'obs', str(date.year) , conv.dt2doy(date) , rnxname)
    return url

def renag_server(stat,date):
    urlserver = "ftp://renag.unice.fr/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , conv.dt2doy(date) , rnxname)
    return url

def orpheon_server(stat,date,user='',passwd=''):
    urlserver = "ftp://renag.unice.fr/"
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , conv.dt2doy(date) , rnxname)
    return url,user,passwd

def ovsg_server(stat,date,user='',passwd=''):
    if dt.datetime(2009,1,1) <= date <= dt.datetime(2014,2,10):
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS-GPSDATA.backtemp_20140210/"
    else:
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS/GPSDATA/"
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , conv.dt2doy(date) , 'rinex' , rnxname)
    return url,user,passwd

def geoaus_server(stat,date):
    """ Geosciences Australia
        ex : ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/2010/10063/ """
    urlserver = "ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/"
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , date.strftime('%y') + conv.dt2doy(date) , rnxname)
    return url

def sonel_server(stat,date):
    """ex : ftp://ftp.sonel.org/gps/data/2015/001/ """
    urlserver = 'ftp://ftp.sonel.org/gps/data/'
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver,str(date.year),conv.dt2doy(date),rnxname)
    return url

def effective_save_dir(parent_archive_dir,stat,date,archtype ='stat'):
    """
    INTERNAL_FUNCTION

    archtype =
        stat
        stat/year
        stat/year/doy
        year/doy
        year/stat
        week/dow
        OR only '/' for a dirty saving in the parent folder
        ... etc ... """
    if archtype == '/':
        return parent_archive_dir

    out_save_dir = parent_archive_dir
    fff = archtype.split('/')
    year = str(date.year)
    doy = conv.dt2doy(date)
    week, dow = conv.dt2gpstime(date)
    for f in fff:
        out_save_dir = os.path.join(out_save_dir,eval(f))
    return out_save_dir


############################################################################
######## PRODUCTS DOWNLOADER
############################################################################


def force_weekly_file_fct(force_weekly_file,sp3clk,day_in):
    if force_weekly_file == False:
        day = day_in

    elif type(force_weekly_file) is str or type(force_weekly_file) is int:
        print("INFO : The weekly file will be downloaded (DoW =",force_weekly_file,")")
        print("       Check force_weekly_file option if you don't want it")
        day = force_weekly_file

    elif sp3clk in ("erp","sum") and force_weekly_file == True:
        print("INFO : The weekly file (DoW = 7) will be downloaded for " + sp3clk.upper())
        day = 7
        
    return day




def orbclk_cddis_server(date,center='igs', sp3clk = 'sp3', repro=0, mgex=False,
                        longname = False, force_weekly_file=False):
    """
    longname is experimental and only for MGEX yet !! (180426)
    """
    urlserver='ftp://cddis.gsfc.nasa.gov/pub/gps/products/'
    if mgex:
        urlserver = urlserver + 'mgex/'
    if repro == 0:
        rep_fldr = ''
    else:
        rep_fldr = 'repro' + str(repro)
    if repro != 0:
        center     = list(center)
        center[-1] = str(repro)
        center = ''.join(center)
    if center in ("cod","cof","co2","cf2") and sp3clk == "sp3":
        print("INFO : CODE orbit extension changed to eph")
        sp3clk = "eph"

    # date definition
    week, day = conv.dt2gpstime(date)

    ## force_weekly_file handeling
    day = force_weekly_file_fct(force_weekly_file_fct,sp3clk,day)

    if not longname: # e.g. gbm19903.sp3.Z
        if not 'igu' in center:
            orbname = center + str(week).zfill(4)  + str(day) +'.'+ sp3clk +'.Z'
            url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)
        else:
            igusuffix = center[-2:]
            center    = center[:-2]
            orbname = center + str(week).zfill(4)  + str(day)  + '_' + igusuffix +'.'+ sp3clk + '.Z'
            url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)
    else: # e.g. COD0MGXFIN_20180580000_01D_05M_ORB.SP3.gz
        if len(center) == 3:
            center = center.upper() + "0"
        else:
            center = center.upper()

        datelong = date.strftime("%Y%j")

        if "SHA" in center:
            if sp3clk == "sp3":
                sp3clk_long = "_01D_15M_ORB.SP3"
            elif sp3clk == "clk":
                sp3clk_long = "_01D_05M_CLK.CLK"
            elif sp3clk == "erp":
                sp3clk_long = "_03D_12H_ERP.ERP"
            elif sp3clk == "bia":
                sp3clk_long = "_01D_01D_OSB.BIA"
            elif sp3clk == "snx":
                sp3clk_long = "_01D_000_SOL.SNX"

            orbname = center + "MGXRAP_" + datelong + "0000"  + sp3clk_long + ".gz"
        else:
            if sp3clk == "sp3":
                sp3clk_long = "_01D_15M_ORB.SP3"
            elif sp3clk == "clk":
                sp3clk_long = "_01D_30S_CLK.CLK"
            elif sp3clk == "erp":
                sp3clk_long = "_03D_12H_ERP.ERP"
            elif sp3clk == "bia":
                sp3clk_long = "_01D_01D_OSB.BIA"
            elif sp3clk == "snx":
                sp3clk_long = "_01D_000_SOL.SNX"

            orbname = center + "MGXFIN_" + datelong + "0000"  + sp3clk_long + ".gz"

        url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)

    return url



def orbclk_ign_server(date,center='igs', sp3clk = 'sp3', repro=0, mgex=False,
                        longname = False,force_weekly_file=False):
    """
    longname is experimental and only for MGEX yet !! (180426)
    Must be merged with orbclk_cddis_server to make a big equivalent fct (180523)
    """
    urlserver='ftp://igs.ign.fr/pub/igs/products/'
    if mgex:
        urlserver = urlserver + 'mgex/'
    if repro == 0:
        rep_fldr = ''
    else:
        rep_fldr = 'repro' + str(repro)
    if repro != 0:
        center     = list(center)
        center[-1] = str(repro)
        center = ''.join(center)
    week, day = conv.dt2gpstime(date)

    ## force_weekly_file handeling
    day = force_weekly_file_fct(force_weekly_file,sp3clk,day)

    if not longname: # e.g. gbm19903.sp3.Z
        if not 'igu' in center:
            orbname = center + str(week).zfill(4)  + str(day) +'.'+ sp3clk +'.Z'
            url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)
        else:
            igusuffix = center[-2:]
            center    = center[:-2]
            orbname = center + str(week).zfill(4)  + str(day)  + '_' + igusuffix +'.'+ sp3clk + '.Z'
            url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)
    else: # e.g. COD0MGXFIN_20180580000_01D_05M_ORB.SP3.gz
        if len(center) == 3:
            center = center.upper() + "0"
        else:
            center = center.upper()

        datelong = date.strftime("%Y%j")

        if "SHA" in center:
            if sp3clk == "sp3":
                sp3clk_long = "_01D_15M_ORB.SP3"
            elif sp3clk == "clk":
                sp3clk_long = "_01D_05M_CLK.CLK"
            elif sp3clk == "erp":
                sp3clk_long = "_03D_12H_ERP.ERP"
            elif sp3clk == "bia":
                sp3clk_long = "_01D_01D_OSB.BIA"

            orbname = center + "MGXRAP_" + datelong + "0000"  + sp3clk_long + ".gz"
        else:
            if sp3clk == "sp3":
                sp3clk_long = "_01D_15M_ORB.SP3"
            elif sp3clk == "clk":
                sp3clk_long = "_01D_30S_CLK.CLK"
            elif sp3clk == "erp":
                sp3clk_long = "_03D_12H_ERP.ERP"
            elif sp3clk == "bia":
                sp3clk_long = "_01D_01D_OSB.BIA"

            orbname = center + "MGXFIN_" + datelong + "0000"  + sp3clk_long + ".gz"

        url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)

    return url

def orbclk_igscb_server(date,center='gfz', sp3clk = 'sp3',repro=0):
    urlserver='ftp://igscb.jpl.nasa.gov/pub/product/'
    if repro == 0:
        rep_fldr = ''
    elif repro == 3:
        rep_fldr = 'repro' + str(repro)
    if repro != 0:
        center     = list(center)
        center[-1] = str(repro)
        center = ''.join(center)
    week, day = conv.dt2gpstime(date)
    if not 'igu' in center:
        orbname = center + str(week).zfill(4) + str(day) +'.' + sp3clk + '.Z'
        url = os.path.join(urlserver,str(week).zfill(4) ,rep_fldr,orbname)
    else:
        igusuffix = center[-2:]
        center    = center[:-2]
        orbname = center + str(week).zfill(4) + str(day)  + '_' + igusuffix +'.' + sp3clk + '.Z'
        url = os.path.join(urlserver,str(week).zfill(4) ,rep_fldr,orbname)
    return url

def orbclk_gfz_local_server(date,center='gfz', sp3clk='sp3',repro=0):
    if repro == 0:
        urlserver = '/dsk/igs_archive/IGS/SAVE_PROD_1d/'
    elif repro == 3:
        urlserver = '/dsk/repro3/ARCHIVE/IGS/SAVE_PROD_1d/'
    else:
        print("ERR : check the repro !!!")
        raise Exception

    week, day = conv.dt2gpstime(date)
    
    orbname = center + str(week).zfill(4) + str(day) +'.' + sp3clk + '.Z'
    url = os.path.join(urlserver,str(week).zfill(4),orbname)
    
    return url

def downloader(url,savedir,force = False,
               check_if_file_already_exists_uncompressed=True):
    """
    general function to download a file

    INTERNAL_FUNCTION
    """
#    print url
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
            print("INFO :", os.path.basename(f) , "already exists locally ;)")
            return None
        
    ##### LOCAL FILE (particular case for GFZ)
    if os.path.isfile(url):
        print("INFO : downloader : the is a local file, a simple copy will be used")
        print("       URL :",url)
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
            print("WARN :",rnxname,"not found on server :(")
            print(url_print)
            return ""
        print("INFO :" , rnxname ," found on server :)")
        data = f.read()
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        outpath = os.path.join(savedir , rnxname)
        with open(outpath, "wb") as code:
            code.write(data)
        return_str = outpath
    else:
        print("ERR : something goes wrong with the URL")
        print("     ", url)
        return_str = ""


    # effective downloading (old version)
#    try:
#        f = urllib2.urlopen(url)
#    except urllib2.HTTPError:
#        print rnxname," not found :("
#        print url
#        return None
#    print rnxname," found :)"
#    data = f.read()
#    if not os.path.exists(savedir):
#        os.makedirs(savedir)
#    with open(os.path.join(savedir , rnxname), "wb") as code:
#        code.write(data)
    return return_str

def start_end_date_easy(start_year,start_doy,end_year,end_doy):
    start = conv.doy2dt(start_year,start_doy)
    end   = conv.doy2dt(end_year,end_doy)
    return start , end



def effective_save_dir_orbit(parent_archive_dir,calc_center,date,
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
    year = str(date.year)
    doy = conv.dt2doy(date)
    week, dow = conv.dt2gpstime(date)

    for f in fff:
        if "wkwwww" in f:
            f_evaluated = "wk" + str(week).zfill(4)
        else:
            f_evaluated = eval(f)
        out_save_dir = os.path.join(out_save_dir,str(f_evaluated))
    return out_save_dir


def downloader_wrap(intup):
    downloader(*intup)
    return None

def multi_downloader_rinex(statdico,archive_dir,startdate,enddate,
                           archtype ='stat',parallel_download=4,
                           sorted_mode=False,user='',passwd=''):
    """
    Parameters
    ----------
    statdico : dict
        a statdico is a dictionary associating Archives Centers to list of stations

        Exemple:
            >>> statdico['archive center 1'] = ['STA1','STA2','STA3', ...]
            >>> statdico['archive center 2'] = ['STA2','STA1','STA4', ...]

        the supported archive center are (july 2015):
            igs (cddis center)

            igs_garner (for the garner center, but not very reliable)

            rgp (St Mandé center)

            rgp_1Hz (all the 24 hourly rinex for the day will be downloaded)

            renag

            ovsg

            unavco

            sonel

            geoaus (Geosciences Australia)

            nav or brdc as archive center allows to download nav files (using 'BRDC' as station name) from the CDDIS server


    archtype : str
        string describing how the archive directory is structured, e.g :

            stat

            stat/year

            stat/year/doy

            year/doy

            year/stat

            week/dow/stat

            ... etc ...


    sorted_mode : bool
        if False:
            using the map multiprocess fct so the download order will
            be scrambled
        if True:
            using the apply multiprocess fct so the download order will be
            in the chrono. order
        The scrambled (False) is better, bc. it doesn't create zombies processes


    user & passwd : str
        user & password for a locked server


    Returns
    -------
    url_list : list of str
        list of URLs

    savedir_list : list of str
        list of downloaded products paths
    """

    curdate = startdate

    pool = mp.Pool(processes=parallel_download)

    urllist = []
    savedirlist = []


    print("generating the list of potential RINEXs ...")
    while curdate <= enddate:
        for netwk , statlis in list(statdico.items()):
            for stat in statlis:
                stat = stat.lower()
                mode1Hz = False

                if netwk == 'igs':
                    url = igs_cddis_server(stat,curdate)
                elif netwk == 'igs_garner':
                    url = igs_garner_server(stat,curdate)
                elif netwk == 'rgp':
                    url = rgp_ign_smn_server(stat,curdate)
                elif netwk == 'rgp_1Hz':
                    urls = rgp_ign_smn_1Hz_server(stat,curdate)
                    mode1Hz = True
                elif netwk == 'renag':
                    url = renag_server(stat,curdate)
                elif netwk == 'orpheon':
                    url = orpheon_server(stat,curdate,user,passwd)
                elif netwk == 'ovsg':
                    url = ovsg_server(stat,curdate,user,passwd)
                elif netwk == 'unavco':
                    url = unavco_server(stat,curdate)
                elif netwk == 'sonel':
                    url = sonel_server(stat,curdate)
                elif netwk == 'geoaus':
                    url = geoaus_server(stat,curdate)
                elif netwk in ('nav' , 'brdc'):
                    url = igs_cddis_nav_server(stat,curdate)
                else:
                    print('WARN : unkwn server dic in the dico, skip ...')
                    continue

                savedir = effective_save_dir(archive_dir,stat,curdate,archtype)

                if not mode1Hz:
                    urllist.append(url)
                    savedirlist.append(savedir)
                else:
                    urllist = urllist + urls
                    savedirlist = savedirlist + [savedir] * len(urls)

        curdate = curdate + dt.timedelta(days=1)

    #savedirlist = [x for (y,x) in sorted(zip(urllist,savedirlist))]
    #urllist     = sorted(urllist)
    
    urllist,savedirlist = utils.sort_binom_list(urllist,savedirlist)

    print(" ... done")
    print(len(urllist),"potential RINEXs")

    if sorted_mode:
        results = [pool.apply_async(downloader , args=(u,sd)) for u,sd in zip(urllist,savedirlist)]
    else:
        _ = pool.map(downloader_wrap,list(zip(urllist,savedirlist)))


    localfiles_lis = []
    skiped_url = 0
    for url , savedir in zip(urllist,savedirlist):
        try:
            localfile = os.path.join(savedir,os.path.basename(url))
            if os.path.isfile(localfile):
                localfiles_lis.append(localfile)
        except:
            # because of a weird error
            # i = p.rfind('/') + 1
            # AttributeError: 'tuple' object has no attribute 'rfind'
            skiped_url += 1
            continue

    pool.close()
    print('INFO : ' , skiped_url , 'returned url skiped because of a weird error, but it is not important :)')
    return localfiles_lis

#    return zip(urllist,savedirlist)


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

    print("INFO : will rename" , longname_filepath_in)

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
        print("ERR : filetype not found for",longname_basename)
    
    ### Compression handeling
    if longname_basename[-3:] == ".gz":
        shortname = shortname + ".gz"
    elif longname_basename[-2:] == ".Z":
        shortname = shortname + ".Z"
        
    shortname_filepath = os.path.join(output_dirname , shortname)
    
    if not force and os.path.isfile(shortname_filepath):
        print("INFO : skip", longname_filepath_in)
        print("     ",shortname_filepath,"already exists")        
        return shortname_filepath
    
    if not dryrun:
        print("INFO : renaming" , longname_filepath_in,"=>",shortname_filepath)
        shutil.copy2(longname_filepath_in , shortname_filepath)

    if rm_longname_file and not dryrun:
        print("INFO : remove " , longname_filepath_in)
        os.remove(longname_filepath_in)

    return shortname_filepath



def rnx_long2short_name(longname_filepath_in):
    """
    MUST BE IMPROVED
    """
    
    longname_basename = os.path.basename(longname_filepath_in)
    longname_dirname  = os.path.dirname(longname_filepath_in)
    
    Longname_basename_splitted = longname_basename.split("_")
    
    datepart_str = Longname_basename_splitted[2]
    yyyy = datepart_str[:4]
    ddd  = datepart_str[4:7]

    shortname_basename = longname_basename[:4].lower() + ddd + "0." + yyyy[2:] + "o"
    
    return os.path.join(longname_dirname,shortname_basename)
    
    
    


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

    if type(calc_center) is str:
        calc_center =  [ ''.join([calc_center]) ] # POURQUOI CETTE LIGNE ?? (150717)
#        calc_center =  [calc_center]

    pool = mp.Pool(processes=parallel_download)
    urllist = []
    savedirlist = []

    if sp3clk == 'clk':
        typ = 'Clocks'
    elif sp3clk == 'clk_30s':
        typ = 'Clocks (30s)'
    elif sp3clk == 'sp3':
        typ = 'Orbits'
    elif sp3clk == 'erp':
        typ = 'ERP'
    elif sp3clk == 'snx':
        typ = 'SINEXs'
    elif sp3clk == 'sum':
        typ = 'SUM files'
    elif sp3clk == 'bia':
        typ = 'ISBs'
    else:
        typ = '????'

    print("generating the list of potential " + typ + " ...")
    
    for cc in calc_center:
            curdate = startdate
            while curdate <= enddate:
                if re.search("igs([0-9][0-9]|yy|YY)P",cc):
                    cc = "igs" + str(curdate.year)[2:] + "P"
                    print("INFO : IGS reference frame snx/ssc, correcting the year : ",cc)

                url = ''
                if archive_center == 'cddis':
                    url = orbclk_cddis_server(curdate,cc,repro=repro,sp3clk=sp3clk,
                                              force_weekly_file=force_weekly_file)
                elif archive_center == 'cddis_mgex':
                    url = orbclk_cddis_server(curdate,cc,repro=repro,sp3clk=sp3clk,
                                              mgex=True,force_weekly_file=force_weekly_file)
                elif archive_center == 'cddis_mgex_longname':
                    url = orbclk_cddis_server(curdate,cc,repro=repro,
                                              sp3clk=sp3clk,mgex=True,longname=True,
                                              force_weekly_file=force_weekly_file)
                elif archive_center == 'ign':
                    url = orbclk_ign_server(curdate,cc,repro=repro,sp3clk=sp3clk,
                                              mgex=False,force_weekly_file=force_weekly_file)
                elif archive_center == 'ign_mgex':
                    url = orbclk_ign_server(curdate,cc,repro=repro,sp3clk=sp3clk,
                                              mgex=True,force_weekly_file=force_weekly_file)
                elif archive_center == 'ign_mgex_longname':
                    url = orbclk_ign_server(curdate,cc,repro=repro,
                                              sp3clk=sp3clk,mgex=True,longname=True,
                                              force_weekly_file=force_weekly_file)
                elif archive_center == 'gfz_local':
                    url = orbclk_gfz_local_server(curdate,cc,repro=repro,
                                              sp3clk=sp3clk)

                else:
                    print('ERR : Wrong archive_center name !!! :' , archive_center)
                urllist.append(url)
                savedir = effective_save_dir_orbit(archive_dir,cc,curdate,archtype)
                savedirlist.append(savedir)
                curdate = curdate + dt.timedelta(days=1)

    savedirlist = [x for (y,x) in sorted(zip(urllist,savedirlist))]
    urllist = sorted(urllist)

    print(" ... done")
    print(len(urllist),"potential " + typ)

#    if sorted_mode:
#        print 'aaaa'
#        results = [pool.apply_async(downloader , args=(u,sd)) for u,sd in zip(urllist,savedirlist)]
#    else:
#        print 'bb'
#        _ = pool.map(downloader_wrap,zip(urllist,savedirlist))
#
#
    if not sorted_mode:
        _ = pool.map(downloader_wrap,list(zip(urllist,savedirlist)))
    else:
         results = [pool.apply_async(downloader , args=(u,sd)) for u,sd in zip(urllist,savedirlist)]

    localfiles_lis = []
    
    if not return_also_uncompressed_files:
        for url , savedir in zip(urllist,savedirlist):
            localfile = os.path.join(savedir,os.path.basename(url))
            if os.path.isfile(localfile):
                localfiles_lis.append(localfile)
    else:
    
        for url , savedir in zip(urllist,savedirlist):

            localfile = os.path.join(savedir,os.path.basename(url))        
            
            Pot_compress_files_list = [localfile]
            Pot_compress_files_list.append(localfile.replace(".gz",""))
            Pot_compress_files_list.append(localfile.replace(".Z",""))
            Pot_compress_files_list = list(set(Pot_compress_files_list))
            
            for potential_exisiting_file in Pot_compress_files_list:
                if os.path.isfile(potential_exisiting_file):
                    localfiles_lis.append(potential_exisiting_file)
    
    pool.close()
    return localfiles_lis

def multi_finder_rinex(main_dir,rinex_types=('o','d','d.Z','d.z'),
                       specific_stats = [] ):
    """
    from a main_dir, find all the rinexs in this folder and his subfolder

    (corresponding to the rinex_types)

    and return a list of the found rinexs

    is very similar with geodetik.rinex_lister and  gins_runner.get_rinex_list

    But this one is the most elaborated , must be used in priority !!!
    """
    files_raw_lis , _ = utils.walk_dir(main_dir)

    yylis = [str(e).zfill(2) for e in list(range(80,100)) + list(range(0,dt.datetime.now().year - 2000 + 1))]

    rinex_lis = []


    for f in files_raw_lis:
        for rnxext in rinex_types:
            for yy in yylis:
                if f.endswith(yy + rnxext):
                    rinex_lis.append(f)


    # CASE FOR specific_stats
    if len(specific_stats) > 0:
        rinex_lis2 = []
        for stat in specific_stats:
            for rnx in rinex_lis:
                if stat in os.path.basename(rnx):
                    rinex_lis2.append(rnx)
        rinex_lis = rinex_lis2

    print('INFO : ' , len(rinex_lis) , 'RINEXs found')

    return rinex_lis


def multi_archiver_rinex(rinex_lis,parent_archive_dir,archtype='stat',
                         move=True,force_mv_or_cp=True):
    """
    from rinex_lis, a list of rinex (generated by the function
    multi_finder_rinex)

    move (if move=True) of copy (if move=False) those rinexs in the
    parent_archive_dir according to the archtype,
    string describing how the archive directory is structured, e.g :
            stat

            stat/year

            stat/year/doy

            year/doy

            year/stat

            week/dow/stat

            ... etc ...
    """

    mv_cnt   = 0
    skip_cnt = 0
    print('INFO : RINEXs as input :' , len(rinex_lis))

    if move:
        mv_fct = utils.move
    else:
        mv_fct = utils.copy2

    for rnx in rinex_lis:
        date = conv.rinexname2dt(rnx)
        rnxname = os.path.basename(rnx)
        stat = rnxname[0:4]
        savedir = effective_save_dir(parent_archive_dir, stat, date, archtype)
        utils.create_dir(savedir)

        if not force_mv_or_cp and os.path.isfile(os.path.join(savedir,rnxname)):
            skip_cnt += 1
            continue
        else:
            mv_fct(rnx,savedir)
            mv_cnt += 1

    print('INFO : RINEXs skiped   :' , skip_cnt, '(because already existing)')
    print('INFO : RINEXs moved    :' , mv_cnt)

    return None


def find_IGS_products_files(parent_dir,File_type,ACs,date_start,date_end=None,
                            recursive_search=True,severe=True,
                            compressed="incl",
                            regex_old_naming = True,
                            regex_new_naming = True,
                            regex_igs_tfcc_naming = True,
                            add_weekly_file = False):
    """
    Find all product files in a parent folder which correspond to file type(s),
    AC(s) and date(s)

    Parameters
    ----------
    parent_dir : str or list of str
        The parent directory (i.e. the archive) where files are stored
        can be a string (path of the archive) or a list of file paths
        (given by the function utils.find_recursive) in order to gain time


    File_type : str or list of str
        File type(s) researched (sp3, erp, clk ...)
        can be a list of string for several file type paths
        or a string like 'sp3' if only one file type is researched

    ACs : 'all' or str or list of str
        AC(s) researched
        can be a list of string for several ACs
        or a string like 'gfz' if only one AC is researched
        if 'all', search for all the ACs

    date_start : dt.datetime or 2-tuple list of int
        begining of the time period researched
        can be a datetime
        or a 2-tuple (wwww,d) e.g. (1990,0)

    date_end : None or dt.datetime or 2-tuple list of int
        end of the time period researched
        can be a datetime
        or a 2-tuple (wwww,d) e.g. (1990,0)
        if None, then only date_start is researched
    
    severe : bool
        If True, raises an exception if something goes wrong

    compressed : str
        How the compressed files are handled
        "incl": include the compressed files
        "only": only consider the compressed files
        "excl": exclude the compressed files
        
    regex_old_naming : bool
        Handle old naming format

    regex_new_naming : bool
        Handle new naming format        
    
    regex_igs_tfcc_naming : bool
        Handle TFCC specific format (for SINEX files)
        
    add_weekly_file : bool
        Also handle the weekly file (day 7)
        Implemented only for the  old naming format (for the moment)

    Returns
    -------
    Files_select_cumul_list : list
        List of files found

    """
    
    ###### Prelim Checks   ##############
    if not os.path.exists(parent_dir):
        print("ERR : parent directory doesn't exist")
        print(parent_dir)
        if severe:
            raise Exception

    ###### Time management ##############
    ###### work internally with datetime
    ### date_start
    if type(date_start) is dt.datetime:
        date_start_ok = date_start
    else:
        date_start_ok = conv.gpstime2dt(*date_start)
    ### date_end
    if not date_end:
        date_end_ok = date_start_ok
    elif type(date_end) is dt.datetime:
        date_end_ok = date_end
    else:
        date_end_ok = conv.gpstime2dt(*date_end)
    ### generate time period with a while loop
    Dates_list = [date_start_ok]
    while Dates_list[-1] < date_end_ok:
        Dates_list.append(Dates_list[-1]  + dt.timedelta(days=1))
        
    ### manage weekly file 
    Dates_wwwwd_list   = [utils.join_improved("",*conv.dt2gpstime(d,outputtype=str)) for d in Dates_list]
    Dates_yyyyddd_list = [utils.join_improved("",*reversed(conv.dt2doy_year(d))) for d in Dates_list]
    
    

    ###### File type / ACs management ##############

    if not utils.is_iterable(File_type):
        File_type = [File_type]
    if not utils.is_iterable(ACs):
        ACs = [ACs]

    ###### General file search management ##############
    if utils.is_iterable(parent_dir):
        FILE_LIST = parent_dir
    elif recursive_search:
        # All the files are listed first
        FILE_LIST = utils.find_recursive(parent_dir,".*",case_sensitive=False)
    else:
        FILE_LIST = glob.glob(parent_dir + "/*")


    ###### Regex Definition ##############
    
    join_regex_and = lambda  L : "(" +  "|".join(L) + ")"
    
    Re_patt_big_stk = []
    
    
    ### compression handeling
    if compressed == "excl":
        re_patt_comp = "$"
    elif compressed == "incl":
        re_patt_comp = "(\.Z|\.gz|)$"
    elif compressed == "only":
        re_patt_comp = "(\.Z|\.gz)$"
    else:
        print("ERR : check 'compressed' keyword (excl,incl, or only)")
        raise Exception
        
    if regex_old_naming: 
        if ACs[0] == "all":
            re_patt_ac = "\w{3}"
        else:
            re_patt_ac = join_regex_and([ac.lower() for ac in ACs])
                 
        if add_weekly_file:
            Dates_wwwwd_list_4old = [ e[:-1] + "(" + e[-1] + "|" + "7)"  for e in Dates_wwwwd_list]
        else:
            Dates_wwwwd_list_4old = Dates_wwwwd_list
        
        re_patt_date   = join_regex_and(Dates_wwwwd_list_4old)
        re_patt_filtyp = join_regex_and(File_type)
        re_patt_big_old_naming = re_patt_ac + re_patt_date + "\." + re_patt_filtyp + re_patt_comp
        Re_patt_big_stk.append(re_patt_big_old_naming)
        
    if regex_new_naming: ### search for new name convention
        if ACs[0] == "all":
            re_patt_ac = "\W{3}"        
        else:
            re_patt_ac = join_regex_and([ac.upper() for ac in ACs])
        re_patt_date   = join_regex_and(["_"+e for e in Dates_yyyyddd_list]) #add _ because it can raise a conflit with the old format
        re_patt_filtyp = join_regex_and([fil.upper() for fil in File_type])
        re_patt_big_new_naming = ".*".join((re_patt_ac,re_patt_date,re_patt_filtyp + re_patt_comp))
        Re_patt_big_stk.append(re_patt_big_new_naming)
        
    if regex_igs_tfcc_naming:
        Dates_yy_list = list(set([str(conv.gpstime2dt(int(e[0:4]),int(e[4])).year)[2:] for e in Dates_wwwwd_list]))
        Dates_wwww_list = list(set([e[:-1] for e in Dates_wwwwd_list]))
        #Dates_wwww_dot_list = [e + "\." for e in Dates_wwww_list]
        re_patt_year = join_regex_and(Dates_yy_list)
        # 2x re_patt_date : because .sum doesn't the day
        re_patt_date = join_regex_and(Dates_wwwwd_list + Dates_wwww_list)
        re_patt_filtyp = "\." +  join_regex_and(File_type)

        re_patt_big_igs_tfcc_naming = "igs" + re_patt_year + "P" + re_patt_date + ".*" + re_patt_filtyp + re_patt_comp
        Re_patt_big_stk.append(re_patt_big_igs_tfcc_naming)

    re_patt_big = join_regex_and(Re_patt_big_stk)
        
    print("INFO : REGEX researched :")    
    print(re_patt_big)
    
    ###### Specific file search management ##############
    Files_select_list = []
    for fil in FILE_LIST:
        if re.search(re_patt_big,os.path.basename(fil),re.IGNORECASE):
            Files_select_list.append(fil)
            
    if len(Files_select_list) == 0:
        print("ERR : no products found")
        if severe:
            raise Exception

    print("INFO : number of files found :",len(Files_select_list))    
    print(re_patt_big)
    
    return Files_select_list

def FTP_downloader(ftp_obj,filename,localdir):   
    localpath = os.path.join(localdir,filename)
    
    if not utils.empty_file_check(localpath):
        print("INFO:",filename,"already exists ;)")
        bool_dl = True
    else:
        try:
            localfile = open(localpath, 'wb')
            ftp_obj.retrbinary('RETR ' + filename, localfile.write, 1024)
            localfile.close()
            bool_dl = True
            print("INFO:",filename,"downloaded :)")
    
        except:
            print("WARN:",localpath,"download failed :(")
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
    

def multi_downloader_orbs_clks_2(archive_dir,startdate,enddate,
                            AC_names = ("wum","cod"),
                            prod_types = ("sp3","clk"),
                            remove_patterns=("ULA",),
                            archtype ='week',
                            new_name_conv = True,
                            parallel_download=4,
                            archive_center='ign',
                            mgex=True,repro=0,sorted_mode=False,
                            return_also_uncompressed_files=True,
                            ftp_download=False,
                            dow_manu=False):
    """
    dow_manu = False, no dow manu, consider the converted dow from the time span, regular case
    dow_manu = None, no dow in the REGEX, the crawler will search only for the week
    dow_manu = 0 or 7: the dow in question    
    """
    
    if mgex:
        mgex_str = "mgex/"
    else:
        mgex_str = ""
        
    if not utils.is_iterable(remove_patterns):
        remove_patterns = [remove_patterns]
    
    if archive_center == "cddis":
        arch_center_main    = 'cddis.gsfc.nasa.gov'
        arch_center_basedir = '/pub/gps/products/' + mgex_str
        
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


    print("INFO : data center used :",archive_center)

    Dates_list = conv.dt_range(startdate,enddate)

    Localfiles_lis = []
    wwww_dir_previous = None
    pool = mp.Pool(processes=parallel_download) 
   

    ## create the FTP object
    ftp = FTP(arch_center_main)
    ftp.login()
    ## create a list of FTP object for multiple downloads (unstable...)
    if ftp_download and parallel_download > 1:
        Ftp_obj_list = [FTP(arch_center_main) for i in range(parallel_download)]
        [f.login() for f in Ftp_obj_list]    
        # define the main obj for crawling
        ftp = Ftp_obj_list[0]
    
    ### check if the pattern of the wished products are in the listed daily files
    for patt_tup in list(itertools.product(Dates_list,AC_names,prod_types)):
        dt_cur , ac_cur , prod_cur = patt_tup
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
               
        print("INFO : ","Search products for day",wwww,dow,"AC/prod",ac_cur,prod_cur)
        wwww_dir = os.path.join(arch_center_basedir,str(wwww))
        print("       Move to:",wwww_dir)
        if wwww_dir_previous != wwww_dir:
            try:
                ftp.cwd(wwww_dir)
            except:
                print("WARN:",wwww_dir,"do not exists, skiping...")
            Files_listed_in_FTP = ftp.nlst()
            wwww_dir_previous = wwww_dir
            if len(Files_listed_in_FTP) == 0:
                print("WARN: no files found in directory",wwww_dir)
                
                
        Files_remote_date_list = []

        pattern_old_nam = ac_cur+".*"+str(wwww)+str(dow)+".*"+prod_cur+"\..*"
        Files = [f for f in Files_listed_in_FTP if re.search(pattern_old_nam,f)]
        
        pattern_new_nam = ""
        
        Files_new_nam = []
        if new_name_conv: ### search for new name convention
            
            if dow is None:
                print("ERR: dow == None and search for new name convention, Error ...")
                raise  Exception()
                
            ac_newnam   = ac_cur.upper()

            doy_newnam  = "".join(reversed(conv.dt2doy_year(conv.gpstime2dt(wwww,dow))))
            prod_newnam = prod_cur.upper()
            
            pattern_new_nam = utils.join_improved(".*",ac_newnam,doy_newnam,prod_newnam)
            pattern_new_nam = ".*" + pattern_new_nam + "\..*"

            Files_new_nam   = [f for f in Files_listed_in_FTP if re.search(pattern_new_nam,f)]
        
        
        print("      ","Regex :",pattern_old_nam,pattern_new_nam)        
        Files = Files + Files_new_nam
        
        if len(Files) == 0:
            print("WARN : ","no product found for",*patt_tup)
        
        Files_remote_date_list = Files_remote_date_list + Files
            
        ### exclude some pattern
        for negpatt in remove_patterns:
            Files_remote_date_list = [e for e in Files_remote_date_list if not re.search(negpatt,e)]

        archive_dir_specif = effective_save_dir_orbit(archive_dir,
                                                      ac_cur,
                                                      dt_cur,
                                                      archtype)
            
        utils.create_dir(archive_dir_specif)
        
        ### Generation of the Download fct inputs
        Files_remote_date_chunck = utils.chunkIt(Files_remote_date_list,
                                                 parallel_download)
        Downld_tuples_list = []
        Potential_localfiles_list = []

        if ftp_download: ### FTP download is not recommended
            for ftpobj , Chunk in zip(Ftp_obj_list,Files_remote_date_chunck):
                for filchunk in Chunk:
                        Potential_localfiles_list.append(os.path.join(archive_dir_specif,filchunk))
                        if parallel_download == 1:
                            Downld_tuples_list.append((ftpobj,filchunk,archive_dir_specif))
                        else:
                            Downld_tuples_list.append((arch_center_main,wwww_dir,
                                                       filchunk,archive_dir_specif))
        else: ### HTML download, recommended
            Downld_tuples_list = itertools.product(["/".join(('ftp://' + arch_center_main,wwww_dir,f)) for f in Files_remote_date_list],[archive_dir_specif])
            [Potential_localfiles_list.append(os.path.join(archive_dir_specif,f)) for f in Files_remote_date_list]

            
        ### Actual Download, FTP download is not recommended
        if ftp_download and parallel_download == 1:
            for tup in Downld_tuples_list:
                FTP_downloader(*tup)
        elif ftp_download and parallel_download > 1:
            _ = pool.map_async(FTP_downloader_wo_objects,Downld_tuples_list)
        elif not ftp_download and parallel_download == 1:
            for tup in Downld_tuples_list:
                downloader_wrap(tup)
        elif not ftp_download and parallel_download > 1:
            _ = pool.map(downloader_wrap,Downld_tuples_list)
        ## Those 2 other methods are unstables
        ##  _ = pool.map_async(downloader_wrap,Downld_tuples_list)
        ##  _ = [pool.apply(downloader_wrap,(tup,)) for tup in Downld_tuples_list]
    

    ### Independent files existence check
    if not return_also_uncompressed_files:
        for localfile in Potential_localfiles_list:
            if os.path.isfile(localfile):
                Localfiles_lis.append(localfile)
    else:
        for localfile in Potential_localfiles_list:
            Pot_compress_name_list = [localfile]
            Pot_compress_name_list.append(localfile.replace(".gz",""))
            Pot_compress_name_list.append(localfile.replace(".Z",""))
            Pot_compress_name_list = list(set(Pot_compress_name_list))
            
            for pot_compress_name in Pot_compress_name_list:
                if os.path.isfile(pot_compress_name):
                    Localfiles_lis.append(pot_compress_name)
    
    pool.close()
    return Localfiles_lis
