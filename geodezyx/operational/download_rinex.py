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
# from ftplib import FTP
# import glob
# import itertools
import multiprocessing as mp
import os
# import pandas as pd 
# import re
# import shutil
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



#_____________ _   _ ________   __  _____                      _                 _
#|  __ \|_   _| \ | |  ____\ \ / / |  __ \                    | |               | |
#| |__) | | | |  \| | |__   \ V /  | |  | | _____      ___ __ | | ___   __ _  __| | ___ _ __
#|  _  /  | | | . ` |  __|   > <   | |  | |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |/ _ \ '__|
#| | \ \ _| |_| |\  | |____ / . \  | |__| | (_) \ V  V /| | | | | (_) | (_| | (_| |  __/ |
#|_|  \_\_____|_| \_|______/_/ \_\ |_____/ \___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|\___|_|



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

def rob_nav_server(stat,date):
    # a privilegier
    urlserver = "ftp://epncb.oma.be/pub/obs/BRDC/"
    #ftp://epncb.oma.be/pub/obs/BRDC/2018/BRDC00GOP_R_20180010000_01D_MN.rnx.gz
    rnxname = "BRDC00GOP_R_" + conv.dt2str(date,'%Y%j') + "0000_01D_MN.rnx.gz"
    url = os.path.join(urlserver , str(date.year) , rnxname)
    return url

def rgp_ign_smn_server(stat,date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , conv.dt2doy(date) , 'data_30' , rnxname)
    return url

def rgp_ign_mlv_server(stat,date):
    urlserver = "ftp://rgpdata.ensg.eu/pub/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , conv.dt2doy(date) , 'data_30' , rnxname)
    return url

def rgp_ign_smn_1Hz_server(stat,date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"

    urls = []

    for h in range(24):
        date_session = date
        date_session = date_session.replace(hour=h)

        log.info('%s session %s' , date_session , h)
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

def uwiseismic_server(stat,date,user='',passwd=''):
    urlserver = "ftp://www2.uwiseismic.com/"
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , 'rinex' , str(date.year) , conv.dt2doy(date) , rnxname)
    return url,user,passwd

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
    
    if len(archtype) > 0 and archtype[0] == "/":
        log.warn("The archive type description starts with a /, remove it to avoid an error")

    out_save_dir = parent_archive_dir
    fff = archtype.split('/')
    #year = str(date.year)
    #doy = conv.dt2doy(date)
    week, dow = conv.dt2gpstime(date)
    for f in fff:
        out_save_dir = os.path.join(out_save_dir,eval(f))
    return out_save_dir


def ens_fr(stat,date):    
    urlserver='ftp://gnss.ens.fr/pub/public/crl/GPS/rinex/'
    rnxname = conv.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , conv.dt2doy(date) , rnxname)
    return url

def multi_downloader_rinex(statdico,archive_dir,startdate,enddate,
                           archtype ='stat',parallel_download=4,
                           sorted_mode=False,user='',passwd='',
                           filter_ftp_crawler=True,
                           path_ftp_crawled_files_save=None,
                           path_ftp_crawled_files_load=None,
                           silent_mode=False,
                           final_archive_for_sup_check=None):
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

            rgp_mlv (Marne la Vallée center)

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
        
    filter_ftp_crawler : bool
        use an improved FTP crawler to find which files actually exist
        to accelerate the download. If path_ftp_crawled_files_load is given, 
        use this loaded list.
        
    path_ftp_crawled_files_save : str
        will save at the given path (directory+filname) in a pickle containing 
        the list of the existing RINEXs found on the server by the FTP crawler.
        It allows to use this list directly if one face a timeout during 
        the download part.
        NB for advanced users: the pickle is a tuple (urllist,savedirlist)
    
    path_ftp_crawled_files_load : str
        load and use the list of the existing RINEXs found on the FTP server,
        generated by a previous run of the FTP crawler (called by 
        multi_downloader_rinex or directly by ftp_files_crawler).
        overrides an internal call of ftp_files_crawler.
        
    silent_mode : bool
        List the available RINEXs without downloading them. 
        Useful only if path_ftp_crawled_files_save is given
        
    final_archive_for_sup_check : str
        The final archive path or a file containing the archived RINEXs in
        their final destination. 
        useful if the final archive is different from archive_dir

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


    log.info("generating the list of potential RINEXs ...")
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
                elif netwk == 'rgp_mlv':
                    url = rgp_ign_mlv_server(stat,curdate)
                elif netwk == 'rgp_1Hz':
                    urls = rgp_ign_smn_1Hz_server(stat,curdate)
                    mode1Hz = True
                elif netwk == 'renag':
                    url = renag_server(stat,curdate)
                elif netwk == 'orpheon':
                    url = orpheon_server(stat,curdate,user,passwd)
                elif netwk == 'uwiseismic':
                    url = uwiseismic_server(stat,curdate,user,passwd)
                elif netwk == 'ovsg':
                    url = ovsg_server(stat,curdate,user,passwd)
                elif netwk == 'unavco':
                    url = unavco_server(stat,curdate)
                elif netwk == 'sonel':
                    url = sonel_server(stat,curdate)
                elif netwk == 'geoaus':
                    url = geoaus_server(stat,curdate)
                elif netwk in ('nav' , 'brdc'):
                    url = rob_nav_server(stat,curdate)
                elif netwk == 'ens_fr':
                    url = ens_fr(stat,curdate)
                else:
                    log.warning('unkwn server dic in the dico, skip ...')
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
    
    try:
        urllist,savedirlist = utils.sort_binom_list(urllist,savedirlist)
    except TypeError as err:
        log.error("unable to sort the URL and the save directory") 
        log.error("TIP: you maybe asked for servers with & without password in the same statdico") 
        raise err

    log.info(" ... done")
    log.info(str(len(urllist)) + " potential RINEXs")

    ### Use of the advanced FTP Crawler
    if filter_ftp_crawler:
        if path_ftp_crawled_files_load:
            urllist,savedirlist = utils.pickle_loader(path_ftp_crawled_files_load)
        else:
            urllist,savedirlist = dlutils.ftp_files_crawler(urllist,savedirlist)
            if path_ftp_crawled_files_save:
                savetup = (urllist,savedirlist)
                utils.pickle_saver(savetup,
                                   full_path=path_ftp_crawled_files_save)
                
    ##### Check if the rinex file already exists in the final archive
    if final_archive_for_sup_check:
        Files_final_arch = utils.find_recursive(final_archive_for_sup_check,
                                                "*")
        Files_final_arch_basename = [os.path.basename(e) for e in Files_final_arch]
        
        urllist_new,savedirlist_new = [],[]
        
        for u,sd in zip(urllist,savedirlist):
            if not os.path.basename(u) in Files_final_arch_basename:
                urllist_new.append(u)
                savedirlist_new.append(sd)
        
        urllist,savedirlist = urllist_new,savedirlist_new
                    
    if not silent_mode:
        if sorted_mode:
            _ = [pool.apply_async(dlutils.downloader , args=(u,sd)) for u,sd in zip(urllist,savedirlist)]
        else:
            _ = pool.map(dlutils.downloader_wrap,list(zip(urllist,savedirlist)))


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
    log.debug(str(skiped_url) +  ' returned url skiped because of a weird error, but it is not important...')
    return localfiles_lis , savedirlist

#    return zip(urllist,savedirlist)


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
    log.info('RINEXs as input : %s' , len(rinex_lis))

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

    log.info('RINEXs skiped :' + str(skip_cnt) + ' (because already exist)')
    log.info('RINEXs moved  :' + str(mv_cnt))

    return None


