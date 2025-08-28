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
#### Import the logger
import logging
import multiprocessing as mp
import os

import geodezyx.operational.download_utils as dlutils
#### geodeZYX modules
from geodezyx import conv, utils

log = logging.getLogger('geodezyx')


##########  END IMPORT  ##########


# _____________ _   _ ________   __  _____                      _                 _
# |  __ \|_   _| \ | |  ____\ \ / / |  __ \                    | |               | |
# | |__) | | | |  \| | |__   \ v /  | |  | | _____      ___ __ | | ___   __ _  __| | ___ _ __
# |  _  /  | | | . ` |  __|   > <   | |  | |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |/ _ \ '__|
# | | \ \ _| |_| |\  | |____ / . \  | |__| | (_) \ v  v /| | | | | (_) | (_| | (_| |  __/ |
# |_|  \_\_____|_| \_|______/_/ \_\ |_____/ \___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|\___|_|


############################################################################
######## RINEX DOWNLOADER
############################################################################


def igs_sopac_server_legacy(stat, date):
    # plante si trop de requete
    urlserver = "ftp://garner.ucsd.edu/pub/rinex/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def igs_cddis_server_legacy(stat, date, user="", passwd=""):
    # a privilegier
    urlserver = "ftp://gdc.cddis.eosdis.nasa.gov/gps/data/daily/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(
        urlserver, str(date.year), conv.dt2doy(date), date.strftime("%y") + "d", rnxname
    )
    if not passwd:
        passwd = "sakic@ipgp.fr"
    return url, user, passwd


def igs_cddis_nav_server_legacy(stat, date):
    # a privilegier
    urlserver = "ftp://cddis.gsfc.nasa.gov/gps/data/daily/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date, "n.Z")
    url = os.path.join(
        urlserver, str(date.year), conv.dt2doy(date), date.strftime("%y") + "n", rnxname
    )
    return url


def nav_bkg_server(stat, date):
    urlserver = "ftp://igs-ftp.bkg.bund.de/IGS/BRDC/"
    # ftp://igs-ftp.bkg.bund.de/IGS/BRDC/2024/082/BRDC00WRD_S_20240820000_01D_MN.rnx.gz
    rnxname = "BRDC00WRD_S_" + conv.dt2str(date, "%Y%j") + "0000_01D_MN.rnx.gz"
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def nav_rob_server_legacy(stat, date):
    urlserver = "ftp://epncb.oma.be/pub/obs/BRDC/"
    # ftp://epncb.oma.be/pub/obs/BRDC/2018/BRDC00GOP_R_20180010000_01D_MN.rnx.gz
    rnxname = "BRDC00GOP_R_" + conv.dt2str(date, "%Y%j") + "0000_01D_MN.rnx.gz"
    url = os.path.join(urlserver, str(date.year), rnxname)
    return url


def rgp_ign_smn_server_legacy(stat, date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), "data_30", rnxname)
    return url


def rgp_ign_mlv_server(stat, date):
    urlserver = "ftp://rgpdata.ensg.eu/pub/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), "data_30", rnxname)
    return url


def rgp_ign_smn_1_hz_server_legacy(stat, date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"

    urls = []

    for h in range(24):
        date_session = date
        date_session = date_session.replace(hour=h)

        log.info("%s session %s", date_session, h)
        rnxname = conv.statname_dt2rinexname(
            stat.lower(), date_session, session_a_instead_of_daily_session=1
        )
        url = os.path.join(
            urlserver, str(date.year), conv.dt2doy(date), "data_1", rnxname
        )

        urls.append(url)

    return urls


def unavco_server_legacy(stat, date):
    urlserver = "ftp://data-out.unavco.org/pub/rinex"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, "obs", str(date.year), conv.dt2doy(date), rnxname)
    return url


def renag_server(stat, date):
    urlserver = "ftp://renag.unice.fr/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def uwiseismic_server(stat, date, user="", passwd=""):
    urlserver = "ftp://www2.uwiseismic.com/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, "rinex", str(date.year), conv.dt2doy(date), rnxname)
    return url, user, passwd


def orpheon_server_legacy(stat, date, user="", passwd=""):
    urlserver = "ftp://renag.unice.fr/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url, user, passwd


def ovsg_server_legacy(stat, date, user="", passwd=""):
    if dt.datetime(2009, 1, 1) <= date <= dt.datetime(2014, 2, 10):
        urlserver = (
            "http://webobs.ovsg.univ-ag.fr/rawdata/GPS-GPSDATA.backtemp_20140210/"
        )
    else:
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS/GPSDATA/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), "rinex", rnxname)
    return url, user, passwd


def geoaus_server_legacy(stat, date):
    """Geosciences Australia
    ex : ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/2010/10063/"""
    urlserver = "ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(
        urlserver, str(date.year), date.strftime("%y") + conv.dt2doy(date), rnxname
    )
    return url


def sonel_server_legacy(stat, date):
    """ex : ftp://ftp.sonel.org/gps/data/2015/001/"""
    urlserver = "ftp://ftp.sonel.org/gps/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def ens_fr_legacy(stat, date):
    urlserver = "ftp://gnss.ens.fr/pub/public/crl/GPS/rinex/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def effective_save_dir(parent_archive_dir, stat, date, archtype="stat"):
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
        ... etc ..."""
    if archtype == "/":
        return parent_archive_dir

    if len(archtype) > 0 and archtype[0] == "/":
        log.warn(
            "The archive type description starts with a /, remove it to avoid an error"
        )

    out_save_dir = parent_archive_dir
    fff = archtype.split("/")
    year = str(date.year)
    doy = conv.dt2doy(date)
    _, _ = year, doy  ## simply to remove the unused linter warning...
    week, dow = conv.dt2gpstime(date)
    for f in fff:
        out_save_dir = os.path.join(out_save_dir, eval(f))
    return out_save_dir


def multi_downloader_rinex(**kwargs):
    log.warn(
        "multi_downloader_rinex is a legacy alias for the newly renamed function download_gnss_rinex"
    )
    return download_gnss_rinex_legacy(**kwargs)


def download_gnss_rinex_legacy(
    statdico,
    archive_dir,
    startdate,
    enddate,
    archtype="stat",
    parallel_download=4,
    sorted_mode=False,
    user="",
    passwd="",
    filter_ftp_crawler=True,
    path_ftp_crawled_files_save=None,
    path_ftp_crawled_files_load=None,
    quiet_mode=False,
    final_archive_for_sup_check=None,
    force=False,
):
    """
    Download GNSS RINEX files from various archive centers.

    Parameters
    ----------
    statdico : dict
        Dictionary associating archive centers to list of stations.
        Format: statdico['archive_center'] = ['STA1', 'STA2', ...]
        
        Supported archive centers:
        - 'igs_cddis' or 'igs' : CDDIS center
        - 'igs_sopac' : SOPAC/UCSD/SIO center (not very reliable)
        - 'rgp' : IGN's RGP St Mandé center
        - 'rgp_mlv' : IGN's RGP Marne la Vallée center  
        - 'rgp_1Hz' : IGN's RGP, downloads all 24 hourly RINEX for the day
        - 'renag' : RENAG network
        - 'ovsg' : OVSG network
        - 'unavco' : UNAVCO center
        - 'sonel' : SONEL network
        - 'geoaus' : Geosciences Australia
        - 'ens_fr' : ENS France
        - 'nav' or 'brdc' : Navigation files from ROB server using GOP files
        - 'nav_rt' or 'brdc_rt' : Real-time navigation files from BKG server
    archive_dir : str
        Parent archive directory where files will be stored.
    startdate : datetime
        Start date of the wished period.
    enddate : datetime  
        End date of the wished period.
    archtype : str, optional
        String describing how the archive directory is structured.
        Examples: 'stat', 'stat/year', 'stat/year/doy', 'year/doy', 
        'year/stat', 'week/dow/stat', etc. Default is "stat".
    parallel_download : int, optional
        Number of parallel downloads. Default is 4.
    sorted_mode : bool, optional
        If False, uses the map multiprocess function so the download order 
        will be scrambled. If True, uses the apply multiprocess function so 
        the download order will be in chronological order. The scrambled 
        (False) is better because it doesn't create zombie processes. 
        Default is False.
    user : str, optional
        Username for locked servers. Default is "".
    passwd : str, optional
        Password for locked servers. Default is "".
    filter_ftp_crawler : bool, optional
        Use an improved FTP crawler to find which files actually exist
        to accelerate the download. If path_ftp_crawled_files_load is given,
        use this loaded list. Default is True.
    path_ftp_crawled_files_save : str, optional
        Path (directory + filename) to save a pickle containing the list of 
        existing RINEXs found on the server by the FTP crawler. Allows using 
        this list directly if facing a timeout during the download part.
        Default is None.
    path_ftp_crawled_files_load : str, optional
        Load and use the list of existing RINEXs found on the FTP server,
        generated by a previous run of the FTP crawler. Overrides an internal 
        call of crawl_ftp_files. Default is None.
    quiet_mode : bool, optional
        List the available RINEXs without downloading them. Useful only if 
        path_ftp_crawled_files_save is given. Default is False.
    final_archive_for_sup_check : str, optional
        The final archive path or a file containing the archived RINEXs in
        their final destination. Useful if the final archive is different 
        from archive_dir. Default is None.
    force : bool, optional
        Force the download even if the file already exists locally. 
        Default is False.

    Returns
    -------
    url_list : list of str
        List of URLs.
    savedir_list : list of str
        List of downloaded products paths.

    Examples
    --------
    >>> statdic = dict()
    >>> statdic['igs_cddis'] = ['ZIMM']
    >>> output_dir = '/home/USER/test_dl_rnx'
    >>> startdate = dt.datetime(2000, 1, 1)
    >>> enddate = dt.datetime(2000, 1, 31)
    >>> geodezyx.operational.download_gnss_rinex_legacy(statdic,
    """
    
    curdate = startdate

    pool = mp.Pool(processes=parallel_download)

    urllist = []
    savedirlist = []

    log.info("generating the list of potential RINEXs ...")

    while curdate <= enddate:
        for netwk, statlis in list(statdico.items()):

            if not utils.is_iterable(statlis):
                log.warning(
                    "%s given in the 'statdico' should be a list and it is not.",
                    statlis,
                )

            for stat in statlis:
                stat = stat.lower()
                mode1hz = False
                secure_ftp = False

                if netwk in ("igs_cddis", "igs"):
                    secure_ftp = True
                    url = igs_cddis_server_legacy(stat, curdate, user, passwd)
                elif netwk == "igs_sopac":
                    url = igs_sopac_server_legacy(stat, curdate)
                elif netwk == "rgp":
                    url = rgp_ign_smn_server_legacy(stat, curdate)
                elif netwk == "rgp_mlv":
                    url = rgp_ign_mlv_server(stat, curdate)
                elif netwk == "rgp_1Hz":
                    urls = rgp_ign_smn_1_hz_server_legacy(stat, curdate)
                    mode1hz = True
                elif netwk == "renag":
                    url = renag_server(stat, curdate)
                elif netwk == "orpheon":
                    url = orpheon_server_legacy(stat, curdate, user, passwd)
                elif netwk == "uwiseismic":
                    url = uwiseismic_server(stat, curdate, user, passwd)
                elif netwk == "ovsg":
                    url = ovsg_server_legacy(stat, curdate, user, passwd)
                elif netwk == "unavco":
                    url = unavco_server_legacy(stat, curdate)
                elif netwk == "sonel":
                    url = sonel_server_legacy(stat, curdate)
                elif netwk == "geoaus":
                    url = geoaus_server_legacy(stat, curdate)
                elif netwk in ("nav", "brdc"):
                    url = nav_rob_server_legacy(stat, curdate)
                elif netwk in ("nav_rt", "brdc_rt"):
                    url = nav_bkg_server(stat, curdate)
                elif netwk == "ens_fr":
                    url = ens_fr_legacy(stat, curdate)
                else:
                    log.warning("unkwn server dic in the dico, skip ...")
                    continue

                savedir = effective_save_dir(archive_dir, stat, curdate, archtype)

                if not mode1hz:
                    urllist.append(url)
                    savedirlist.append(savedir)
                else:
                    urllist = urllist + urls
                    savedirlist = savedirlist + [savedir] * len(urls)

        curdate = curdate + dt.timedelta(days=1)

    # savedirlist = [x for (y,x) in sorted(zip(urllist,savedirlist))]
    # urllist     = sorted(urllist)

    try:
        urllist, savedirlist = utils.sort_binom_list(urllist, savedirlist)
    except TypeError as err:
        log.error("unable to sort the URL and the save directory")
        log.error(
            "TIP: you maybe asked for servers with & without password in the same statdico"
        )
        raise err

    log.info(" ... done")
    log.info(str(len(urllist)) + " potential RINEXs")

    ### Use of the advanced FTP Crawler
    if filter_ftp_crawler:
        if path_ftp_crawled_files_load:  ## if the previous files are loaded
            urllist, savedirlist = utils.pickle_loader(path_ftp_crawled_files_load)
        else:  ## regular case
            urllist, savedirlist = dlutils.ftp_files_crawler_legacy(
                urllist, savedirlist, secure_ftp=secure_ftp
            )
            if path_ftp_crawled_files_save:
                savetup = (urllist, savedirlist)
                utils.pickle_saver(savetup, full_path=path_ftp_crawled_files_save)

    ##### Check if the rinex file already exists in the final archive
    if final_archive_for_sup_check and not force:
        Files_final_arch = utils.find_recursive(final_archive_for_sup_check, "*")
        Files_final_arch_basename = [os.path.basename(e) for e in Files_final_arch]

        urllist_new, savedirlist_new = [], []

        for u, sd in zip(urllist, savedirlist):
            if not os.path.basename(u) in Files_final_arch_basename:
                urllist_new.append(u)
                savedirlist_new.append(sd)

        urllist, savedirlist = urllist_new, savedirlist_new

    if not quiet_mode:
        if not secure_ftp:  ### all the servers except CDDIS
            if sorted_mode:
                _ = [
                    pool.apply_async(dlutils.downloader, args=(u, sd, force))
                    for u, sd in zip(urllist, savedirlist)
                ]
            else:
                forcelis = [force] * len(urllist)
                _ = pool.map(
                    dlutils.downloader_wrap, list(zip(urllist, savedirlist, forcelis))
                )
        else:  ## secure FTP i.e. CDDIS
            ftp_obj, _ = dlutils.ftp_objt_create(
                secure_ftp_inp=secure_ftp,
                host=url[0].split("/")[2],
                user=url[1],
                passwd=url[2],
            )

            for iurl, isavedir in zip(urllist, savedirlist):
                localpath, bool_dl = dlutils.ftp_downld_mono(ftp_obj, iurl[0], isavedir)

    localfiles_lis = []
    skiped_url = 0
    for url, savedir in zip(urllist, savedirlist):
        try:
            if type(url) is tuple:
                url0 = url[0]
            else:
                url0 = url
            localfile = os.path.join(savedir, os.path.basename(url0))
            if os.path.isfile(localfile):
                localfiles_lis.append(localfile)
        except Exception as e:
            # because of a weird error
            # i = p.rfind('/') + 1
            # AttributeError: 'tuple' object has no attribute 'rfind'
            skiped_url += 1
            log.warning(e)
            continue

    pool.close()
    if skiped_url > 0:
        log.debug(str(skiped_url) + " returned url skipped because of an Exception")
    return localfiles_lis, savedirlist


#    return zip(urllist,savedirlist)


def rnx_long2short_name(longname_filepath_in):
    """
    MUST BE IMPROVED
    """

    longname_basename = os.path.basename(longname_filepath_in)
    longname_dirname = os.path.dirname(longname_filepath_in)

    Longname_basename_splitted = longname_basename.split("_")

    datepart_str = Longname_basename_splitted[2]
    yyyy = datepart_str[:4]
    ddd = datepart_str[4:7]

    shortname_basename = longname_basename[:4].lower() + ddd + "0." + yyyy[2:] + "o"

    return os.path.join(longname_dirname, shortname_basename)


def multi_archiver_rinex(
    rinex_lis, parent_archive_dir, archtype="stat", move=True, force_mv_or_cp=False
):
    """
    Move or copy RINEX files to an archive directory structure.
    
    Takes a list of RINEX files and organizes them in the parent_archive_dir 
    according to the specified archive type structure.

    Parameters
    ----------
    rinex_lis : list of str
        List of RINEX file paths (typically generated by multi_finder_rinex).
    parent_archive_dir : str
        Parent directory where files will be archived.
    archtype : str, optional
        String describing how the archive directory is structured.
        Examples: 'stat', 'stat/year', 'stat/year/doy', 'year/doy', 
        'year/stat', 'week/dow/stat', etc. Default is "stat".
    move : bool, optional
        If True, move the files. If False, copy the files. Default is True.
    force_mv_or_cp : bool, optional
        If True, force move/copy even if the file already exists in the 
        destination. Default is False.

    Returns
    -------
    None

    Notes
    -----
    The function logs the number of files processed, skipped, and moved/copied.
    Directory structure is created automatically
    """

    mv_cnt = 0
    skip_cnt = 0
    log.info("RINEXs as input : %s", len(rinex_lis))

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

        if not force_mv_or_cp and os.path.isfile(os.path.join(savedir, rnxname)):
            skip_cnt += 1
            continue
        else:
            mv_fct(rnx, savedir)
            mv_cnt += 1

    log.info("RINEXs skiped :" + str(skip_cnt) + " (because already exist)")
    log.info("RINEXs moved  :" + str(mv_cnt))

    return None
