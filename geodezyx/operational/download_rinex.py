#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 23/04/2024 21:31:45

@author: psakic
"""
import datetime as dt
import glob
import logging
import os
import pathlib
import re

import pandas as pd

import geodezyx.operational.download_utils as dlutils
from geodezyx import conv

log = logging.getLogger(__name__)


def _rnx_obs_rgx(stat, date):
    rnx2rgx = conv.statname_dt2rinexname(stat.lower(), date, rnxtype=".*")
    rnx3rgx = conv.statname_dt2rinexname_long(
        stat,
        date,
        country="...",
        data_source=".",
        file_period="01D",
        data_freq="...",
        data_type=".O",
        format_compression=".*",
    )
    return rnx2rgx, rnx3rgx


def _rnx_nav_rgx(stat, date, sys=".", data_source="."):
    rnx2rgx = conv.statname_dt2rinexname(stat.lower(), date, rnxtype=".*")
    rnx3rgx = conv.statname_dt2rinexname_long(
        stat,
        date,
        country="...",
        data_source=data_source,
        file_period="01D",
        data_freq="",
        data_type=sys + "N",
        format_compression=".*",
    )
    
    return rnx2rgx, rnx3rgx


def _generic_server(stat, date, urlserver):
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(stat, date)
    ### generate urls
    urldir = str(os.path.join(urlserver, str(date.year), conv.dt2doy(date)))
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic, key 2 and 3 are for rinex version
    urldic = dict()
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic


def igs_sopac_server(stat, date):
    # plante si trop de requete
    urlserver = "ftp://garner.ucsd.edu/pub/rinex/"
    urldic = _generic_server(stat, date, urlserver)
    return urldic


def igs_cddis_server(stat, date):
    # plante si trop de requete
    urlserver = "ftp://gdc.cddis.eosdis.nasa.gov/gps/data/daily/"

    # we can not use _generic_server here because of the specific server path structure

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(stat, date)

    ### generate urls
    urldir = os.path.join(
        urlserver, str(date.year), conv.dt2doy(date), date.strftime("%y") + "d"
    )
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic, key 2 and 3 are for rinex version
    urldic = {}
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic


def igs_ign_server(stat, date):
    # plante si trop de requete
    urlserver = "ftp://igs.ign.fr/pub/igs/data/"
    urldic = _generic_server(stat, date, urlserver)
    return urldic


def igs_ign_ensg_server(stat, date):
    # plante si trop de requete
    urlserver = "ftp://igs.ensg.eu/pub/igs/data/"
    urldic = _generic_server(stat, date, urlserver)
    return urldic


def nav_rob_server(stat, date):
    urlserver = "ftp://epncb.oma.be/pub/obs/BRDC/"

    # can not use _generic_server here because of the specific server path structure / file name

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_nav_rgx(stat, date)  ### NAV RNX HERE !!!

    ### generate urls
    urldir = os.path.join(urlserver, str(date.year))  ## NO DOY FOR THIS ONE !!!
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic, key 2 and 3 are for rinex version
    urldic = {}
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic


def sonel_server(stat, date):
    urlserver = "ftp://ftp.sonel.org/gps/data/"
    urldic = _generic_server(stat, date, urlserver)
    return urldic


def euref_server(stat, date):
    urlserver = "ftp://epncb.oma.be/pub/obs/"

    # can not use _generic_server here because of upper case RINEX 2 names
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(stat, date)
    rnx2rgx = rnx2rgx.upper()

    ### generate urls
    urldir = str(os.path.join(urlserver, str(date.year), conv.dt2doy(date)))
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic, key 2 and 3 are for rinex version
    urldic = dict()
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic

def nav_bkg_server(stat, date):
    urlserver = "ftp://igs-ftp.bkg.bund.de/IGS/BRDC/"
    # ftp://igs-ftp.bkg.bund.de/IGS/BRDC/2024/082/BRDC00WRD_S_20240820000_01D_MN.rnx.gz
    
    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_nav_rgx(stat, date, sys="M", data_source="S") ### NAV RNX HERE !!!

    ### generate urls
    urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date))
    rnx3url = os.path.join(urldir, rnx3rgx)
    
    urldic = dict()
    urldic[3] = rnx3url

    return urldic


############ not adapted yet after april 2024 mods
def igs_cddis_nav_server_legacy(stat, date):
    # table_proto privilegier
    urlserver = "ftp://cddis.gsfc.nasa.gov/gps/data/daily/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date, "n.Z")
    url = os.path.join(
        urlserver, str(date.year), conv.dt2doy(date), date.strftime("%y") + "n", rnxname
    )
    return url



def rgp_ign_smn_server_legacy(stat, date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), "data_30", rnxname)
    return url


def rgp_ign_mlv_server_legacy(stat, date):
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
            stat.lower(), date_session, session_a_instead_of_daily_session=True
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


def renag_server_legacy(stat, date):
    urlserver = "ftp://renag.unice.fr/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def uwiseismic_server_legacy(stat, date, user="", passwd=""):
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


def ens_fr_legacy(stat, date):
    urlserver = "ftp://gnss.ens.fr/pub/public/crl/GPS/rinex/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def _server_select(datacenter, site, curdate):
    mode1hz = False
    secure_ftp = False
    urldic = dict()
    if datacenter in ("igs_cddis", "igs"):
        urldic = igs_cddis_server(site, curdate)
        secure_ftp = True
    elif datacenter == "igs_sopac":
        urldic = igs_sopac_server(site, curdate)
    elif datacenter == "igs_ign":
        urldic = igs_ign_server(site, curdate)
    elif datacenter == "igs_ign_ensg":
        urldic = igs_ign_ensg_server(site, curdate)
    elif datacenter == "sonel":
        urldic = sonel_server(site, curdate)
    elif datacenter == "euref":
        urldic = euref_server(site, curdate)
    elif datacenter in ("nav", "brdc"):
        urldic = nav_rob_server(site, curdate)
    elif datacenter in ('nav_rt', 'brdc_rt'):
        urldic = nav_bkg_server(site, curdate)
    # elif datacenter == 'rgp':
    #     urldic = rgp_ign_smn_server_legacy(site, curdate)
    # elif datacenter == 'rgp_mlv':
    #     urldic = rgp_ign_mlv_server(site, curdate)
    # elif datacenter == 'rgp_1Hz':
    #     urls = rgp_ign_smn_1_hz_server_legacy(site, curdate)
    #     mode1hz = True
    # elif datacenter == 'renag':
    #     urldic = renag_server(site, curdate)
    # elif datacenter == 'orpheon':
    #     urldic = orpheon_server_legacy(site, curdate)
    # elif datacenter == 'uwiseismic':
    #     urldic = uwiseismic_server(site, curdate)
    # elif datacenter == 'ovsg':
    #     urldic = ovsg_server_legacy(site, curdate)
    # elif datacenter == 'unavco':
    #     urldic = unavco_server_legacy(site, curdate)
    # elif datacenter == 'geoaus':
    #     urldic = geoaus_server_legacy(site, curdate)
    # elif datacenter == 'ens_fr':
    #     urldic = ens_fr_legacy(site, curdate)
    else:
        log.warning("unkwn server dic in the dico, skip ...")
        return None, None, None

    return urldic, secure_ftp, mode1hz


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


def _rnx_regex_in_dir(rnx_regex, dir_files_list):
    r = re.compile(rnx_regex)
    l = list(filter(r.search, dir_files_list))
    if len(l) == 0:
        return None
    else:
        return l[0]


def ftp_files_crawler(
    table,
    secure_ftp=False,
    user=None,
    passwd=None,
    path_ftp_crawled_files_save=None,
    path_all_ftp_files_save=None,
    force=False,
):
    """
    filter the table with download_gnss_rinex with an
    optimized FTP crawl
    """

    def _save_crawled_files(table_inp):
        if path_ftp_crawled_files_save:
            table_inp.to_csv(path_ftp_crawled_files_save)
        return None

    def _get_and_save_all_ftp_files(all_ftp_files_stk_inp):
        if all_ftp_files_stk_inp:
            all_ftp_files_out = pd.concat(all_ftp_files_stk_inp)
            all_ftp_files_out.reset_index(drop=True, inplace=True)
        else:
            all_ftp_files_out = pd.Series([], dtype=str)

        if path_all_ftp_files_save:
            all_ftp_files_out.to_csv(path_all_ftp_files_save)

        return all_ftp_files_out

    ### rename the columns
    if user or passwd:
        loginftp = True
    else:
        loginftp = False

    ### Do the correct split for the URLs
    table_use = table.copy()

    #### Initialisation of the 1st variables for the loop
    prev_row_ftpobj = table_use.iloc[0]
    prev_row_cwd = table_use.iloc[0]
    ftp_files_list = []
    all_ftp_files_stk = []
    local_files_lis = []
    count_loop = 0  # restablish the connexion after count_nmax loops (avoid freezing)
    count_nmax = 100
    ftpobj = None

    for irow, row in table_use.iterrows():

        if row["crawled"]:
            continue

        #### do a local check
        if (prev_row_cwd["outdir"] != row["outdir"]) or irow == 0:
            local_files_lis = glob.glob(row["outdir"] + "/*")
        rnxlocal = _rnx_regex_in_dir(row["rnxrgx"], local_files_lis)
        if rnxlocal and not force:
            log.info("%s already exists locally ;)", os.path.basename(rnxlocal))
            table_use.loc[irow, "rnxnam"] = ""
            continue

        count_loop = count_loop + 1  #### must be after local file check
        ####### we recreate a new FTP object if the host URL is not the same
        if (
            row["host"] != prev_row_ftpobj["host"]
            or count_loop > count_nmax
            or count_loop == 1
        ):

            if ftpobj:  ## close previous FTP object
                ftpobj.close()

            ftpobj, _ = dlutils.ftp_objt_create(
                secure_ftp_inp=secure_ftp,
                # chdir='/',
                host=prev_row_ftpobj["host"],
                user=user,
                passwd=passwd,
            )
            prev_row_ftpobj = row
            if count_loop > count_nmax:
                count_loop = 0
                _save_crawled_files(table_use)
                _get_and_save_all_ftp_files(all_ftp_files_stk)

        ####### we recreate a new file list if the date path is not the same
        if (prev_row_cwd["dir"] != row["dir"]) or irow == 0:
            log.info("chdir " + row["dir"])
            ftpobj.cwd("/")

            try:  #### we try to change for the right folder
                ftpobj.cwd(row["dir"])
            except Exception as e:  #### If not possible, then no file in the list
                log.warning("unable to chdir to %s, exception %s", row["dir"], e)
                ftp_files_list = []

            ftp_files_list = dlutils._ftp_dir_list_files(ftpobj)
            ftp_files_series = pd.Series(ftp_files_list)
            ftp_files_series = (
                os.path.join("ftp://", row["host"], row["dir"]) + "/" + ftp_files_series
            )
            all_ftp_files_stk.append(ftp_files_series)

            prev_row_cwd = row

        ####### we check if the files is avaiable
        rnx_return = _rnx_regex_in_dir(row["rnxrgx"], ftp_files_list)

        if rnx_return:
            table_use.loc[irow, "rnxnam"] = rnx_return
            log.info(rnx_return + " found on server :)")
        else:
            table_use.loc[irow, "rnxnam"] = ""
            log.warning(row["rnxrgx"] + " not found on server :(")
        table_use.loc[irow, "crawled"] = True

    rnx_ok = table_use["rnxnam"].str.len().astype(bool)

    join_funal_url = lambda e: os.path.join("ftp://", e["host"], e["dir"], e["rnxnam"])
    table_use.loc[rnx_ok, "url_true"] = table_use.loc[rnx_ok].apply(
        join_funal_url, axis=1
    )

    _save_crawled_files(table_use)
    all_ftp_files = _get_and_save_all_ftp_files(all_ftp_files_stk)

    if ftpobj:
        ftpobj.close()

    return table_use, all_ftp_files


def download_gnss_rinex(
    statdico,
    archive_dir,
    startdate,
    enddate,
    archtype="stat",
    parallel_download=4,
    user="anonymous",
    passwd="anonymous@isp.com",
    path_ftp_crawled_files_save=None,
    path_ftp_crawled_files_load=None,
    skip_crawl=False,
    path_all_ftp_files_save=None,
    quiet_mode=False,
    final_archive_for_sup_check=None,
    force=False,
    get_rnx2=True,
    get_rnx3=True,
):
    """
    Parameters
    ----------
    statdico : dict
        a statdico is a dictionary associating Archives Centers to list of stations

        Exemple:
            >>> statdico['archive center 1'] = ['STA1','STA2','STA3', ...]
            >>> statdico['archive center 2'] = ['STA2','STA1','STA4', ...]

        the supported archive center are (april 2024):
            igs_cddis or igs (CDDIS data center)

            igs_sopac (for the SOPAC/UCSD/SIO data center, but not very reliable)

            igs_ign (IGN's data center, main server at St Mandé)

            igs_ign_ensg (IGN's data center, secondary server at ENSG, Marne-la-Vallée)

            sonel

            euref (EPN data center hosted at ROB)

            nav or brdc as archive center allows to download nav files (using 'BRDC' as station name)
            from the ROB server, using GOP files

            nav_rt or brdc_rt as archive center allows to download *real time* nav files
            from the BKG server

            ***** not reimplemented yet *****
            rgp (IGN's RGP St Mandé center)

            rgp_mlv (IGN's RGP Marne la Vallée center)

            rgp_1Hz (IGN's RGP, all the 24 hourly rinex for the day will be downloaded)

            renag

            ovsg

            unavco

            sonel

            geoaus (Geosciences Australia)

            ens_fr
            ***** not reimplemented yet *****

    archive_dir : str
        the root directory on your local drive were to store the RINEXs

    archtype : str
        string describing how the archive directory is structured, e.g :

            stat

            stat/year

            stat/year/doy

            year/doy

            year/stat

            week/dow/stat

            ... etc ...

    user & passwd : str
        user & password for a secure server

    path_ftp_crawled_files_save : str
        will save at the given path (directory+filename) in a CSV file containing
        the list of the existing RINEXs found on the server by the FTP crawler.
        It allows to use this list directly if one face a timeout during
        the download part.

    path_ftp_crawled_files_load : str
        load and use the list of the existing RINEXs found on the FTP server,
        generated by a previous run of the FTP crawler (called by
        download_gnss_rinex or directly by ftp_files_crawler).
        a new call of ftp_files_crawler can be bypassed with
        `skip_crawl`.

    skip_crawl : bool
        when a list of the existing RINEXs found on the FTP server
        is provided with path_ftp_crawled_files_load, skip
        a new call of ftp_files_crawler.
        Then, we assume the provided list as a complete one,
        ready for the download step

    path_all_ftp_files_save : str
        will save at the given path (directory+filename) in a CSV file
        ALL the remote files found on the FTP server.

    quiet_mode : bool
        List the available RINEXs without downloading them.
        Useful only if path_ftp_crawled_files_save is given

    final_archive_for_sup_check : str
        The final archive path or a file containing the archived RINEXs in
        their final destination.
        useful if the final archive is different from archive_dir
        *** not re-implemented yet ***

    force : bool
        Force the download even if the file already exists locally

    get_rnx2 & get_rnx3 : bool
        limit the search/download to RINEX2 (short names) and/or
        RINEX3 (long names) depending on the boolean given

    Returns
    -------
    url_list : list of str
        list of URLs

    savedir_list : list of str
        list of downloaded products paths

    Minimal exemple
    ---------------
        >>> statdic = dict()
        >>> statdic['igs_cddis'] = ['ZIMM','tlse']
        >>> archive_dir = '/home/USER/test_dl_rnx'
        >>> startdate = dt.datetime(2020,1,1)
        >>> enddate = dt.datetime(2020,1,31)
        >>> geodezyx.operational.download_gnss_rinex(statdic, archive_dir, startdate, enddate)
    """

    date_range = conv.dt_range(startdate, enddate)

    table_proto = []

    for k, v in statdico.items():
        datacenter = k
        site_lis = v

    for date in date_range:
        for site in site_lis:
            urldic, secure_ftp, mode1hz = _server_select(datacenter, site, date)
            if not urldic:
                continue
            outdir = effective_save_dir(archive_dir, site, date, archtype)
            for rnxver, rnxurl in urldic.items():
                if rnxver == 2 and not get_rnx2:
                    continue
                if rnxver == 3 and not get_rnx3:
                    continue
                table_proto.append((date, site, outdir, rnxver, rnxurl))

    table = pd.DataFrame(
        table_proto, columns=["date", "site", "outdir", "ver", "url_theo"]
    )
    table["crawled"] = False
    table["url_true"] = None
    table["rnxnam"] = ""

    urlpathobj = table["url_theo"].apply(pathlib.Path)
    table["rnxrgx"] = urlpathobj.apply(lambda p: p.name)
    table["host"] = urlpathobj.apply(lambda p: p.parts[1])
    table["dir"] = urlpathobj.apply(lambda p: os.path.join(*p.parts[2:-1]))

    if path_ftp_crawled_files_load:
        table = pd.read_csv(path_ftp_crawled_files_load)

    if skip_crawl:
        table_crawl = table
    else:
        table_crawl, files_all = ftp_files_crawler(
            table,
            secure_ftp=secure_ftp,
            user=user,
            passwd=passwd,
            path_ftp_crawled_files_save=path_ftp_crawled_files_save,
            path_all_ftp_files_save=path_all_ftp_files_save,
            force=force,
        )

    #### get only the valid (true) url
    table_dl = table_crawl.loc[table_crawl["url_true"].dropna().index]

    if len(table_dl) == 0:
        log.error(
            "no valid RINEX URL found/selected on the FTP server, check your inputs"
        )

    if not quiet_mode and len(table_dl) > 0:
        dlutils.ftp_download_frontend(
                table_dl["url_true"].values,
                table_dl["outdir"].values,
                parallel_download=parallel_download,
                secure_ftp=secure_ftp,
                user=user,
                passwd=passwd,
                force=force,
                )

    return None
