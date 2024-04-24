#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 23/04/2024 21:31:45

@author: psakic
"""
import re
import glob
import pathlib
import time

import numpy as np
import pandas as pd

from geodezyx import conv, utils
import datetime as dt
import os

import geodezyx.operational.download_utils as dlutils

import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool

import logging

log = logging.getLogger('geodezyx')


def _rnx_rgx(stat, date):
    rnx2rgx = conv.statname_dt2rinexname(stat.lower(), date, rnxtype="*")
    rnx3rgx = conv.statname_dt2rinexname_long(stat,
                                              date,
                                              country="...",
                                              data_source=".",
                                              file_period="01D",
                                              data_freq="...",
                                              data_type=".O",
                                              format_compression='.*')
    return rnx2rgx, rnx3rgx


def igs_sopac_server(stat, date):
    # plante si trop de requete
    urlserver = "ftp://garner.ucsd.edu/pub/rinex/"

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_rgx(stat, date)

    ### generate urls
    urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date))
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic, key 2 and 3 are for rinex version
    urldic = {}
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic


def igs_cddis_server(stat, date):
    # plante si trop de requete
    urlserver = "ftp://gdc.cddis.eosdis.nasa.gov/gps/data/daily/"

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_rgx(stat, date)

    ### generate urls
    urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date), date.strftime('%y') + 'd')
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic, key 2 and 3 are for rinex version
    urldic = {}
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic


############ not adapted yet
def igs_cddis_nav_server(stat, date):
    # table_proto privilegier
    urlserver = "ftp://cddis.gsfc.nasa.gov/gps/data/daily/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date, 'n.Z')
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), date.strftime('%y') + 'n', rnxname)
    return url


def nav_bkg_server(stat, date):
    urlserver = "ftp://igs-ftp.bkg.bund.de/IGS/BRDC/"
    #ftp://igs-ftp.bkg.bund.de/IGS/BRDC/2024/082/BRDC00WRD_S_20240820000_01D_MN.rnx.gz
    rnxname = "BRDC00WRD_S_" + conv.dt2str(date, '%Y%j') + "0000_01D_MN.rnx.gz"
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def nav_rob_server(stat, date):
    urlserver = "ftp://epncb.oma.be/pub/obs/BRDC/"
    #ftp://epncb.oma.be/pub/obs/BRDC/2018/BRDC00GOP_R_20180010000_01D_MN.rnx.gz
    rnxname = "BRDC00GOP_R_" + conv.dt2str(date, '%Y%j') + "0000_01D_MN.rnx.gz"
    url = os.path.join(urlserver, str(date.year), rnxname)
    return url


def rgp_ign_smn_server(stat, date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), 'data_30', rnxname)
    return url


def rgp_ign_mlv_server(stat, date):
    urlserver = "ftp://rgpdata.ensg.eu/pub/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), 'data_30', rnxname)
    return url


def rgp_ign_smn_1Hz_server(stat, date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"

    urls = []

    for h in range(24):
        date_session = date
        date_session = date_session.replace(hour=h)

        log.info('%s session %s', date_session, h)
        rnxname = conv.statname_dt2rinexname(stat.lower(), date_session,
                                             session_a_instead_of_daily_session=1)
        url = os.path.join(urlserver, str(date.year), conv.dt2doy(date),
                           'data_1', rnxname)

        urls.append(url)

    return urls


def unavco_server(stat, date):
    urlserver = 'ftp://data-out.unavco.org/pub/rinex'
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, 'obs', str(date.year), conv.dt2doy(date), rnxname)
    return url


def renag_server(stat, date):
    urlserver = "ftp://renag.unice.fr/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def uwiseismic_server(stat, date, user='', passwd=''):
    urlserver = "ftp://www2.uwiseismic.com/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, 'rinex', str(date.year), conv.dt2doy(date), rnxname)
    return url, user, passwd


def orpheon_server(stat, date, user='', passwd=''):
    urlserver = "ftp://renag.unice.fr/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url, user, passwd


def ovsg_server(stat, date, user='', passwd=''):
    if dt.datetime(2009, 1, 1) <= date <= dt.datetime(2014, 2, 10):
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS-GPSDATA.backtemp_20140210/"
    else:
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS/GPSDATA/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), 'rinex', rnxname)
    return url, user, passwd


def geoaus_server(stat, date):
    """ Geosciences Australia
        ex : ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/2010/10063/ """
    urlserver = "ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), date.strftime('%y') + conv.dt2doy(date), rnxname)
    return url


def sonel_server(stat, date):
    """ex : ftp://ftp.sonel.org/gps/data/2015/001/ """
    urlserver = 'ftp://ftp.sonel.org/gps/data/'
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def ens_fr(stat, date):
    urlserver = 'ftp://gnss.ens.fr/pub/public/crl/GPS/rinex/'
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def effective_save_dir(parent_archive_dir, stat, date, archtype='stat'):
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
    year = str(date.year)
    doy = conv.dt2doy(date)
    _, _ = year, doy  ## simply to remove the unused linter warning...
    week, dow = conv.dt2gpstime(date)
    for f in fff:
        out_save_dir = os.path.join(out_save_dir, eval(f))
    return out_save_dir


def _rnx_regex_in_dir(rnx_regex, dir_files_list):
    r = re.compile(rnx_regex)
    l = list(filter(r.match, dir_files_list))
    if len(l) == 0:
        return None
    else:
        return l[0]


def ftp_files_crawler(table, secure_ftp=False, user=None, passwd=None,
                      path_ftp_crawled_files_save=None):
    """
    filter urllist,savedirlist generated with download_gnss_rinex with an
    optimized FTP crawl

    """

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
    ftp_files_all_stk = []
    count_loop = 0  # restablish the connexion after count_nmax loops (avoid freezing)
    count_nmax = 100

    for irow, row in table_use.iterrows():
        count_loop = count_loop + 1

        #### do a local check
        rnxlocal = _rnx_regex_in_dir(row['rnxrgx'], glob.glob(row['outdir'] + "/*"))
        if rnxlocal:
            log.info("%s already exists locally ;)", rnxlocal)
            table_use.loc[irow, 'rnxnam'] = ''
            continue

        ####### we recreate a new FTP object if the host URL is not the same
        if row['host'] != prev_row_ftpobj['host'] or count_loop > count_nmax or count_loop == 1:
            ftpobj, _ = dlutils.ftp_objt_create(secure_ftp_inp=secure_ftp,
                                                host=prev_row_ftpobj['host'],
                                                user=user,
                                                passwd=passwd)
            prev_row_ftpobj = row
            if count_loop > count_nmax:
                count_loop = 0
                if path_ftp_crawled_files_save:
                    table_use.to_csv(path_ftp_crawled_files_save)

        ####### we recreate a new file list if the date path is not the same
        if (prev_row_cwd["dir"] != row["dir"]) or irow == 0:
            log.info("chdir " + row["dir"])
            ftpobj.cwd("/")

            try:  #### we try to change for the right folder
                ftpobj.cwd(row['dir'])
            except:  #### If not possible, then no file in the list
                ftp_files_list = []

            ftp_files_list = dlutils._ftp_dir_list_files(ftpobj)
            ftp_files_series = pd.Series(ftp_files_list)
            ftp_files_series = os.path.join('ftp://', row['host'], row['dir']) + ftp_files_series
            ftp_files_all_stk.append(ftp_files_series)

            prev_row_cwd = row

        ####### we check if the files is avaiable

        rnx_return = _rnx_regex_in_dir(row["rnxrgx"], ftp_files_list)

        if rnx_return:
            table_use.loc[irow, 'rnxnam'] = rnx_return
            log.info(rnx_return + " found on server :)")
        else:
            table_use.loc[irow, 'rnxnam'] = ''
            log.warning(row["rnxrgx"] + " not found on server :(")

    rnx_ok = table_use['rnxnam'].str.len().astype(bool)

    join_funal_url = lambda e: os.path.join('ftp://', e['host'], e['dir'], e['rnxnam'])
    table_use.loc[rnx_ok, 'url_true'] = table_use.loc[rnx_ok].apply(join_funal_url, axis=1)

    ftp_files_all = pd.concat(ftp_files_all_stk)

    if path_ftp_crawled_files_save:
        table_use.to_csv(path_ftp_crawled_files_save)

    ftpobj.close()

    return table_use , ftp_files_all

def ftp_download_frontend(urllist,
                          savedirlist,
                          parallel_download=1,
                          secure_ftp=False,
                          user="anonymous",
                          passwd="anonymous@anonymous.com",
                          force=True):

    urlpathobj = pd.Series(urllist).apply(pathlib.Path)
    host_use = urlpathobj.apply(lambda p: p.parts[1]).unique()[0]

    ftpobj_main, ftpobj_lis = dlutils.ftp_objt_create(secure_ftp_inp=secure_ftp,
                                                      host=host_use,
                                                      parallel_download=parallel_download,
                                                      user=user,
                                                      passwd=user)

    #for url, savedir in zip(urllist, savedirlist):
    #    localpath, bool_dl = dlutils.ftp_downloader(ftpobj_main, url, savedir)

    ftpobj_mp_lis = ftpobj_lis * int(np.ceil(len(urllist) / parallel_download))

    if len(ftpobj_mp_lis) < len(urllist):
        log.warning("less FTP objects than URL for parallel download, contact the main developper")

    # alternative solution
    # https://stackoverflow.com/questions/62903886/why-is-giving-me-this-error-typeerror-cannot-pickle-io-textiowrapper-object
    pool = ThreadPool(parallel_download)  #mp.Pool(processes=parallel_download)

    _ = pool.map(dlutils.ftp_downloader_wrap, list(zip(ftpobj_mp_lis,
                                                       urllist,
                                                       savedirlist)))
    return


if __name__ == "__main__":
    startdate = dt.datetime(2021, 2, 1)
    enddate = dt.datetime(2021, 2, 5)

    date_range = conv.dt_range(startdate, enddate)

    archive_dir = '/home/psakicki/aaa_FOURBI/test_dl_rnx'
    archtype = 'year/doy/stat'

    user = 'anonymous'
    passwd = 'toto@toto.com'
    secure_ftp = True
    parallel_download = 1

    table_proto = []
    site_lis = ['zimm', 'tlse', 'abmf']
    for date in date_range:
        for site in site_lis:
            urldic = igs_cddis_server(site, date)
            outdir = effective_save_dir(archive_dir, site, date, archtype)
            for rnxver, rnxurl in urldic.items():
                table_proto.append((date, site, outdir, rnxver, rnxurl))

    table = pd.DataFrame(table_proto, columns=['date', 'site', 'outdir', 'ver', 'url_theo'])
    table['url_true'] = None
    table['rnxnam'] = ''

    urlpathobj = table["url_theo"].apply(pathlib.Path)
    table["rnxrgx"] = urlpathobj.apply(lambda p: p.name)
    table["host"] = urlpathobj.apply(lambda p: p.parts[1])
    table["dir"] = urlpathobj.apply(lambda p: os.path.join(*p.parts[2:-1]))

    path_ftp_crawled_files_save = "/home/psakicki/aaa_FOURBI/crawl"
    table_use = ftp_files_crawler(table, secure_ftp=secure_ftp,
                                  user=user, passwd=passwd,
                                  path_ftp_crawled_files_save=path_ftp_crawled_files_save)

    #### get only the valid (true) url
    table_dl = table_use.loc[table_use['url_true'].dropna().index]

    if len(table_dl) == 0:
        log.error("no valid RINEX URL fond on the FTP server, check your inputs")

    ftp_download_frontend(table_dl['url_true'].values,
                          table_dl['outdir'].values,
                          parallel_download=parallel_download,
                          secure_ftp=secure_ftp,
                          user=user,
                          passwd=passwd)
