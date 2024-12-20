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
# External modules
import datetime as dt
# from ftplib import FTP
import glob
# Import the logger
import logging
# import itertools
# import multiprocessing as mp
import os
import re

import pandas as pd

# geodeZYX modules
from geodezyx import conv
from geodezyx import utils

log = logging.getLogger(__name__)

##########  END IMPORT  ##########


def rinex_finder(
    main_dir,
    short_name=True,
    long_name=True,
    gfz_godc_name=True,
    compressed=None,
    specific_sites=[],
    start_epoch=None,
    end_epoch=None,
):
    """
    Parameters
    ----------
    main_dir : str
        main directory where the RINEXs are stored.
    short_name : bool, optional
        check if the pattern matches a short name RINEX. The default is True.
    long_name : bool, optional
        check if the pattern matches a long name RINEX. The default is True.
    gfz_godc_name: bool, optional
        check if the pattern matches a GFZ's GODC (GNSS Operational Data Center)
        internal long name RINEX. The default is True.
    compressed : bool or None
        check if the pattern matches a compressed RINEX (True) or not (False)
        if None, does not matter (return both compressed or not)
    specific_sites : list, optional
        Filter only those specific sites. The default is [].
    start_epoch and end_epoch : datetime, optional
        Filter the RINEXs between those two epochs (included)
        Can be for instance
        `start_epoch=dt.datetime(2021,1,1)` and
        `end_epoch=dt.datetime(2021,12,31)`
        The default is None.

    Returns
    -------
    Files_rnx_lis : list
        Found RINEXs list.


    Notes
    -----

    is very similar with geodetik.rinex_lister,  gins_runner.get_rinex_list,
    operational.multi_finder_rinex

    But this one is the most recent and elaborated (July 2022),
    must be used in priority !!!

    """

    Files_raw_lis, _ = utils.walk_dir(main_dir)

    Files_rnx_lis = []

    for f in Files_raw_lis:
        fbase = os.path.basename(f)
        regex_match = conv.rinex_regex_search_tester(
            fbase,
            short_name=short_name,
            long_name=long_name,
            gfz_godc_name=gfz_godc_name,
            compressed=compressed,
        )
        if regex_match:
            Files_rnx_lis.append(f)

    # SECOND FILTERING IF specific_sites LIST IS DEFINED
    if len(specific_sites) > 0:
        Files_rnx_lis_tmp = []
        for site in specific_sites:
            for rnx in Files_rnx_lis:
                rnx_bn = os.path.basename(rnx)
                if site.lower() in rnx_bn or site.upper() in rnx_bn:
                    Files_rnx_lis_tmp.append(rnx)
        Files_rnx_lis = Files_rnx_lis_tmp

    # THIRD FILTERING IF start or end epoch are defined
    if start_epoch or end_epoch:
        if not start_epoch:
            start_epoch = dt.datetime(1980, 1, 1)
        if not end_epoch:
            end_epoch = dt.datetime(2099, 1, 1)
        Files_rnx_lis_tmp = []
        Dates_rnx = [conv.rinexname2dt(rnx) for rnx in Files_rnx_lis]

        for rnx, date in zip(Files_rnx_lis, Dates_rnx):
            if start_epoch <= date and date <= end_epoch:
                Files_rnx_lis_tmp.append(rnx)
        Files_rnx_lis = Files_rnx_lis_tmp

    Files_rnx_lis = list(sorted(Files_rnx_lis))
    log.info(str(len(Files_rnx_lis)) + " RINEXs found")

    return Files_rnx_lis


def read_rinex_list_table(rnx_list_inp):
    """
    Generate a Table from a RINEX list

    Parameters
    ----------
    p : str, list or iterable
        RINEX list. Can be a Python iterable (a list) or a path to a list file
        (a string)

    Returns
    -------
    DF : DataFrame
        a RINEX Table.

    Note
    ----
    From script
    .../geodezyx_toolbox_PS_perso_scripts/IPGP_OVS/rinex_lister/rinex_list_2019_2022_mk01.py

    """

    if utils.is_iterable(rnx_list_inp):
        DF = pd.DataFrame(rnx_list_inp)
    else:
        DF = pd.read_csv(rnx_list_inp, header=None)

    DF.columns = ["path"]
    DF["name"] = DF["path"].apply(os.path.basename)
    DF["site"] = DF["name"].str[:4]
    DF["date"] = DF["name"].apply(conv.rinexname2dt)
    DF["sd"] = list(zip(*(DF["site"], DF["date"])))

    return DF


# def date_filter(DFin,strt,end):
#     return DFin[(DFin.date >= strt) & (DFin.date < end)].copy()

# def clean_w_year_path(DFin):

#     DF = DFin.copy()
#     DF["yr_path"] = DF.path.str.extract("([0-9]{4})").astype(int)
#     BOOL = DF["yr_path"] == DF["date"].dt.year
#     DF = DF[BOOL]

#     print("before/after clean_w_year_path",len(DFin),len(DF))

#     return DF


def multi_finder_rinex(
    main_dir, rinex_types=("o", "d", "d.Z", "d.z"), specific_stats=[]
):
    """
    from a main_dir, find all the rinexs in this folder and his subfolder
    (corresponding to the rinex_types)
    and return a list of the found rinexs

    Designed for RINEX-2 / short names only

    is very similar with geodetik.rinex_lister,  gins_runner.get_rinex_list,
    operational.rinex_finder

    operational.rinex_finder must be used in priority !!! (July 2022)
    """

    log.warning("multi_finder_rinex depreciated, use rinex_finder instead!!")

    files_raw_lis, _ = utils.walk_dir(main_dir)

    yylis = [
        str(e).zfill(2)
        for e in list(range(80, 100))
        + list(range(0, dt.datetime.now().year - 2000 + 1))
    ]

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

    log.info(str(len(rinex_lis)) + " RINEXs found")

    return rinex_lis


def find_IGS_products_files(
    parent_dir,
    File_type,
    ACs,
    date_start,
    date_end=None,
    recursive_search=True,
    severe=True,
    compressed="incl",
    regex_old_naming=True,
    regex_new_naming=True,
    regex_igs_tfcc_naming=True,
    add_weekly_file=False,
    add_hourly_file=False,
):
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

    add_hourly_file : bool
        Also handle ultra rapide hourly file
        Implemented only for the new naming format

    Returns
    -------
    Files_select_cumul_list : list
        List of files found

    """

    ###### Prelim Checks   ##############
    if not os.path.exists(parent_dir):
        log.error("parent directory doesn't exist")
        log.error(parent_dir)
        if severe:
            raise Exception

    ###### Time management ##############
    # work internally with datetime
    # date_start
    if type(date_start) is dt.datetime:
        date_start_ok = date_start
    else:
        date_start_ok = conv.gpstime2dt(*date_start)
    # date_end
    if not date_end:
        date_end_ok = date_start_ok
    elif type(date_end) is dt.datetime:
        date_end_ok = date_end
    else:
        date_end_ok = conv.gpstime2dt(*date_end)
    # generate time period with a while loop
    Dates_list = [date_start_ok]

    if add_hourly_file:
        deltat = dt.timedelta(seconds=3600)
    else:
        deltat = dt.timedelta(days=1)

    while Dates_list[-1] < date_end_ok:
        Dates_list.append(Dates_list[-1] + deltat)

    # manage weekly file
    Dates_wwwwd_list = [
        utils.join_improved("", *conv.dt2gpstime(d, outputtype=str)) for d in Dates_list
    ]

    if add_hourly_file:
        Dates_yyyyddd_list = [
            utils.join_improved(
                "", *reversed(conv.dt2doy_year(d)), str(d.hour).zfill(2)
            )
            for d in Dates_list
        ]
    else:
        Dates_yyyyddd_list = [
            utils.join_improved("", *reversed(conv.dt2doy_year(d))) for d in Dates_list
        ]

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
        FILE_LIST = utils.find_recursive(parent_dir, ".*", case_sensitive=False)
    else:
        FILE_LIST = glob.glob(parent_dir + "/*")

    ###### Regex Definition ##############

    def join_regex_and(L):
        return "(" + "|".join(L) + ")"

    Re_patt_big_stk = []

    # compression handeling
    if compressed == "excl":
        re_patt_comp = "$"
    elif compressed == "incl":
        re_patt_comp = r"(\.Z|\.gz|)$"
    elif compressed == "only":
        re_patt_comp = r"(\.Z|\.gz)$"
    else:
        log.error("check 'compressed' keyword (excl,incl, or only)")
        raise Exception

    if regex_old_naming:
        if ACs[0] == "all":
            re_patt_ac = r"\w{3}"
        else:
            re_patt_ac = join_regex_and([ac.lower() for ac in ACs])

        if add_weekly_file:
            Dates_wwwwd_list_4old = [
                e[:-1] + "(" + e[-1] + "|" + "7)" for e in Dates_wwwwd_list
            ]
        else:
            Dates_wwwwd_list_4old = Dates_wwwwd_list

        re_patt_date = join_regex_and(Dates_wwwwd_list_4old)
        re_patt_filtyp = join_regex_and(File_type)
        re_patt_big_old_naming = (
            re_patt_ac + re_patt_date + r"\." + re_patt_filtyp + re_patt_comp
        )
        Re_patt_big_stk.append(re_patt_big_old_naming)

    if regex_new_naming:  # search for new name convention
        if ACs[0] == "all":
            re_patt_ac = r"\W{3}"
        else:
            re_patt_ac = join_regex_and([ac.upper() for ac in ACs])
        # add _ because it can raise a conflit with the old format
        re_patt_date = join_regex_and(["_" + e for e in Dates_yyyyddd_list])
        re_patt_filtyp = join_regex_and([fil.upper() for fil in File_type])
        re_patt_big_new_naming = ".*".join(
            (re_patt_ac, re_patt_date, re_patt_filtyp + re_patt_comp)
        )
        Re_patt_big_stk.append(re_patt_big_new_naming)

    if regex_igs_tfcc_naming:
        Dates_yy_list = list(
            set(
                [
                    str(conv.gpstime2dt(int(e[0:4]), int(e[4])).year)[2:]
                    for e in Dates_wwwwd_list
                ]
            )
        )
        Dates_wwww_list = list(set([e[:-1] for e in Dates_wwwwd_list]))
        # Dates_wwww_dot_list = [e + "\." for e in Dates_wwww_list]
        re_patt_year = join_regex_and(Dates_yy_list)
        # 2x re_patt_date : because .sum doesn't the day
        re_patt_date = join_regex_and(Dates_wwwwd_list + Dates_wwww_list)
        re_patt_filtyp = r"\." + join_regex_and(File_type)

        re_patt_big_igs_tfcc_naming = (
            "igs"
            + re_patt_year
            + "P"
            + re_patt_date
            + ".*"
            + re_patt_filtyp
            + re_patt_comp
        )
        Re_patt_big_stk.append(re_patt_big_igs_tfcc_naming)

    re_patt_big = join_regex_and(Re_patt_big_stk)

    #log.info("REGEX researched: %s",re_patt_big)

    ###### Specific file search management ##############
    Files_select_list = []
    for fil in FILE_LIST:
        if re.search(re_patt_big, os.path.basename(fil), re.IGNORECASE):
            Files_select_list.append(fil)

    if len(Files_select_list) == 0:
        log.error("no products found")
        if severe:
            raise Exception

    log.info("%i files found matching REGEX %s", len(Files_select_list), re_patt_big)

    return Files_select_list
