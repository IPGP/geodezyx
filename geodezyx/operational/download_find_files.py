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
import glob
# Import the logger
import logging
import os
import re
import pandas as pd

# geodeZYX modules
from geodezyx import conv, utils

log = logging.getLogger('geodezyx')


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
    Find RINEX files in a specified directory and filter them based on various criteria.

    Parameters
    ----------
    main_dir : str
        Main directory where the RINEX files are stored.
        The directory can contain wildcards ('*', '?', etc.) and date aliases 
        ('%Y', '%j', etc.). If the main_dir contains date aliases and both 
        start_epoch and end_epoch are defined, you can indicate a more precise 
        directory (e.g., main_dir = "/path/to/data/*/%Y/%j/").
        Note: The day level is the maximum resolution for the wildcard.
    short_name : bool, optional
        Check if the pattern matches a short name RINEX. Default is True.
    long_name : bool, optional
        Check if the pattern matches a long name RINEX. Default is True.
    gfz_godc_name : bool, optional
        Check if the pattern matches a GFZ's GODC (GNSS Operational Data Center) 
        internal long name RINEX. Default is True.
    compressed : bool or None, optional
        Check if the pattern matches a compressed RINEX (True) or not (False).
        If None, does not matter (returns both compressed and uncompressed).
        Default is None.
    specific_sites : list of str, optional
        Filter only those specific sites. Default is [].
    start_epoch : datetime, optional
        Start date for filtering the RINEX files. Default is None.
    end_epoch : datetime, optional
        End date for filtering the RINEX files. Default is None.

    Returns
    -------
    files_rnx_lis : list of str
        List of found RINEX files.

    Notes
    -----
    This function is very similar to geodetik.rinex_lister, gins_runner.get_rinex_list, 
    and operational.multi_finder_rinex. However, this one is the most recent and 
    elaborated (July 2022) and should be used in priority.
    """

    # If main_dir contains a wildcard and both start_epoch and end_epoch are defined, search for files in the date range
    if "%" in main_dir and start_epoch and end_epoch:
        files_raw_lis = []
        for epo in conv.dt_range(start_epoch, end_epoch, day_step=1):
            main_dir_abs = os.path.abspath(epo.strftime(main_dir))
            files_raw_lis_epo, _ = utils.walk_dir(main_dir_abs)
            files_raw_lis.extend(files_raw_lis_epo)
    else:
        main_dir_abs = os.path.abspath(main_dir)
        files_raw_lis, _ = utils.walk_dir(main_dir_abs)

    files_rnx_lis = []

    # Filter files based on the provided naming patterns
    for f in files_raw_lis:
        fbase = os.path.basename(f)
        regex_match = conv.rinex_regex_search_tester(
            fbase,
            short_name=short_name,
            long_name=long_name,
            gfz_godc_name=gfz_godc_name,
            compressed=compressed,
        )
        if regex_match:
            files_rnx_lis.append(f)

    # Second filtering if specific_sites list is defined
    if len(specific_sites) > 0:
        files_rnx_lis_tmp = []
        for site in specific_sites:
            for rnx in files_rnx_lis:
                rnx_bn = os.path.basename(rnx)
                if utils.is_in_str(rnx_bn, site.upper(), site[0:4].lower(), site[0:4].upper()):
                    files_rnx_lis_tmp.append(rnx)
        files_rnx_lis = files_rnx_lis_tmp

    # Third filtering if start or end epoch are defined
    if start_epoch or end_epoch:
        if not start_epoch:
            start_epoch = dt.datetime(1980, 1, 1)
        if not end_epoch:
            end_epoch = dt.datetime(2099, 1, 1)
        files_rnx_lis_tmp = []
        dates_rnx = [conv.rinexname2dt(rnx) for rnx in files_rnx_lis]

        for rnx, date in zip(files_rnx_lis, dates_rnx):
            if start_epoch <= date <= end_epoch:
                files_rnx_lis_tmp.append(rnx)
        files_rnx_lis = files_rnx_lis_tmp

    # Sort the final list of RINEX files
    files_rnx_lis = list(sorted(files_rnx_lis))
    log.info(str(len(files_rnx_lis)) + " RINEXs found")

    return files_rnx_lis


def read_rinex_list_table(rnx_list_inp):
    """
    Generate a Table from a RINEX list

    Parameters
    ----------
    rnx_list_inp : str, list or iterable
        RINEX list. Can be a Python iterable (a list) or a path to a list file
        (a string)

    Returns
    -------
    df : DataFrame
        a RINEX Table.

    Note
    ----
    From script
    .../geodezyx_toolbox_PS_perso_scripts/IPGP_OVS/rinex_lister/rinex_list_2019_2022_mk01.py

    """

    if utils.is_iterable(rnx_list_inp):
        df = pd.DataFrame(rnx_list_inp)
    else:
        df = pd.read_csv(rnx_list_inp, header=None)

    df.columns = ["path"]
    df["name"] = df["path"].apply(os.path.basename)
    df["site"] = df["name"].str[:4]
    df["date"] = df["name"].apply(conv.rinexname2dt)
    df["sd"] = list(zip(*(df["site"], df["date"])))

    return df


# def date_filter(df_inp,strt,end):
#     return df_inp[(df_inp.date >= strt) & (df_inp.date < end)].copy()

# def clean_w_year_path(df_inp):

#     DF = df_inp.copy()
#     DF["yr_path"] = DF.path.str.extract("([0-9]{4})").astype(int)
#     BOOL = DF["yr_path"] == DF["date"].dt.year
#     DF = DF[BOOL]

#     print("before/after clean_w_year_path",len(df_inp),len(DF))

#     return DF


def find_igs_products_files(
        parent_dir,
        file_type,
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


    file_type : str or list of str
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
    dates_wwwwd_list = [
        utils.join_improved("", *conv.dt2gpstime(d, outputtype=str)) for d in Dates_list
    ]

    if add_hourly_file:
        dates_yyyyddd_list = [
            utils.join_improved(
                "", *reversed(conv.dt2doy_year(d)), str(d.hour).zfill(2)
            )
            for d in Dates_list
        ]
    else:
        dates_yyyyddd_list = [
            utils.join_improved("", *reversed(conv.dt2doy_year(d))) for d in Dates_list
        ]

    ###### File type / ACs management ##############

    if not utils.is_iterable(file_type):
        file_type = [file_type]
    if not utils.is_iterable(ACs):
        ACs = [ACs]

    ###### General file search management ##############
    if utils.is_iterable(parent_dir):
        file_list = parent_dir
    elif recursive_search:
        # All the files are listed first
        file_list = utils.find_recursive(parent_dir, ".*", regex=True)
    else:
        file_list = glob.glob(parent_dir + "/*")

    ###### Regex Definition ##############

    def join_regex_and(L):
        return "(" + "|".join(L) + ")"

    re_patt_big_stk = []

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
            dates_wwwwd_list_4old = [
                e[:-1] + "(" + e[-1] + "|" + "7)" for e in dates_wwwwd_list
            ]
        else:
            dates_wwwwd_list_4old = dates_wwwwd_list

        re_patt_date = join_regex_and(dates_wwwwd_list_4old)
        re_patt_filtyp = join_regex_and(file_type)
        re_patt_big_old_naming = (
                re_patt_ac + re_patt_date + r"\." + re_patt_filtyp + re_patt_comp
        )
        re_patt_big_stk.append(re_patt_big_old_naming)

    if regex_new_naming:  # search for new name convention
        if ACs[0] == "all":
            re_patt_ac = r"\W{3}"
        else:
            re_patt_ac = join_regex_and([ac.upper() for ac in ACs])
        # add _ because it can raise a conflit with the old format
        re_patt_date = join_regex_and(["_" + e for e in dates_yyyyddd_list])
        re_patt_filtyp = join_regex_and([fil.upper() for fil in file_type])
        re_patt_big_new_naming = ".*".join(
            (re_patt_ac, re_patt_date, re_patt_filtyp + re_patt_comp)
        )
        re_patt_big_stk.append(re_patt_big_new_naming)

    if regex_igs_tfcc_naming:
        dates_yy_list = list(
            set(
                [
                    str(conv.gpstime2dt(int(e[0:4]), int(e[4])).year)[2:]
                    for e in dates_wwwwd_list
                ]
            )
        )
        dates_wwww_list = list(set([e[:-1] for e in dates_wwwwd_list]))
        # Dates_wwww_dot_list = [e + "\." for e in Dates_wwww_list]
        re_patt_year = join_regex_and(dates_yy_list)
        # 2x re_patt_date : because .sum doesn't the day
        re_patt_date = join_regex_and(dates_wwwwd_list + dates_wwww_list)
        re_patt_filtyp = r"\." + join_regex_and(file_type)

        re_patt_big_igs_tfcc_naming = (
                "igs"
                + re_patt_year
                + "p"
                + re_patt_date
                + ".*"
                + re_patt_filtyp
                + re_patt_comp
        )
        re_patt_big_stk.append(re_patt_big_igs_tfcc_naming)

    re_patt_big = join_regex_and(re_patt_big_stk)

    # log.info("REGEX researched: %s",re_patt_big)

    ###### Specific file search management ##############
    files_select_list = []
    for fil in file_list:
        if re.search(re_patt_big, os.path.basename(fil), re.IGNORECASE):
            files_select_list.append(fil)

    if len(files_select_list) == 0:
        log.error("no products found")
        if severe:
            raise Exception

    log.info("%i files found matching REGEX %s", len(files_select_list), re_patt_big)

    return files_select_list

##### FUNCTION GRAVEYARD #####

# def multi_finder_rinex(
#     main_dir, rinex_types=("o", "d", "d.Z", "d.z"), specific_stats=[]
# ):
#     """
#     from a main_dir, find all the rinexs in this folder and his subfolder
#     (corresponding to the rinex_types)
#     and return a list of the found rinexs
#
#     Designed for RINEX-2 / short names only
#
#     is very similar with geodetik.rinex_lister,  gins_runner.get_rinex_list,
#     operational.rinex_finder
#
#     operational.rinex_finder must be used in priority !!! (July 2022)
#     """
#
#     log.warning("multi_finder_rinex depreciated, use rinex_finder instead!!")
#
#     files_raw_lis, _ = utils.walk_dir(main_dir)
#
#     yylis = [
#         str(e).zfill(2)
#         for e in list(range(80, 100))
#         + list(range(0, dt.datetime.now().year - 2000 + 1))
#     ]
#
#     rinex_lis = []
#
#     for f in files_raw_lis:
#         for rnxext in rinex_types:
#             for yy in yylis:
#                 if f.endswith(yy + rnxext):
#                     rinex_lis.append(f)
#
#     # CASE FOR specific_stats
#     if len(specific_stats) > 0:
#         rinex_lis2 = []
#         for stat in specific_stats:
#             for rnx in rinex_lis:
#                 if stat in os.path.basename(rnx):
#                     rinex_lis2.append(rnx)
#         rinex_lis = rinex_lis2
#
#     log.info(str(len(rinex_lis)) + " RINEXs found")
#
#     return rinex_lis
