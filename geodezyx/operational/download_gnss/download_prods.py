#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic
This sub-module of geodezyx.operational contains functions to download
gnss data and products from distant IGS servers.
it can be imported directly with:
from geodezyx import operational
The geodezyx toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License
Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

########## BEGIN IMPORT ##########
#### External modules
import itertools

#### Import the logger
import logging
import os
import pandas as pd

#### geodeZYX modules
from geodezyx import conv, utils
from . import download_utils as dlutils

log = logging.getLogger("geodezyx")

##########  END IMPORT  ##########


#  _____               _            _         _____                      _                 _
# |  __ \             | |          | |       |  __ \                    | |               | |
# | |__) | __ ___   __| |_   _  ___| |_ ___  | |  | | _____      ___ __ | | ___   __ _  __| | ___ _ __
# |  ___/ '__/ _ \ / _` | | | |/ __| __/ __| | |  | |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |/ _ \ '__|
# | |   | | | (_) | (_| | |_| | (__| |_\__ \ | |__| | (_) \ v  v /| | | | | (_) | (_| | (_| |  __/ |
# |_|   |_|  \___/ \__,_|\__,_|\___|\__|___/ |_____/ \___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|\___|_|


############################################################################
######## PRODUCTS DOWNLOADER
############################################################################


from .servers_prods import _server_select_prods


def _regex_prods(dt_cur, ac_cur, prod_cur, new_name_conv=True, dow_manu=False):
    """
    Generate a regex pattern to match GNSS product files for both old and
    new naming conventions.

    Parameters
    ----------
    dt_cur : datetime
        The date for which to build the regex.
    ac_cur : str
        Analysis center name (e.g., "wum", "cod", "IGS0OPSRAP").
    prod_cur : str
        Product type (e.g., "sp3", "clk", "erp").
    new_name_conv : bool, optional
        Also search for the new long naming convention. Default is True.
    dow_manu : int, bool, or None, optional
        Manual day-of-week override. Default is False (use computed DOW).

    Returns
    -------
    str
        A combined regex pattern (old|new convention) to match product files.
    str
        The DOW string used (for logging).
    """
    wwww, dow = conv.dt2gpstime(dt_cur)

    # Manage the cases of manual DOW
    if type(dow_manu) is int:
        dow = dow_manu
    elif dow_manu is None:
        dow = ""
    elif dow_manu is False:
        pass
    else:
        dow = str(dow_manu)

    # Old naming convention regex: e.g. wum22380.sp3.Z
    ptrn_oldnam = (
        ac_cur.lower() + ".*" + str(wwww) + str(dow) + ".*" + prod_cur.lower() + r"\..*"
    )

    # New naming convention regex: e.g. WUM0MGXFIN_20200150000_01D_15M_ORB.SP3.gz
    ptrn_newnam = ""
    if new_name_conv:
        ac_newnam = ac_cur.upper()
        doy_newnam = "".join(reversed(conv.dt2doy_year(dt_cur))) + str(
            dt_cur.hour
        ).zfill(2)
        prod_newnam = prod_cur.upper()

        ptrn_newnam = utils.join_improved(".*", ac_newnam, doy_newnam, prod_newnam)
        ptrn_newnam = ".*" + ptrn_newnam + r"\..*"

    # Combine both patterns with OR
    if ptrn_newnam:
        combined_rgx = "(" + ptrn_oldnam + "|" + ptrn_newnam + ")"
    else:
        combined_rgx = ptrn_oldnam

    return combined_rgx, str(dow)


def gen_table_prods(
    dates_list,
    ac_names,
    prod_types,
    archive_dir,
    archive_center,
    archtype="week",
    new_name_conv=True,
    mgex=False,
    repro=0,
    dow_manu=False,
):
    """
    Generate a crawl table for GNSS product downloads.

    Analogous to ``gen_crawl_table`` in ``download_rinex.py``, this function
    builds a pandas DataFrame describing all the remote files to search for,
    with their host, directory, regex pattern, and local output directory.

    Parameters
    ----------
    dates_list : list of datetime
        List of dates to download products for.
    ac_names : tuple of str
        Analysis center names.
    prod_types : tuple of str
        Product types (e.g., "sp3", "clk").
    archive_dir : str
        Parent local archive directory.
    archive_center : str
        Name of the IGS archive/data center.
    archtype : str, optional
        Local archive directory structure. Default is "week".
    new_name_conv : bool, optional
        Also handle the new naming convention. Default is True.
    mgex : bool, optional
        Get MGEX products. Default is False.
    repro : int, optional
        Reprocessing campaign number. Default is 0.
    dow_manu : int, bool, or None, optional
        Manual DOW override. Default is False.

    Returns
    -------
    pd.DataFrame
        Table with download metadata for each product file, with columns:
        date, ac, prod, outdir, host, dir, regex, protocol, crawled,
        ok_dwl, ok_loc, url_true, filename.
    """
    host, basedir, protocol, secure_ftp = _server_select_prods(
        archive_center, mgex, repro
    )

    if host is None:
        log.error("Could not resolve archive center: %s", archive_center)
        return pd.DataFrame()

    repro_str = "repro" + str(repro) + "/" if repro else ""

    table_proto = []

    for dt_cur, ac_cur, prod_cur in itertools.product(dates_list, ac_names, prod_types):
        wwww, dow = conv.dt2gpstime(dt_cur)

        # Build remote directory path
        remote_dir = os.path.join(basedir, str(wwww), repro_str)
        # Remove trailing slash for consistency
        remote_dir = remote_dir.rstrip("/")

        # Build regex pattern
        regex, dow_str = _regex_prods(dt_cur, ac_cur, prod_cur, new_name_conv, dow_manu)

        # Build local output directory
        outdir = dlutils.effective_save_dir_orbit(archive_dir, ac_cur, dt_cur, archtype)

        table_proto.append(
            (dt_cur, ac_cur, prod_cur, outdir, host, remote_dir, regex, protocol)
        )

    cols = ["date", "ac", "prod", "outdir", "host", "dir", "regex", "protocol"]
    # Create DataFrame
    table = pd.DataFrame(table_proto, columns=cols)

    # Add status columns
    table["crawled"] = False
    table["ok_dwl"] = False
    table["ok_loc"] = False
    table["url_true"] = None
    table["filename"] = ""

    return table


def download_gnss_products(
    archive_dir,
    startdate,
    enddate=None,
    ac_names=("grg", "cod"),
    prod_types=("sp3", "clk"),
    remove_patterns=("ULA",),
    archtype="week",
    new_name_conv=True,
    parallel_download=1,
    archive_center="ign",
    mgex=False,
    repro=0,
    return_also_uncompressed_files=True,
    dow_manu=False,
    quiet_mode=False,
    force=False,
    path_ftp_crawled_files_save=None,
    path_ftp_crawled_files_load=None,
    skip_crawl=False,
    path_all_ftp_files_save=None,
):
    """
    Download GNSS products from different IGS data centers.

    Uses a unified table-based crawl approach shared with
    ``download_gnss_rinex``: first builds a crawl table, then crawls
    the FTP server to locate available files (checking local existence
    first), and finally downloads the matched files.

    Parameters
    ----------
    archive_dir : str
        the parent directory where the products will be stored.
    startdate : datetime or list of datetime
        the start date in regular calendar date.
    enddate : datetime or None
        the end date in regular calendar date.
        If None, only the startdate is be considered as a date list.
    ac_names : tuple, optional
        the names of the wished analysis centers.
        It also control the product's lattency with the new naming convention:
        simply add it completly in the AC name e.g. IGS0OPSRAP
        The default is ("grg","cod").
    prod_types : tuple, optional
        the wished products.
        The default is ("sp3","clk").
    remove_patterns : tuple, optional
        the patterns you want to exclude.
        The default is ("ULA",).
    archtype : str, optional
        structure of the local archive sub-directories.
        see `effective_save_dir_orbit` function for more details.
        The default is 'week'.
        an alternative can be 'year/doy'.
    new_name_conv : bool, optional
        Also handle the new name convention. The default is True.
    parallel_download : int, optional
        control parallel download.
        The default is 4.
    archive_center : str, optional
        name of the IGS's archive/data center. The default is 'ign'.
    mgex : bool, optional
        get MGEX products. The default is False.
    repro : int, optional
        get repro products. The default is 0 i.e. operational products.
    sorted_mode : bool, optional
        sort the download or not. The default is False.
    return_also_uncompressed_files : bool, optional
        in the final list output, return also already downloaded and
        uncompressed files. The default is True.
    ftp_download : bool, optional
        kept for backward compatibility. The default is False.
    dow_manu : int or bool or None, optional
        Control the download for weekly files.
        dow_manu = False: use computed DOW (regular case).
        dow_manu = None: no DOW in regex, search only by week.
        dow_manu = 0 or 7: specific DOW.
        The default is False.
    quiet_mode : bool, optional
        List the available products without downloading them.
        Useful with ``path_ftp_crawled_files_save``. Default is False.
    force : bool, optional
        Force re-download even if files exist locally. Default is False.
    path_ftp_crawled_files_save : str, optional
        Save the crawled files table as CSV at this path. Default is None.
    path_ftp_crawled_files_load : str, optional
        Load a previously saved crawl table from this path. Default is None.
    skip_crawl : bool, optional
        Skip the FTP crawl step (use loaded table as-is). Default is False.
    path_all_ftp_files_save : str, optional
        Save ALL remote files found on the FTP server as CSV. Default is None.

    Returns
    -------
    list
        list of the local files's paths.

    Note
    ----
    The new naming convention has been fully adopted since GPS Week 2238-0.

    To control the lattency with the new naming convention,
    simply add it completly in the AC name e.g. IGS0OPSRAP.
    """

    if not utils.is_iterable(remove_patterns):
        remove_patterns = list(remove_patterns)

    # Handle CDDIS specifics
    if archive_center == "cddis":
        parallel_download = 1
        log.info("cddis as data center, no parallel download forced")

    log.info("data center used : %s", archive_center)
    log.info("mgex/repro : %s/%s", mgex, repro)

    # Build dates list
    if enddate:
        dates_list = conv.dt_range(startdate, enddate)
    elif utils.is_iterable(startdate):
        dates_list = startdate
    else:
        dates_list = [startdate]

    ###################################################################
    ########### Build or load the crawl table
    if path_ftp_crawled_files_load:
        table_ini = pd.read_csv(path_ftp_crawled_files_load)
    else:
        table_ini = gen_table_prods(
            dates_list,
            ac_names,
            prod_types,
            archive_dir,
            archive_center,
            archtype=archtype,
            new_name_conv=new_name_conv,
            mgex=mgex,
            repro=repro,
            dow_manu=dow_manu,
        )

    if len(table_ini) == 0:
        log.error("No product entries generated for the given criteria.")
        return []

    log.info("Crawl table: %d entries generated", len(table_ini))

    ###################################################################
    ########### FTP Crawl (with local file check)
    if skip_crawl:
        table_use = table_ini
        files_all = pd.Series([], dtype=str)
        files_loc = pd.Series([], dtype=str)
    else:
        table_use, files_all, files_loc = dlutils.crawl_ftp_files(
            table_ini,
            sftp="auto",
            path_ftp_crawled_files_save=path_ftp_crawled_files_save,
            path_all_ftp_files_save=path_all_ftp_files_save,
            force=force,
            all_files_mode=True,  # products: get ALL matches (old+new conv)
            exclude_patterns=list(remove_patterns) if remove_patterns else None,
        )

    ###################################################################
    ########### Download
    # Get only the valid (non-null) URLs
    table_dl = table_use.loc[table_use["url_true"].dropna().index]

    out_tup_lis = []

    if len(table_dl) == 0:
        log.error(
            "no valid product URL found/selected on the FTP server, check your inputs"
        )
    elif quiet_mode:
        log.warning("quiet mode, no download was performed")
    else:
        # Expand semicolon-joined URLs (produced by all_files_mode) into individual URLs.
        # Each row in table_dl may hold multiple ';'-separated URLs for a given
        # date/AC/product combination.  ftp_downld_front / ftp_downld_mono only
        # handles one URL at a time, so we must split them out here.
        url_list_exp = []
        outdir_list_exp = []
        protocol_list_exp = []
        for url_val, outdir_val, prot_val in zip(
            table_dl["url_true"].values,
            table_dl["outdir"].values,
            (table_dl["protocol"] == "sftp").values,
        ):
            parts = str(url_val).split(";")
            url_list_exp.extend(parts)
            outdir_list_exp.extend([outdir_val] * len(parts))
            protocol_list_exp.extend([prot_val] * len(parts))

        out_tup_lis = dlutils.ftp_downld_front(
            url_list_exp,
            outdir_list_exp,
            parallel_download=parallel_download,
            secure_ftp=protocol_list_exp,
            force=force,
        )

    ###################################################################
    ########### Collect results
    # Add local paths to the output
    if len(files_loc) > 0:
        loc_tup_lis = [(f, True) for f in files_loc]
        out_tup_lis_fin = loc_tup_lis + out_tup_lis
    else:
        out_tup_lis_fin = out_tup_lis

    log.info(
        "Product files fetched: total: %d, downloaded: %d, already here: %d",
        len(out_tup_lis_fin),
        len(out_tup_lis),
        len(files_loc),
    )

    # Return plain list of paths for backward compatibility
    localfiles_lis = []
    for item in out_tup_lis_fin:
        if isinstance(item, tuple):
            fpath, ok = item
            if ok and os.path.isfile(fpath):
                localfiles_lis.append(fpath)
        elif isinstance(item, str) and os.path.isfile(item):
            localfiles_lis.append(item)

    # Also check for uncompressed variants if requested
    if return_also_uncompressed_files:
        extra_files = []
        for fpath in localfiles_lis:
            for variant in [fpath.replace(".gz", ""), fpath.replace(".Z", "")]:
                if variant != fpath and os.path.isfile(variant):
                    extra_files.append(variant)
        localfiles_lis = list(set(localfiles_lis + extra_files))

    localfiles_lis = list(sorted(localfiles_lis))

    return localfiles_lis


def multi_downloader_orbs_clks_2(**kwargs):
    log.warning(
        "multi_downloader_orbs_clks_2 is a legacy alias for the newly renamed "
        "function download_gnss_products"
    )
    return download_gnss_products(**kwargs)
