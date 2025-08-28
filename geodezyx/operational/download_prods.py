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
import itertools
#### Import the logger
import logging
import multiprocessing as mp
import os
import re
import shutil

import numpy as np

#### geodeZYX modules
from geodezyx import conv, utils
import geodezyx.operational.download_utils as dlutils
log = logging.getLogger('geodezyx')

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


def download_gnss_products(
    archive_dir,
    startdate,
    enddate,
    AC_names=("wum", "cod"),
    prod_types=("sp3", "clk"),
    remove_patterns=("ULA",),
    archtype="week",
    new_name_conv=True,
    parallel_download=4,
    archive_center="ign",
    mgex=False,
    repro=0,
    sorted_mode=False,
    return_also_uncompressed_files=True,
    ftp_download=False,
    dow_manu=False,
):
    """
    Download GNSS products from different IGS data centers

    Parameters
    ----------
    archive_dir : str
        the parent directory where the products will be stored.
    startdate : datetime
        the start date in regular calendar date.
    enddate : datetime
        the end date in regular calendar date..
    AC_names : tuple, optional
        the names of the wished analysis centers.
        It also control the product's lattency with the new naming convention:
        simply add it completly in the AC name e.g. IGS0OPSRAP
        The default is ("wum","cod").
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
        an alternatiove can be 'year/doy'.
    new_name_conv : bool, optional
        Also handle the new name convention. The default is True.
    parallel_download : int, optional
        control parallel download.
        The default is 4.
    archive_center : TYPE, optional
        name of the IGS's archive/data center. The default is 'ign'.
    mgex : bool, optional
        get MGEX products. The default is False.
    repro : int, optional
        get repro products. The default is 0 i.e. operational products.
    sorted_mode : bool, optional
        sort the download or not. The default is False.
    return_also_uncompressed_files : bool, optional
        in the final list output, return also already downloaded and
        uncompressedfiles. The default is True.
    ftp_download : bool, optional
        DESCRIPTION. The default is False.
    dow_manu : int or bool or None, optional
        Control the download for weekly files
        dow_manu = False, no dow manu,
        consider the converted dow from the time span, regular case
        dow_manu = None,
        no dow in the REGEX, the crawler will search only for the week
        dow_manu = 0 or 7,
        the dow in question
        The default is False.


    Returns
    -------
    list
        list of the local files's paths.

    Note
    ----
    The new naming convention has been fully adopted since GPS Week 2238-0

    to control the lattency with the new naming convention,
    simply add it completly in the AC name e.g. IGS0OPSRAP

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

    log.info("data center used : %s", archive_center)
    log.info("mgex/repro : %s/%s", mgex, repro)

    if archive_center == "cddis":
        arch_center_main = "gdc.cddis.eosdis.nasa.gov"
        arch_center_basedir = "/pub/gps/products/" + mgex_str
        ftp_download = True
        secure_ftp = True
        parallel_download = 1
        log.info("cddis as data center, FTP and no parallel download forced")

    elif archive_center == "cddis_glonass":
        arch_center_main = "cddis.gsfc.nasa.gov"
        arch_center_basedir = "/pub/glonass/products/" + mgex_str

    elif archive_center == "ign":
        arch_center_main = "igs.ign.fr"
        arch_center_basedir = "/pub/igs/products/" + mgex_str

    elif archive_center == "ign_iono":
        arch_center_main = "igs-rf.ign.fr"
        arch_center_basedir = "/pub/"

    elif archive_center == "ensg":
        arch_center_main = "igs.ensg.ign.fr"
        arch_center_basedir = "/pub/igs/products/" + mgex_str

    elif archive_center == "whu":
        arch_center_main = "igs.gnsswhu.cn"
        arch_center_basedir = "/pub/gps/products/" + mgex_str

    elif archive_center == "ign_rf":
        arch_center_main = "igs-rf.ign.fr"
        arch_center_basedir = "/pub/" + mgex_str

    elif archive_center == "ensg_rf":
        arch_center_main = "igs-rf.ensg.ign.fr"
        arch_center_basedir = "/pub/" + mgex_str

    elif archive_center == "acc_xpr_mgex_cmb":
        ##### DO NOT WORK !!!!!
        ## and this archive is not mainteed anyway (2023-01)
        arch_center_main = "http://igsacc.s3-eu-central-1.amazonaws.com/products/mgex/final/2069/igm20694.sp3.Z"
        ftp_download = False
        log.info("ACC experimental mgex combi. as data center, HTTP download forced")

    dates_list = conv.dt_range(startdate, enddate)

    wwww_dir_previous = None
    if parallel_download > 1:
        pool = mp.Pool(processes=parallel_download)

    ###################################################################
    ########### Remote file search

    potential_localfiles_list_all = []

    ### check if the pattern of the wished products are in the listed daily files
    for ipatt_tup, patt_tup in enumerate(
        list(itertools.product(dates_list, AC_names, prod_types))
    ):
        dt_cur, ac_cur, prod_cur = patt_tup
        wwww, dow = conv.dt2gpstime(dt_cur)

        #### Manage the cases of manual DOW
        if type(dow_manu) is int:
            dow = dow_manu
        elif dow_manu is None:
            dow = ""
        elif dow_manu is False:
            pass
        else:
            dow = str(dow_manu)

        log.info(
            "*** Search prods. for %s-%s, AC/prod: %s/%s", wwww, dow, ac_cur, prod_cur
        )
        wwww_dir = os.path.join(arch_center_basedir, str(wwww), repro_str)

        n_ftp_ask = 500  ## Max interrogation of the FTP server to avoid
        ## potential errors
        ## An new FTP instance is created if above it

        if np.mod(ipatt_tup, n_ftp_ask) == 0:
            log.info("Create a new FTP instance")

            ftp, ftp_obj_list = dlutils.ftp_objt_create(
                secure_ftp_inp=secure_ftp, host=arch_center_main
            )

        if wwww_dir_previous != wwww_dir or np.mod(ipatt_tup, n_ftp_ask) == 0:
            log.info("Move to: %s", wwww_dir)
            try:
                ftp.cwd(wwww_dir)
            except:
                log.warning("%s do not exists, skiping...", wwww_dir)
            files_listed_in_ftp = ftp.nlst()
            wwww_dir_previous = wwww_dir
            if len(files_listed_in_ftp) == 0:
                log.warning("no files found in directory %s", wwww_dir)

        files_remote_date_list = []

        pattern_old_nam = (
            ac_cur.lower()
            + ".*"
            + str(wwww)
            + str(dow)
            + ".*"
            + prod_cur.lower()
            + r"\..*"
        )
        files = [f for f in files_listed_in_ftp if re.search(pattern_old_nam, f)]

        pattern_new_nam = ""

        files_new_nam = []
        if new_name_conv:  ### search for new name convention

            if dow is None:
                log.error("dow == None and search for new name convention, Error ...")
                raise Exception()

            ac_newnam = ac_cur.upper()

            doy_newnam = "".join(reversed(conv.dt2doy_year(dt_cur))) + str(
                dt_cur.hour
            ).zfill(2)
            prod_newnam = prod_cur.upper()

            pattern_new_nam = utils.join_improved(
                ".*", ac_newnam, doy_newnam, prod_newnam
            )
            pattern_new_nam = ".*" + pattern_new_nam + r"\..*"

            files_new_nam = [
                f for f in files_listed_in_ftp if re.search(pattern_new_nam, f)
            ]

        log.info("Regex : %s %s", pattern_old_nam, pattern_new_nam)
        files = files + files_new_nam

        if len(files) == 0:
            log.warning("no product found for" + " %s" * len(patt_tup), *patt_tup)
            # log.warning("found files: %s",files_listed_in_ftp)

        files_remote_date_list = files_remote_date_list + files

        ### exclude some pattern
        for negpatt in remove_patterns:
            files_remote_date_list = [
                e for e in files_remote_date_list if not re.search(negpatt, e)
            ]

        ###################################################################
        ########### Download
        archive_dir_specif = dlutils.effective_save_dir_orbit(
            archive_dir, ac_cur, dt_cur, archtype
        )

        if len(files_remote_date_list) > 0:
            utils.create_dir(archive_dir_specif)

        ### Generation of the Download fct inputs
        files_remote_date_chunck = utils.chunkIt(
            files_remote_date_list, parallel_download
        )
        downld_tuples_list = []
        potential_localfiles_list = []

        if ftp_download:  ### FTP Download
            for ftpobj, Chunk in zip(ftp_obj_list, files_remote_date_chunck):
                for filchunk in Chunk:
                    potential_localfiles_list.append(
                        os.path.join(archive_dir_specif, filchunk)
                    )
                    if parallel_download == 1:
                        downld_tuples_list.append(
                            (ftpobj, filchunk, archive_dir_specif)
                        )
                    else:
                        downld_tuples_list.append(
                            (arch_center_main, wwww_dir, filchunk, archive_dir_specif)
                        )
        else:  ### HTTP download
            downld_tuples_list = itertools.product(
                [
                    "/".join(("ftp://" + arch_center_main, wwww_dir, f))
                    for f in files_remote_date_list
                ],
                [archive_dir_specif],
            )
            [
                potential_localfiles_list.append(os.path.join(archive_dir_specif, f))
                for f in files_remote_date_list
            ]

        potential_localfiles_list_all = (
            potential_localfiles_list_all + potential_localfiles_list
        )

        ### Actual Download
        if ftp_download and parallel_download == 1:
            for tup in downld_tuples_list:
                dlutils.ftp_downld_core(*tup)
        elif ftp_download and parallel_download > 1:
            _ = pool.map_async(dlutils.ftp_downloader_wo_objects, downld_tuples_list)
        elif not ftp_download and parallel_download == 1:
            for tup in downld_tuples_list:
                dlutils.downloader_wrap(tup)
        elif not ftp_download and parallel_download > 1:
            _ = pool.map(dlutils.downloader_wrap, downld_tuples_list)

    ###################################################################
    ########### Final Independent files existence check

    localfiles_lis = []
    if not return_also_uncompressed_files:
        pot_locfiles_list_use = potential_localfiles_list_all
    else:
        pot_locfiles_list_use = []
        for localfile in potential_localfiles_list_all:
            pot_compress_name_list = [localfile]
            pot_compress_name_list.append(localfile.replace(".gz", ""))
            pot_compress_name_list.append(localfile.replace(".Z", ""))
            pot_compress_name_list = list(set(pot_compress_name_list))

            pot_locfiles_list_use = pot_locfiles_list_use + pot_compress_name_list

    for pot_localfile in pot_locfiles_list_use:
        if os.path.isfile(pot_localfile):
            localfiles_lis.append(pot_localfile)

    if parallel_download > 1:
        pool.close()
    return localfiles_lis


def multi_downloader_orbs_clks_2(**kwargs):
    log.warning(
        "multi_downloader_orbs_clks_2 is a legacy alias for the newly renamed function download_gnss_products"
    )
    return download_gnss_products(**kwargs)


def orbclk_long2short_name(
    longname_filepath_in,
    rm_longname_file=False,
    center_id_last_letter=None,
    center_manual_short_name=None,
    force=False,
    dryrun=False,
    output_dirname=None,
):
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
    longname_dirname = os.path.dirname(longname_filepath_in)

    if not output_dirname:
        output_dirname = longname_dirname

    center = longname_basename[:3]

    if center_manual_short_name:
        center = center_manual_short_name
    elif center_id_last_letter:
        center_as_list = list(center)
        center_as_list[-1] = center_id_last_letter
        center = "".join(center_as_list)

    yyyy = int(longname_basename.split("_")[1][:4])
    doy = int(longname_basename.split("_")[1][4:7])

    day_dt = conv.doy2dt(yyyy, doy)

    wwww, dow = conv.dt2gpstime(day_dt)

    shortname_prefix = center.lower() + str(wwww).zfill(4) + str(dow)

    ### Type handeling
    if "SP3" in longname_basename:
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

    shortname_filepath = os.path.join(output_dirname, shortname)

    if not force and os.path.isfile(shortname_filepath):
        log.info("skip " + longname_filepath_in)
        log.info(shortname_filepath + " already exists")
        return shortname_filepath

    if not dryrun:
        log.info("renaming " + longname_filepath_in + " => " + shortname_filepath)
        shutil.copy2(longname_filepath_in, shortname_filepath)

    if rm_longname_file and not dryrun:
        log.info("remove " + longname_filepath_in)
        os.remove(longname_filepath_in)

    return shortname_filepath


#  ______                _   _                _____                                         _
# |  ____|              | | (_)              / ____|                                       | |
# | |__ _   _ _ __   ___| |_ _  ___  _ __   | |  __ _ __ __ ___   _____ _   _  __ _ _ __ __| |
# |  __| | | | '_ \ / __| __| |/ _ \| '_ \  | | |_ | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
# | |  | |_| | | | | (__| |_| | (_) | | | | | |__| | | | (_| |\ v /  __/ |_| | (_| | | | (_| |
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \_____|_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
#                                                                        __/ |
#                                                                       |___/


def multi_downloader_orbs_clks(
    archive_dir,
    startdate,
    enddate,
    calc_center="igs",
    sp3clk="sp3",
    archtype="year/doy",
    parallel_download=4,
    archive_center="ign",
    repro=0,
    sorted_mode=False,
    force_weekly_file=False,
    return_also_uncompressed_files=True,
):
    """
    Download IGS products. Can manage MGEX products too.

    Parameters
    ----------
    archive_dir : str
        Parent archive directory where files will be stored.
    startdate : datetime
        Start of the wished period.
    enddate : datetime
        End of the wished period.
    calc_center : str or list of str, optional
        Calc center can be a string or a list, describing the calc center.
        Examples: 'igs', 'grg', 'cod', 'jpl', etc. Default is "igs".
    sp3clk : str, optional
        Product type. Can handle: 'clk', 'clk_30s', 'sp3', 'snx', 'sum', 
        'erp', 'bia'. Default is "sp3".
    archtype : str, optional
        String describing how the archive directory is structured.
        Examples: 'stat', 'stat/year', 'stat/year/doy', 'year/doy', 
        'year/stat', 'week/dow/stat', etc. Default is "year/doy".
    parallel_download : int, optional
        Number of parallel downloads. Default is 4.
    archive_center : str, optional
        Server of download, "regular" IGS or MGEX. Can handle: 'cddis', 
        'cddis_mgex', 'cddis_mgex_longname', 'ign', 'ign_mgex', 
        'ign_mgex_longname', 'gfz_local'. Default is "ign".
    repro : int, optional
        Number of the IGS reprocessing (0 = routine processing). Default is 0.
    sorted_mode : bool, optional
        If False, uses the map multiprocess function so the download order 
        will be scrambled. If True, uses the apply multiprocess function so 
        the download order will be in chronological order. The scrambled 
        (False) is better because it doesn't create zombie processes. 
        Default is False.
    force_weekly_file : bool, optional
        Force download of weekly files. Default is False.
    return_also_uncompressed_files : bool, optional
        Include already downloaded and uncompressed files in the final list 
        output. Default is True.

    Returns
    -------
    localfiles_lis : list of str
        List of downloaded products paths.

    See Also
    --------
    download_gnss_products : Newer implementation of this function.

    Notes
    -----
    This function can manage MGEX products by setting the appropriate 
    archive_center
    """

    log.error(
        "multi_downloader_orbs_clks IS DISCONTINUED, use multi_downloader_orbs_clks_2"
    )
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
#                 if re.search("igs([0-9][0-9]|yy|YY)p",cc):
#                     cc = "igs" + str(curdate.year)[2:] + "p"
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
#                 savedir = dlutils.effective_save_dir_orbit(output_dir,
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
