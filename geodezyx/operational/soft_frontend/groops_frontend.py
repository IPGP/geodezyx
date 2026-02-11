#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 10:54:49 2023

@author: psakicki
"""

import datetime as dt
#### Import the logger
import logging
import os
import time
import urllib
from threading import Thread

import numpy as np

from geodezyx import utils, conv, operational

log = logging.getLogger('geodezyx')


def log_subprocess(pipe, logger=None, file=None, file2=None):
    """
    Intern fuction for subprocess_frontend2
    to write the stdout/err in the console logger + a logfile
    """
    for line in iter(pipe.readline, b""):  # b'\n'-separated lines
        line_clean = line.decode().strip()
        if file:
            file.write(line_clean + "\n")
        if file2:
            file2.write(line_clean + "\n")
        if logger:
            logger(line_clean)
    return


def subprocess_frontend2(
    cmd_in,
    save_log=True,
    log_dir=None,
    log_name_out="out.log",
    log_name_err="err.log",
    logname_timestamp=True,
    err_also_in_outfile=True,
    logger_objt_level_out=None,
    logger_objt_level_err=None,
):
    """
    A generic frontend to run an extern command through subprocess
    and write the stdout and stderr outputs in log files.

    Parameters
    ----------
    cmd_in : str
        The subprocess command.
    save_log : str, optional
        export as log the stdout/stderr in files. The default is True.
    log_dir : str, optional
        directory where the logs will be stored. The default is None.
    log_name_out : str, optional
        filename of the stdout log. The default is "out.log".
    log_name_err : str, optional
        filename of the stderr log. The default is "err.log".
    logname_timestamp : str, optional
        add a timestamp as prefix. The default is True.
    logger_objt_level_out : method of a Logger object, optional
        set the logger level of the stdout messages.
        can be logger.info or logger.debug for instance.
        The default is None (logger.info per default then).
    logger_objt_level_err : method of a Logger object, optional
        set the logger level of the stderr messages.
        can be logger.error or logger.critical for instance.
        The default is None (logger.error per default then).

    Returns
    -------
    exitcode : int
        the exit code of the command.

    Notes
    -----
    Inspired by:
    https://stackoverflow.com/questions/21953835/run-subprocess-and-print-output-to-logging
    And for the threading:
    https://stackoverflow.com/questions/6809590/merging-a-python-scripts-subprocess-stdout-and-stderr-while-keeping-them-disti
    """

    from subprocess import Popen, PIPE

    now = utils.get_timestamp()

    #### manage the paths of the output logs
    if save_log:
        if not log_dir:
            log_dir = os.getcwd()

        if logname_timestamp:
            prefix = now + "_"
        else:
            prefix = ""

        err_file = open(log_dir + "/" + prefix + log_name_err, "w+")
        out_file = open(log_dir + "/" + prefix + log_name_out, "w+")

        if err_also_in_outfile:
            err_in_outfile = out_file
        else:
            err_in_outfile = None

    else:
        out_file = None
        err_file = None
        err_in_outfile = None

    ## set logger level
    if logger_objt_level_out:
        loggerout = logger_objt_level_out
    else:
        loggerout = log.info

    if logger_objt_level_err:
        loggererr = logger_objt_level_err
    else:
        loggererr = log.error

    ## arguments must be splitted
    if type(cmd_in) is str:
        cmd_split = cmd_in.split()
    else:
        cmd_split = cmd_in

    ####### Run the command here  #####################
    process = Popen(cmd_split, stdout=PIPE, stderr=PIPE)
    ##################################################

    ### Get out/err simultaneously
    stdout_thread = Thread(
        target=log_subprocess, args=(process.stdout, loggerout, out_file)
    )
    stderr_thread = Thread(
        target=log_subprocess,
        args=(process.stderr, loggererr, err_file, err_in_outfile),
    )

    stderr_thread.start()
    stdout_thread.start()

    while stdout_thread.is_alive() or stderr_thread.is_alive():
        pass  ### do nothing while the threads are runing

    exitcode = process.wait()  # 0 means success

    if save_log:
        out_file.close()
        err_file.close()

    return exitcode


def groops_basic_runner(
    xml_cfg_path="",
    global_var_dict=dict(),
    xml_var_dict=dict(),
    dry_run=False,
    verbosity="ERROR",
    log_dir=None,
    groops_bin_path="/opt/softs_gnss/groops/bin/groops",
):
    """


    Parameters
    ----------
    xml_cfg_path : str, optional
        the XML config file for GROOPS. The default is "".
    global_var_dict : dict, optional
        A dictionnary to change the global variables values, like:
        global_var_dict["global_var_name"] = new_value
        The default is dict().
    xml_var_dict : dict, optional
        A dictionnary to change the values in the config XML file.
        Not implemented yet.
        The default is dict().
    dry_run : bool, optional
        If True print the command but do not run it . The default is False.
    verbosity : str, optional
        verbosity level of the console logger
        use keywords CRITICAL, ERROR, WARNING, INFO, DEBUG
        The default is 'ERROR'.
    log_dir : str, optional
        If provided, directory where the logs are stored. The default is None.
    groops_bin_path : TYPE, optional
        Path of the GROOPS bin.
        The default is '/opt/softs_gnss/groops/bin/groops'.

    Returns
    -------
    None.

    """

    global_args_str = ""
    for key, val in global_var_dict.items():
        global_args_str = global_args_str + " ".join(
            (" --global", key + "=" + str(val))
        )

    command = " ".join((groops_bin_path, global_args_str, xml_cfg_path))
    log.info("groops command:")
    log.info("%s", command)

    xml_cfg_bn = os.path.basename(xml_cfg_path)

    loglevel_prev = logging.getLevelName(log.level)
    log.setLevel(verbosity)

    if not log_dir:
        save_log = False
    else:
        save_log = True

    if not dry_run:
        subprocess_frontend2(
            command,
            save_log=True,
            log_dir=log_dir,
            log_name_out=xml_cfg_bn + ".out.log",
            log_name_err=xml_cfg_bn + ".err.log",
            logname_timestamp=True,
            logger_objt_level_out=log.debug,
        )

    log.setLevel(loglevel_prev)
    return


def vmf_tropo_downloader(
    output_dir,
    startdate=dt.datetime(2019, 1, 1),
    enddate=dt.datetime(2019, 1, 1),
    model="VMF3",
    version="OP",
):

    Dates_range = conv.dt_range(startdate, enddate)

    url_root = "https://vmf.geo.tuwien.ac.at"
    dir_root = output_dir
    subdir = "trop_products/GRID/1x1/"

    Files_hour_list = []
    for currdate in Dates_range:
        subdir_date = os.path.join(
            subdir, model, model + "_" + version, str(currdate.year)
        )

        url_date = os.path.join(url_root, subdir_date)
        dir_date = os.path.join(dir_root, subdir_date)
        utils.create_dir(dir_date)

        date_str = conv.dt2str(currdate, "%Y%m%d")

        for hourext in ("H00", "H06", "H12", "H18"):
            filename_hour = model + "_" + date_str + "." + hourext
            url_hour = os.path.join(url_date, filename_hour)
            outfile_hour = os.path.join(dir_date, filename_hour)

            log.info("download %s to %s", url_hour, outfile_hour)

            Files_hour_list.append(outfile_hour)
            if not utils.empty_file_check(outfile_hour):
                log.info("%s exists, download skipped", outfile_hour)
                continue
            else:
                DownObjt = urllib.request.urlretrieve(url_hour, outfile_hour)
    return Files_hour_list


# %%


def _search_xml_groops_global(xmlpath, global_label):
    from xml.etree import ElementTree as et

    xmltree = et.parse(xmlpath)

    for element in xmltree.getroot().find("global"):
        if element.attrib["label"] == global_label:
            elementok = element
            return elementok.text
    return None


###############################################################################
###############################################################################


def groops_ppp_full_runner(
    rinex_path,
    project_name,
    igs_ac_10char,
    cfg_files_dict,
    log_root_dir,
    vmf_tropo_root_dir,
    prods_gnss_root_dir,
    cfg_files_root_dir,
    sitelogs_root_dir,
    groops_bin_path="/opt/softs_gnss/groops/bin/groops",
    verbosity="ERROR",
):
    """
    High level function to run a GROOPS's PPP job
    This function download IGS's products, convert them, and run the PPP job
    see the Notes below for more details

    Parameters
    ----------
    rinex_path : str
        the path of the RINEX file to process.
    project_name : str
        a personalized name for your processing project
    igs_ac_10char : str
        the 10 char. ID for the IGS AC you want to use
        e.g. 'COD0OPSFIN'
        can handle operational (OPS) and MGEX (MGX) lines
    cfg_files_dict : dict
        a dictionary controlling the conversion/processing steps
        and the corresponding config files.
    log_root_dir : str
        directory path where the frontend logs will be written
    vmf_tropo_root_dir : str
        directory path where the VMF3/ECMWF grids will be stored
        (auto download)
    prods_gnss_root_dir : str
        directory path where the IGS products (not converted) will be stored
        (auto download)
    cfg_files_root_dir : str
        directory path where the config files are stored.
    sitelogs_root_dir : str
        directory path where the sitelogs files are stored.
    groops_bin_path : str, optional
        Path of the GROOPS bin.
        The default is '/opt/softs_gnss/groops/bin/groops'.
    verbosity : str, optional
        verbosity level of the console logger for the
        GROOPS's messages
        use keywords CRITICAL, ERROR, WARNING, INFO, DEBUG
        The default is 'ERROR'.
        This verbosity level does not apply to this function

    Note
    ----
    **the config files must be checked and edited manually**
    to fit your environnement config
    it is the config files which contains the most useful parameters
    use `groopsGui` to help you

    prototype for config files are in:
    `.../geodezyx/000_exemples/groops_frontend/configfiles/`

    do not forget to update on a regular basis GROOPS's `data` folder:
    ` https://ftp.tugraz.at/outgoing/ITSG/groops/data.zip`

    Returns
    -------
    None.

    """

    ###############################################################################
    ######## Set python fct variables

    #### preliminary test, GROOPS bin is here?
    if not os.path.isfile(groops_bin_path):
        log.critical("GROOPS bin is not here: %s", groops_bin_path)
        return

    #### Internal debug variables
    log.setLevel(logging.INFO)
    debug_sleep_time = 1
    dry_run = False

    #### add a slash to avoid crash
    cfg_files_root_dir = cfg_files_root_dir + "/"

    prods_gnss_dir = os.path.join(prods_gnss_root_dir, igs_ac_10char)

    #### Determine if products are MGEX or not
    if igs_ac_10char[4:7] == "MGX":
        mgex = True
    else:
        mgex = False
    #### Determine if products are repro or not
    if (igs_ac_10char[4] == "R") and (igs_ac_10char[6] == "3"):
        repro = 3
        mgex = False
    else:
        repro = 0

    #### Define project name
    project_name_use = project_name + "_" + igs_ac_10char

    #### Define rinex-linked variables
    site4char = os.path.basename(rinex_path)[:4].lower()
    date_rnx = conv.rinexname2dt(rinex_path)
    date_rnx_mjd = int(conv.dt2mjd(date_rnx))
    date_rin_ymd = conv.dt2str(date_rnx, "%Y-%m-%d")
    date_rin_doy = utils.join_improved("-", *list(reversed(conv.dt2doy_year(date_rnx))))
    date_rin_dow = utils.join_improved("-", *conv.dt2gpstime(date_rnx))

    #### Create associated log directory
    log_dir = os.path.join(
        log_root_dir,
        "_".join((utils.get_timestamp(), project_name_use, site4char, date_rin_ymd)),
    )
    utils.create_dir(log_dir)

    ###############################################################################
    ######## Initial print
    log.info("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    log.info("****** GROOPS RUNNER Start ****************************************")
    log.info("** RINEX: %s", os.path.basename(rinex_path))
    log.info("** AC: %s", igs_ac_10char)
    log.info("** date: %s, MJD: %s", date_rin_ymd, date_rnx_mjd)
    log.info("** year-doy: %s, GPS week-dow: %s", date_rin_doy, date_rin_dow)

    ###############################################################################
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ######## VMF TROPO

    ###############################################################################
    ######## Download Tropo
    log.info("****** VMF troposphere grids download *****************************")

    VMF_files = vmf_tropo_downloader(
        vmf_tropo_root_dir, startdate=date_rnx - dt.timedelta(days=0), enddate=date_rnx
    )

    ###############################################################################
    ######## Convert Tropo
    log.info("****** VMF troposphere grids conversion ***************************")

    xml_cfg_path = cfg_files_root_dir + cfg_files_dict["convTropo"]

    global_var_dict = dict()
    global_var_dict["timeStart"] = date_rnx_mjd
    global_var_dict["timeEnd"] = date_rnx_mjd
    global_var_dict["groopsInpVmfGridFile1"] = VMF_files[-4]
    global_var_dict["groopsInpVmfGridFile2"] = VMF_files[-3]
    global_var_dict["groopsInpVmfGridFile3"] = VMF_files[-2]
    global_var_dict["groopsInpVmfGridFile4"] = VMF_files[-1]

    groops_basic_runner(
        xml_cfg_path,
        global_var_dict,
        log_dir=log_dir,
        dry_run=dry_run,
        verbosity=verbosity,
        groops_bin_path=groops_bin_path,
    )

    ###############################################################################
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ######## PRODUCTS

    ###############################################################################
    ######## Check if converted products exists
    log.info("****** Converted GNSS products existance check ********************")
    #### XMLSEARCH
    # conv_prod_root_dir = "/home/ovsgnss/020_CALC/groops_process/021_conv_igs_prods"
    xml_cfg_path = cfg_files_root_dir + cfg_files_dict["convProds_orbit"]
    conv_prod_root_dir = _search_xml_groops_global(xml_cfg_path, "outIgsProdsDir")
    conv_prod_root_dir = os.path.dirname(conv_prod_root_dir)
    conv_prod_dir = os.path.join(conv_prod_root_dir, igs_ac_10char, date_rin_ymd)

    All_conv_prods = utils.find_recursive(conv_prod_dir, "*dat")
    Clk_conv_prods = utils.find_recursive(conv_prod_dir, "clock*dat")
    Orb_conv_prods = utils.find_recursive(conv_prod_dir, "orbit*dat")

    N_conv_prods = np.array([len(l) for l in [Clk_conv_prods, Orb_conv_prods]])
    if np.all(
        N_conv_prods > 12
    ):  ## 12 is arbitrary and all the prods must fit the test
        download_prods = False
        log.info(
            "Products download/conversion skipped, %s converted files found in %s",
            len(All_conv_prods),
            conv_prod_dir,
        )
    else:
        download_prods = True
        log.info(
            "Products download/conversion will start, no files found in %s",
            conv_prod_dir,
        )

    ###############################################################################
    ######## Download Products

    if download_prods:
        log.info("****** GNSS products download *********************************")

        download_fct = operational.multi_downloader_orbs_clks_2

        Prod_types = ["sp3", "clk", "bia", "obx"]
        AC_names = [igs_ac_10char]
        Prods_out_list = download_fct(
            prods_gnss_dir,
            date_rnx,
            date_rnx,
            AC_names,
            Prod_types,
            archtype="week",
            archive_center="cddis",
            parallel_download=1,
            repro=repro,
            mgex=mgex,
        )

        if not Prods_out_list:
            log.warning("No downloaded products found locally, abort")
            raise Exception

        ###### Check the downloaded products (IGS AC ID)
        Igs_ac_10char_list_out = [os.path.basename(e)[:10] for e in Prods_out_list]
        Igs_ac_10char_list_out = list(set(Igs_ac_10char_list_out))

        if len(Igs_ac_10char_list_out) > 1:
            log.warning(
                "Several downloaded product types found: %s", Igs_ac_10char_list_out
            )

        if Prods_out_list and (igs_ac_10char != Igs_ac_10char_list_out[0]):
            log.warning(
                "Product type downloaded != requested: %s",
                Igs_ac_10char_list_out[0],
                igs_ac_10char,
            )

    ###############################################################################
    ######## Convert Products: clocks
    if download_prods:
        log.info("****** GNSS products conversion: clocks *******************")

        xml_cfg_path = cfg_files_root_dir + cfg_files_dict["convProds_clock"]

        global_var_dict = dict()
        global_var_dict["timeStart"] = date_rnx_mjd
        global_var_dict["timeEnd"] = date_rnx_mjd + 1
        global_var_dict["igsAC10Char"] = igs_ac_10char
        global_var_dict["inpIgsProdsDir"] = prods_gnss_dir
        groops_basic_runner(
            xml_cfg_path,
            global_var_dict,
            log_dir=log_dir,
            dry_run=dry_run,
            verbosity=verbosity,
            groops_bin_path=groops_bin_path,
        )

        time.sleep(debug_sleep_time)

    ###############################################################################
    ######## Convert Products: bias
    if download_prods:
        log.info("****** GNSS products conversion: bias *******************")

        xml_cfg_path = cfg_files_root_dir + cfg_files_dict["convProds_bias"]

        groops_basic_runner(
            xml_cfg_path,
            global_var_dict,
            log_dir=log_dir,
            dry_run=dry_run,
            verbosity=verbosity,
            groops_bin_path=groops_bin_path,
        )

        time.sleep(debug_sleep_time)

    ###############################################################################
    ######## Convert Products: orbits
    if download_prods:
        log.info("****** GNSS products conversion: orbits *******************")

        xml_cfg_path = cfg_files_root_dir + cfg_files_dict["convProds_orbit"]

        groops_basic_runner(
            xml_cfg_path,
            global_var_dict,
            log_dir=log_dir,
            dry_run=dry_run,
            verbosity=verbosity,
            groops_bin_path=groops_bin_path,
        )

        time.sleep(debug_sleep_time)

    ###############################################################################
    ######## Convert Products: alternative attitude
    if download_prods:
        log.info("****** GNSS products conversion: alternative attitude *****")

        xml_cfg_path = cfg_files_root_dir + cfg_files_dict["convProds_alt_attitude"]

        groops_basic_runner(
            xml_cfg_path,
            global_var_dict,
            log_dir=log_dir,
            dry_run=dry_run,
            verbosity=verbosity,
            groops_bin_path=groops_bin_path,
        )

        time.sleep(debug_sleep_time)

    ###############################################################################
    ######## Convert Products: attitude
    if download_prods:
        log.info("****** GNSS products conversion: attitude *****")

        xml_cfg_path = cfg_files_root_dir + cfg_files_dict["convProds_attitude"]

        groops_basic_runner(
            xml_cfg_path,
            global_var_dict,
            log_dir=log_dir,
            dry_run=dry_run,
            verbosity=verbosity,
            groops_bin_path=groops_bin_path,
        )

        time.sleep(debug_sleep_time)

    ###############################################################################
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ######## RINEX & STATION INFO

    ###############################################################################
    ######## convert sitelog > station info

    log.info("****** sitelog > station info conversion **************************")

    xml_cfg_path = cfg_files_root_dir + cfg_files_dict["convSite_sitelog"]

    sitelogs = utils.find_recursive(
        sitelogs_root_dir, ".*" + site4char + ".*log", case_sensitive=False
    )
    if not sitelogs:
        log.error("no sitelog found, sitelog conversion step is skipped")
    else:
        sitelog_path = sitelogs[-1]
        log.info("sitelog used for conversion: %s", sitelog_path)

        global_var_dict = dict()
        # global_var_dict["timeStart"] = date_rnx_mjd
        # global_var_dict["timeEnd"]   = date_rnx_mjd + 1
        global_var_dict["site4Char"] = site4char
        global_var_dict["inpSitelogMainDir"] = os.path.dirname(sitelog_path)
        global_var_dict["sitelogFile"] = os.path.basename(sitelog_path)

        groops_basic_runner(
            xml_cfg_path,
            global_var_dict,
            log_dir=log_dir,
            dry_run=dry_run,
            verbosity=verbosity,
            groops_bin_path=groops_bin_path,
        )

    ###############################################################################
    ######## Convert RINEX

    log.info("****** RINEX observation conversion *******************************")

    xml_cfg_path = cfg_files_root_dir + cfg_files_dict["convSite_rnxObs"]

    global_var_dict = dict()
    global_var_dict["timeStart"] = date_rnx_mjd
    global_var_dict["timeEnd"] = date_rnx_mjd + 1
    global_var_dict["site4Char"] = site4char
    global_var_dict["inpGnssObsRnxFile"] = rinex_path

    groops_basic_runner(
        xml_cfg_path,
        global_var_dict,
        log_dir=log_dir,
        dry_run=dry_run,
        verbosity=verbosity,
        groops_bin_path=groops_bin_path,
    )

    ###############################################################################
    ######## Edit station list

    log.info("****** Station list edition ***************************************")

    #### XMLSEARCH
    xml_cfg_path = cfg_files_root_dir + cfg_files_dict["gnssProcessing"]
    station_list_dir = _search_xml_groops_global(xml_cfg_path, "groopsStationListDir")
    station_list_file = _search_xml_groops_global(xml_cfg_path, "groopsStationListFile")
    station_list_path = os.path.join(station_list_dir, station_list_file)

    F = open(station_list_path, "w+")
    F.write(site4char)
    F.close()

    log.info("site %s written in station list file %s", site4char, station_list_path)

    ###############################################################################
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ######## Run Processing

    log.info("****** Processing run *********************************************")

    xml_cfg_path = cfg_files_root_dir + cfg_files_dict["gnssProcessing"]

    global_var_dict = dict()
    global_var_dict["igsAC10Char"] = igs_ac_10char
    global_var_dict["projectName"] = project_name_use
    global_var_dict["timeStart"] = date_rnx_mjd
    global_var_dict["timeEnd"] = date_rnx_mjd + 1

    groops_basic_runner(
        xml_cfg_path,
        global_var_dict,
        log_dir=log_dir,
        dry_run=dry_run,
        verbosity=verbosity,
        groops_bin_path=groops_bin_path,
    )

    log.info("****** GROOPS RUNNER End ******************************************")
    log.info("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

    return


### ===========================================================================
### ===========================================================================
