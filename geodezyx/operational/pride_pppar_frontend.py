#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 18:01:49 2024

@author: psakic
"""

import datetime as dt

#### Import the logger
import logging
import multiprocessing as mp
import os
import re
import shutil
import subprocess

import hatanaka

from geodezyx import conv
from geodezyx import files_rw
from geodezyx import operational
from geodezyx import utils

log = logging.getLogger(__name__)

import re


def remove_regex_reserved_characters(input_string):
    """
    Removes all REGEX reserved characters from the input string.

    Parameters
    ----------
    input_string : str
        The string to be cleaned of REGEX reserved characters.

    Returns
    -------
    str
        The cleaned string with all REGEX reserved characters removed.
    """
    # Define a pattern that matches all REGEX reserved characters
    regex_reserved_chars = r"[.*+?^${}()|[\]\\]"

    # Use re.sub to replace all occurrences of the reserved characters with an empty string
    cleaned_string = re.sub(regex_reserved_chars, "", input_string)

    return cleaned_string

def run_command(command):
    """
    Runs a shell command and captures both stdout and stderr.

    Parameters
    ----------
    command : str
        The shell command to be executed.

    Notes
    -----
    This function uses subprocess.Popen to run the command in a new process.
    It continuously reads and prints stdout and stderr until the process finishes.
    The function prints the return code of the process once it completes.
    """
    # Run the command and capture both stdout and stderr
    process = subprocess.Popen(
        [command],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        executable="/bin/bash",
        shell=True,
    )

    # Continuously read and print stdout and stderr
    while True:
        # Read a line from stdout
        stdout_line = process.stdout.read().decode("utf-8")
        if stdout_line:
            print("STDOUT: %s", stdout_line.strip())
        # Read a line from stderr
        # stderr_line = process.stderr.read1().decode("utf-8")
        # if stderr_line:
        #    print(f"STDERR: {stderr_line.strip()}", end='', flush=True)

        # Check if the process has finished
        return_code = process.poll()
        if return_code is not None:
            print("RETURN CODE: %s", return_code)
            break


def dl_brdc_pride_pppar(prod_parent_dir, date_list):
    """
    Downloads BRDC (Broadcast Ephemeris) files for PRIDE PPPAR from a given directory and date list.

    Parameters
    ----------
    prod_parent_dir : str
        The parent directory where the products are stored.
    date_list : iterable of datetime
        The list of dates for which the BRDC files are to be downloaded.

    Returns
    -------
    list
        A list of downloaded BRDC files.

    Notes
    -----
    This function first rounds the dates in the date list to the nearest day and removes duplicates.
    It then downloads the BRDC files for each unique date using the operational.download_gnss_rinex function.
    The downloaded files are appended to a list which is returned at the end.
    """
    brdc_lis = []
    ######### BROADCAST
    date_list_uniq = [conv.round_dt(d, "1D", "floor") for d in date_list]
    date_list_uniq = sorted(list(set(date_list_uniq)))

    for date in date_list_uniq:
        brdc = operational.download_gnss_rinex(
            {"nav_rt": ["BRDC"]},
            prod_parent_dir,
            date,
            date,
            archtype="year/doy",
            parallel_download=1,
            force=False,
        )
        brdc_lis.append(brdc)

    return brdc_lis


def dl_prods_pride_pppar(prod_parent_dir, date_list, prod_ac_name):
    """
    Downloads GNSS products for PRIDE PPPAR from a given directory and date list.

    Parameters
    ----------
    prod_parent_dir : str
        The parent directory where the products are stored.
    date_list : iterable of datetime
        The list of dates for which the GNSS products are to be downloaded.
    prod_ac_name : str
        The name of the analysis center providing the products.

    Returns
    -------
    list
        A list of downloaded GNSS products.

    Notes
    -----
    This function downloads various GNSS products such as orbits, clocks, biases, etc.
    It iterates over specified data centers and attempts to download the products.
    If at least 5 products are found, the function stops further downloads.
    """

    dl_prods_fct = operational.download_gnss_products

    ######### ORBITS CLOCKS ETC...
    for data_center in ("ign", "whu"):  ##'whu'
        if "MGX" in prod_ac_name:
            mgex = True
        else:
            mgex = False

        prods = dl_prods_fct(
            prod_parent_dir,
            min(date_list),
            max(date_list),
            AC_names=(prod_ac_name,),
            prod_types=("sp3", "clk", "bia", "obx", "erp"),
            remove_patterns=("ULA",),
            archtype="year/doy",
            new_name_conv=True,
            parallel_download=1,
            archive_center=data_center,
            mgex=mgex,
            repro=0,
            sorted_mode=False,
            return_also_uncompressed_files=True,
            ftp_download=False,
            dow_manu=False,
        )

        if len(prods) >= 5:
            log.info("enougth products found: %s", len(prods))
            break

    return prods


def get_right_brdc(brdc_lis_inp, tmp_dir_inp):
    """
    Selects the appropriate BRDC (Broadcast Ephemeris) file from a list and unzips it.

    Parameters
    ----------
    brdc_lis_inp : list
        A list of BRDC file paths.
    tmp_dir_inp : str
        The directory where the unzipped BRDC file will be stored.

    Returns
    -------
    tuple
        A tuple containing the original BRDC file path and the unzipped BRDC file path.
    """
    if len(brdc_lis_inp) == 1:  ## normal case
        brdc_ori = brdc_lis_inp[0]
        brdc_unzip = files_rw.unzip_gz_Z(brdc_ori, out_gzip_dir=tmp_dir_inp)
    elif len(brdc_lis_inp) == 0:
        brdc_ori = None
        brdc_unzip = None
        log.warning("no brdc. found")
    elif len(brdc_lis_inp) > 1:
        log.warning("several brdc found, keep the last one")
        log.warning(brdc_lis_inp)
        brdc_ori = brdc_lis_inp[-1]
        brdc_unzip = files_rw.unzip_gz_Z(brdc_ori, out_gzip_dir=tmp_dir_inp)
    else:
        #### this should never happend
        brdc_ori = None
        brdc_unzip = None
        pass

    return brdc_ori, brdc_unzip


def get_best_latency(prod_lis_inp):
    """
    Selects the best latency from a list of product files.

    Internal function for get_right_prod.

    Parameters
    ----------
    prod_lis_inp : list
        A list of product file paths.

    Returns
    -------
    list
        The best latency found in the list.
        But can be several if the latency is the same for several products.

    """
    latency_priority_lis = ["FIN", "RAP", "ULT", "NRT"]

    out_prods_lis = []

    for lat in latency_priority_lis:
        for prod in prod_lis_inp:
            if lat in prod:
                out_prods_lis.append(prod)
                break

    if len(out_prods_lis) == 0:
        log.warning("no prod. found with a known latency")
        out_prods_lis = prod_lis_inp

    return out_prods_lis


def get_right_prod(prod_lis_inp, tmp_dir_inp, prod_name, default_fallback):
    """
    Selects the appropriate product file from a list and unzips it if necessary.

    Parameters
    ----------
    prod_lis_inp : list
        A list of product file paths.
    tmp_dir_inp : str
        The directory where the unzipped product file will be stored.
    prod_name : str
        The name of the product.
    default_fallback : bool
        If True, uses default values if products are not found.

    Returns
    -------
    tuple
        A tuple containing the unzipped product file path and the original product file path.
    """
    if len(prod_lis_inp) == 1:  ## normal case
        prod_ori = prod_lis_inp[0]
        prod_out = files_rw.unzip_gz_Z(prod_ori, out_gzip_dir=tmp_dir_inp)
    elif len(prod_lis_inp) == 0 and default_fallback:
        log.warning("no prod. %s found, fallback to 'Default' in cfg file", prod_name)
        prod_out = "Default"
        prod_ori = "Default"
    elif len(prod_lis_inp) == 0 and not default_fallback:
        log.warning(
            "no prod. %s found, no fallback to Default set (default_fallback option), aborting...",
            prod_name
        )
        prod_out = None
        prod_ori = None
    elif len(prod_lis_inp) > 1:
        log.warning("several prod found, search for the best latency")
        log.warning(prod_lis_inp)
        prod_lis_bst = get_best_latency(prod_lis_inp)

        if len(prod_lis_bst) == 1:
            prod_ori = prod_lis_bst[0]
        else:
            log.warning("several prod. found with the same latency, keep the last one")
            log.warning(prod_lis_bst)
            prod_ori = prod_lis_bst[-1]

        prod_out = files_rw.unzip_gz_Z(prod_ori, out_gzip_dir=tmp_dir_inp)
    else:
        #### this should never happend
        prod_out = None
        prod_ori = None
        pass

    return prod_out, prod_ori


def pride_pppar_runner_mono(
        rnx_path,
        cfg_template_path,
        prod_ac_name,
        prod_parent_dir,
        tmp_dir,
        cfg_dir,
        run_dir,
        cfg_prefix="pride_pppar_cfg_1a",
        mode="K",
        options_dic={},
        bin_dir=None,
        force=False,
        dl_prods=False,
        default_fallback=False,
        dl_prods_only=False,
        clean_run_dir=True
):
    """
    Runs the PRIDE PPPAR process for a single RINEX file.

    Parameters
    ----------
    rnx_path : str
        The path to the RINEX file.
    cfg_template_path : str
        The path to the configuration template file.
    prod_ac_name : str
        The name of the analysis center providing the products.
    prod_parent_dir : str
        The parent directory where the products are stored.
    tmp_dir : str
        The temporary directory for intermediate files.
    cfg_dir : str
        The directory for configuration files.
    run_dir : str
        The directory where the run results will be stored.
    cfg_prefix : str, optional
        The prefix for the configuration file name. Default is "pride_pppar_cfg_1a".
    mode : str, optional
        The mode for the PRIDE PPPAR process. Default is "K".
    options_dic : dict, optional
        Additional options for the PRIDE PPPAR process. Default is an empty dictionary.
    bin_dir : str, optional
        The directory where the PRIDE PPPAR binaries are located. Default is None.
    force : bool, optional
        If True, forces the process to run even if logs already exist. Default is False.
    dl_prods : bool, optional
        If True, downloads the necessary products. Default is False.
    default_fallback : bool, optional
        If True, uses default values if products are not found. Default is False.
    dl_prods_only : bool, optional
        If True, only downloads the products and exits. Default is False.
    clean_run_dir : bool, optional
        If True, removes temporary files inside the run directory. Default is True.

    Returns
    -------
    None
    """
    if not bin_dir:
        bin_dir = os.path.join(os.environ["HOME"], ".PRIDE_PPPAR_BIN")

    srt = conv.rinexname2dt(rnx_path)

    doy, year = conv.dt2doy_year(srt)
    rnx_file = os.path.basename(rnx_path)
    hourmin_str = str(srt.hour).zfill(2) + str(srt.minute).zfill(2)

    site = rnx_file[:4].upper()

    prod_ac_name_no_regex_char = remove_regex_reserved_characters(prod_ac_name)

    ########## DEFINE DIRECTORIES
    tmp_dir_use = os.path.join(tmp_dir, year, doy)  ### tmp is year/doy only, no site,
    ### because the 'common' dir must be the same
    cfg_dir_use = os.path.join(cfg_dir, year, doy, site)
    run_dir_use = os.path.join(
        run_dir, mode, prod_ac_name_no_regex_char, site
    )  ### pdp3 add year/doy by itself
    run_dir_ope = os.path.join(run_dir_use, year, doy)
    run_dir_fin = run_dir_ope + "_" + hourmin_str

    ########### CHECK IF LOGS ALREADY EXIST
    logs_existing = utils.find_recursive(run_dir_fin, "log*" + site.lower())
    if len(logs_existing) > 0 and not force:
        log.info("log exists for %s in %s, skip", rnx_file, run_dir_fin)
        return None
    else:
        pass

    utils.create_dir(tmp_dir_use)
    utils.create_dir(cfg_dir_use)
    utils.create_dir(run_dir_use)

    ########### DOWNLOAD PRODUCTS + BROADCASTS
    if dl_prods:
        _ = dl_prods_pride_pppar(prod_parent_dir, [srt], prod_ac_name)
        _ = dl_brdc_pride_pppar(prod_parent_dir, [srt])

    if dl_prods_only:
        log.info("products downloaded, exiting (dl_prods_only is activated.)")
        return None

    ########### UNCOMPRESS RINEX
    rnx_bnm = os.path.basename(rnx_path)
    if "crx" in rnx_bnm or "d.Z" in rnx_bnm or "d.gz" in rnx_bnm:
        rnx_path_tmp = shutil.copy(rnx_path, tmp_dir_use)
        rnx_path_use = str(hatanaka.decompress_on_disk(rnx_path_tmp, delete=True))
    else:
        rnx_path_use = rnx_path

    ########### MOVE BRDC IN TMP (so pdp3 can handle it)
    ## we must also rename the BKG brdc to fake it as the IGS one
    def _find_unzip_brdc():
        """
        Finds and unzips the BRDC file for the given day.

        Returns
        -------
        tuple
            The path to the unzipped BRDC file and the original BRDC file.
        """
        # 1ST STEP: find the brdc
        brdc_pattern = "*BRDC00WRD_S_" + year + doy + "*gz"
        prod_dir_year_doy = os.path.join(prod_parent_dir, year, doy)
        brdc_lis = utils.find_recursive(prod_dir_year_doy, brdc_pattern.strip())

        # 2ND STEP: unzip the brdc
        brdc_ori, brdc_unzip = get_right_brdc(brdc_lis, tmp_dir_use)

        ##### 3RD STEP: rename the brdc
        if not brdc_unzip:
            brdc_out = None
        else:
            # generate the final name for the brdc
            brdc_igs_nam = brdc_unzip.replace("BRDC00WRD_S_",
                                              "BRDC00IGS_R_")
            # best case: the final renamed brdc is already here
            if os.path.isfile(brdc_igs_nam):
                brdc_out = brdc_igs_nam
            # if the final renamed brdc is not here, we rename the unzipped one
            elif brdc_unzip and os.path.isfile(brdc_unzip):
                os.rename(brdc_unzip, brdc_igs_nam)
                brdc_out = brdc_igs_nam
            # a weird case where the unzipped brdc is missing. This should not happend
            else:
                brdc_out = brdc_unzip

        return brdc_out, brdc_ori

    ########### GENERATE CONFIG FILE
    def _find_unzip_prod(prod):
        """
        Finds and unzips the right products for a given day.

        Parameters
        ----------
        prod : str
            The type of product to find and unzip.

        Returns
        -------
        tuple
            The path to the unzipped product file and the original product file.
        """
        if "ULT" in prod_ac_name:
            add_hourly_file = True
        else:
            add_hourly_file = False

        find_prod_epoch_ini = srt

        smart_ultra = True
        if ("ULT" in prod_ac_name or "NRT" in prod_ac_name) and smart_ultra:
            delta_epoch_max = 24  # => [0 ... 23]
        else:
            delta_epoch_max = 1  # => [0]

        find_prods = operational.find_IGS_products_files

        prod_lis = []

        #### this loop is to find former ULTRA if the latest is missing
        # it works for ULTRA only, for other latencies, delta_epoch_max = 0
        for i_d_epo in range(delta_epoch_max):

            find_prod_epoch = find_prod_epoch_ini - dt.timedelta(seconds=3600 * i_d_epo)

            if smart_ultra and i_d_epo > 0:
                log.warning(
                    "no prods found for epoch %s, but we try %i hour before (%s)",
                    find_prod_epoch_ini,
                    i_d_epo,
                    find_prod_epoch,
                )

            prod_lis = find_prods(
                prod_parent_dir,
                [prod],
                [prod_ac_name],
                find_prod_epoch,
                severe=False,
                regex_old_naming=False,
                regex_igs_tfcc_naming=False,
                add_hourly_file=add_hourly_file,
            )
            if len(prod_lis) > 0:
                break

        #### differents behaviors depending on the number of files found
        prod_out, prod_ori = get_right_prod(prod_lis, tmp_dir_use, prod, default_fallback)

        return prod_out, prod_ori

    def _change_value_in_cfg(lines_file_inp, key_inp, val_inp, basename_on_val_inp=True):
        """
        Changes the values in the configuration file.

        Parameters
        ----------
        lines_file_inp : list of str
            The lines of the configuration file.
        key_inp : str
            The key to change.
        val_inp : str
            The new value for the key.
        basename_on_val_inp : bool, optional
            If True, uses the basename of the value. Default is True.
        """

        if not val_inp:
            log.warning("no value for %s, keep the Default value", key_inp)
            return None

        if basename_on_val_inp:
            val_use = os.path.basename(val_inp)
        else:
            val_use = val_inp

        for il, lin in enumerate(lines_file_inp):
            f = lin.split("=")
            if key_inp in f[0]:
                f[1] = val_use
                lines_file_inp[il] = "= ".join(f) + "\n"

    with open(cfg_template_path) as fil:
        cfg_lines = fil.readlines()

    ### find and unzip the right products
    sp3_path, _ = _find_unzip_prod("SP3")
    bia_path, _ = _find_unzip_prod("BIA")
    clk_path, _ = _find_unzip_prod("CLK")
    obx_path, _ = _find_unzip_prod("OBX")
    erp_path, _ = _find_unzip_prod("ERP")
    brdc_path, _ = _find_unzip_brdc()

    if (
            not any((sp3_path, bia_path, clk_path, obx_path, erp_path))
            and not default_fallback
    ):
        log.error(
            "an essential prod. at least is missing and no fallback to Default is set, abort"
        )
        return None

    ### change the values in the config file template
    _change_value_in_cfg(cfg_lines, "Product directory", tmp_dir_use, basename_on_val_inp=False)
    _change_value_in_cfg(cfg_lines, "Satellite orbit", sp3_path)
    _change_value_in_cfg(cfg_lines, "Satellite clock", clk_path)
    _change_value_in_cfg(cfg_lines, "ERP", erp_path)
    _change_value_in_cfg(cfg_lines, "Quaternions", obx_path)
    _change_value_in_cfg(cfg_lines, "Code/phase bias", bia_path)

    date_prod_midfix = year + doy + hourmin_str

    cfg_name = cfg_prefix + "_" + prod_ac_name_no_regex_char + "_" + date_prod_midfix
    cfg_path = os.path.join(cfg_dir_use, cfg_name)

    ### write the good config file
    with open(cfg_path, "w+") as fout:
        for l in cfg_lines:
            fout.write(l)

    ##### RUN THE STUFF
    os.environ["PATH"] += ":" + bin_dir

    cmd = " ".join(("pdp3", "--config", cfg_path, "--mode", mode, rnx_path_use))
    log.info(cmd)

    os.chdir(run_dir_use)  ## not run_dir_use, pdp3 will goes by itselft to yyyy/doy

    ## run the command ####
    run_command(cmd)
    #######################

    if clean_run_dir:
        run_dir_files = utils.find_recursive(run_dir_ope, "*")
        idel_files = 0
        for f in run_dir_files:
            if not re.match("[a-z]{3}_[0-9]{7}_.{4}", os.path.basename(f)):
                # log.info("removing %s", f)
                os.remove(f)
                idel_files += 1
        log.info("%s tmp files in run_dir removed", idel_files)

    # handle the cases where the run_dir_fin already exists
    if os.path.isdir(run_dir_fin):
        # the normal case where the run_dir_fin already exists, and we want to force the deletion
        if force:
            shutil.rmtree(run_dir_fin)
        # the unusual case where the run_dir_fin already exists, we did not want to force the deletion,
        # and then the process should have been skipped, but no log has been found inside but we keep it anyway
        else:
            log.warning("run_dir_fin %s already exists, but no log found inside", run_dir_fin)
            timstp = utils.get_timestamp()
            run_dir_fin_old = run_dir_fin + "_" + timstp
            log.warning("renaming the exisiting final dir as %s", os.path.basename(run_dir_fin))
            os.rename(run_dir_fin, run_dir_fin_old)

    ### FINAL rename the run_dir_fin to its final name run_dir_fin (with hourmin)
    os.rename(run_dir_ope, run_dir_fin)

    return None


def pride_pppar_mp_wrap(kwargs_inp):
    try:
        operational.pride_pppar_runner_mono(**kwargs_inp)
    except Exception as e:
        log.error("%s raised, RINEX is skipped: %s",
                  type(e).__name__,
                  kwargs_inp['rnx_path'])
        # pass
        # raise
    return None


def pride_pppar_runner(rnx_path_list,
                       cfg_template_path,
                       prod_ac_name,
                       prod_parent_dir,
                       tmp_dir,
                       cfg_dir,
                       run_dir,
                       multi_process=1,
                       cfg_prefix='pride_pppar_cfg_1a',
                       mode='K',
                       options_dic={},
                       bin_dir=None,
                       force=False,
                       dl_prods=False,
                       default_fallback=False,
                       dl_prods_only=False,
                       clean_run_dir=True):

    date_list = sorted(list(set([conv.rinexname2dt(rnx) - dt.timedelta(seconds=0) for rnx in rnx_path_list])))

    ## when we do multi_process, we download the products first, all at once, to avoid conflicts
    if dl_prods or dl_prods_only:
        if dl_prods_only and not dl_prods_only:
            log.warning("dl_prods_only is activated, but not dl_prods. Products download is forced anyway.")
        _ = dl_prods_pride_pppar(prod_parent_dir, date_list, prod_ac_name)
        _ = dl_brdc_pride_pppar(prod_parent_dir, date_list)
        if dl_prods_only:
            log.info("products downloaded, exiting (dl_prods_only is activated.)")
            return None

    kwargs_list = []
    for rnx_path in rnx_path_list:
        kwargs = {'rnx_path': rnx_path,
                  'cfg_template_path': cfg_template_path,
                  'prod_ac_name': prod_ac_name,
                  'prod_parent_dir': prod_parent_dir,
                  'tmp_dir': tmp_dir,
                  'cfg_dir': cfg_dir,
                  'run_dir': run_dir,
                  'cfg_prefix': cfg_prefix,
                  'bin_dir': bin_dir,
                  'mode': mode,
                  'options_dic': options_dic,
                  'force': force,
                  'dl_prods': dl_prods,
                  'default_fallback': default_fallback,
                  'clean_run_dir': clean_run_dir}

        kwargs_list.append(kwargs)

    if multi_process > 1:
        log.info("multiprocessing: %d cores used", multi_process)

    pool = mp.Pool(processes=multi_process)
    results_raw = [pool.apply_async(pride_pppar_mp_wrap, args=(x,)) for x in kwargs_list]
    results = [e.get() for e in results_raw]

    return results
