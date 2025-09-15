#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 28/02/2025 18:52:44

@author: psakic
"""
import re
import subprocess

import geodezyx.operational.gins_runner.gins_common as gynscmn
import geodezyx.operational.gins_runner.gins_dirs_gen as gynsgen
import geodezyx.operational.gins_runner.gins_dirs_run as gynsrun
import geodezyx.operational.gins_runner.gins_bd_update as gynsbdu
import geodezyx.operational.gins_runner.gins_orbclk_concat as gynsorbclkcat
import datetime as dt

from geodezyx import utils, operational

import multiprocessing as mp
import os
import shutil
import glob
import time

from geodezyx import files_rw, conv

spotgins_wrap = None
spotgins_wrap2 = None

#### Import the logger
import logging

log = logging.getLogger("geodezyx")


LAST_VERSION_VALIDE = "VALIDE_25_1"
#DIR_SPOTGINS_DEFAULT = ""DIR_SPOTGINS_G20_GE.yml""
DIR_SPOTGINS_DEFAULT = "DIR_SPOTGINS_G20_GE_VALIDE_25_1.yml"

def spotgins_run(
    rnxs_path_inp,
    results_folder_inp,
    nprocs=8,
    version=LAST_VERSION_VALIDE,
    const="GE",
    director_generik_path_inp=None,
    director_name_prefix_inp="",
    stations_file_inp=None,
    oceanload_file_inp=None,
    options_prairie_file_inp=None,
    stations_master_file_inp=None,
    force=False,
    no_archive=False,
    no_clean_tmp=False,
    no_updatebd=False,
    no_concat_orb_clk=False,
    verbose=False,
    updatebd_login="",
):
    """
    Run the SPOTGINS process on the provided RINEX paths using multiprocessing.

    Parameters
    ----------
    rnxs_path_inp : list
        List of input RINEX file paths.
    results_folder_inp : str
        Path to the archive folder.
    nprocs : int, optional
        Number of processes to use for multiprocessing. Default is 8.
    version : str, optional
        Version of the GINS software to use.
         Default is the global variable LAST_VERSION_VALIDE.
    const : str, optional
        GNSS constellation to use. Default is "G".
    director_generik_path_inp : str, optional
        Path to the generic director file. Default is None.
    director_name_prefix_inp : str, optional
        Prefix for the director name. Default is an empty string.
    stations_file_inp : str, optional
        Path to the stations file. Default is None.
    oceanload_file_inp : str, optional
        Path to the ocean load file. Default is None.
    options_prairie_file_inp : str, optional
        Path to the options prairie file. Default is None.
    stations_master_file_inp : str, optional
        Path to the stations master file. Default is None.
    force : bool, optional
        Force the generation and run of directors even if they already exist. Default is False.
    no_archive : bool, optional
        If True, do not archive the results. Default is False.
    no_clean_tmp : bool, optional
        If True, do not clean the temporary files. Default is False.
    no_updatebd : bool, optional
        If True, do not update the BDGINS repository. Default is False.
    no_concat_orb_clk : bool, optional
        If True, do not concatenate the orbit and clock files prior to the main processing.
         Default is False.
    verbose : bool, optional
        If True, enable verbose logging. Default is False.
    updatebd_login : str
        The login to connect to the remote `tite` GINS server to update the database.
        We assume that SSH keys have been exchanged to automatize the connexion
    Returns
    -------
    None
    """

    ##### sort the input list
    rnxs_path_use = utils.listify(rnxs_path_inp)
    rnxs_path_use = operational.rinex_table_from_list(
        rnxs_path_use, site9_col=True, size_col=False
    )
    rnxs_path_use = rnxs_path_use["path"].to_list()

    if len(rnxs_path_use) == 0:
        log.error("No input RINEX files to process, skip")
        return

    ##### get the paths of the files needed for SPOTGINS
    dirgen_use, stfi_use, oclo_use, opra_use, siteid9_use = get_spotgins_files(
        director_generik_path_inp,
        stations_file_inp,
        oceanload_file_inp,
        options_prairie_file_inp,
        stations_master_file_inp,
    )

    ##### Need timespan
    if not no_updatebd or not no_concat_orb_clk:
        rnxs_dates = [conv.rinexname2dt(e) for e in rnxs_path_use]
        rnxs_dates = [e for e in rnxs_dates if not e is None]
        date_min, date_max = (
            min(rnxs_dates) - dt.timedelta(days=1),
            max(rnxs_dates) + dt.timedelta(days=1),
        )
    else:
        date_min = dt.datetime(2000, 5, 3)
        date_max = None
    ##### Update the database ################
    if not no_updatebd:
        gynsbdu.bdgins_update(
            date_srt=date_min, date_end=date_max, dir_bdgins="", login=updatebd_login
        )

    ##### concatenate hor/orb ################
    if not no_concat_orb_clk:
        gynsorbclkcat.concat_orb_clk(date_min, date_max, nprocs=nprocs)

    ##### Multi-processing Wrapper ################
    global spotgins_wrap

    def spotgins_wrap(rnx_mono_path_inp):
        ######## QUICK ARCHIVING CHECK ############# # Check if the solution is already archived
        if not force and check_arch_sol(
            rnx_mono_path_inp, results_folder_inp, verbose=verbose
        ):
            log.info(f"Solution already archived for {rnx_mono_path_inp}, skip")
            return None

        ######## DIRECTORS GENERATION ########
        dirr, tmp_folder = gynsgen.gen_dirs_rnxs(
            rnx_paths_inp=rnx_mono_path_inp,
            director_generik_path=dirgen_use,
            director_name_prefix="SPOTGINS",
            stations_file=stfi_use,
            oceanload_file=oclo_use,
            auto_stations_file=False,
            auto_oceanload=False,
            options_prairie_file=opra_use,
            auto_interval=False,
            force=force,
            verbose=verbose,
            sites_id9=siteid9_use,
        )

        ######## DIRECTORS RUN ###############
        const_use = const_adapt(const, dirr, verbose=verbose)
        # const_use = "GE" # ASG use always GE
        opt_gins_90_use = "-const " + const_use
        try:
            gynsrun.run_directors(
                dirr,
                opts_gins_90=opt_gins_90_use,
                version=version,
                cmd_mode="exe_gins_dir",
                force=force,
                verbose=verbose,
                sleep_time_max=nprocs * 10**-1,  # eq 128 procs -> 12.8s max sleep time
            )
        except Exception as e:
            log.error(f"run_directors error: {e}")
            rm_prov_listing(dirr)

        ######## ARCHIVING ####################
        if not no_archive:
            archiv_gins_run(dirr, results_folder_inp, verbose=verbose)

        ######## CLEANING ####################
        if not no_clean_tmp and os.path.exists(tmp_folder):
            if verbose:
                log.info(f"Cleaning temporary folder {tmp_folder}")
            shutil.rmtree(tmp_folder)
        return None

    ##### END Multi-processing Wrapper END ################

    global spotgins_wrap2

    def spotgins_wrap2(*args):
        try:
            return spotgins_wrap(*args)
        except Exception as e:
            log.error(f"spotgins_wrap error: {e}")
            return None  # or something that won't break your pipeline

    pool = mp.Pool(processes=nprocs)
    try:
        res_raw = pool.map(spotgins_wrap2, rnxs_path_use, chunksize=1)
    except Exception as e:
        log.error("pool.map error: %s", e)
    pool.close()

    nrnxs = len(rnxs_path_use)
    log.info("multi run end, all of the %i RINEXs have been procssed", nrnxs)

    return


def get_spotgins_files(
    director_generik_path_inp,
    stations_file_inp,
    oceanload_file_inp,
    options_prairie_file_inp,
    stations_master_file_inp,
):
    """
    Retrieve the paths of the files needed for the SPOTGINS process.

    Parameters
    ----------
    director_generik_path_inp : str
        Path to the generic director file.
    stations_file_inp : str
        Path to the stations file.
    oceanload_file_inp : str
        Path to the ocean load file.
    options_prairie_file_inp : str
        Path to the options prairie file.
    stations_master_file_inp : str
        Path to the stations master file.

    Returns
    -------
    tuple
        A tuple containing the paths to the director file,
        stations file, ocean load file, options prairie file,
        and site IDs.
    """
    if not gynscmn.get_spotgins_path() and (
        not stations_file_inp or not oceanload_file_inp or not options_prairie_file_inp
    ):
        raise ValueError(
            "SPOTGINS path not set ($SPOTGINS_DIR) and stations_file/oceanload_file/options_prairie_file not provided"
        )
    else:
        sptgns_path = gynscmn.get_spotgins_path()

    if director_generik_path_inp:
        dirgen_use = director_generik_path_inp
    else:
        dirgen_use = os.path.join(
            sptgns_path, "metadata", "directeur", DIR_SPOTGINS_DEFAULT
        )

    if stations_file_inp:
        stfi_use = stations_file_inp
    else:
        stfi_use = os.path.join(sptgns_path, "metadata", "stations", "station_file.dat")

    if oceanload_file_inp:
        oclo_use = oceanload_file_inp
    else:
        oclo_use = os.path.join(
            sptgns_path, "metadata", "oceanloading", "load_fes2014b_cf.spotgins"
        )

    if options_prairie_file_inp:
        opra_use = options_prairie_file_inp
    else:
        opra_use = os.path.join(
            sptgns_path, "metadata", "directeur", "options_prairie_static"
        )

    if stations_master_file_inp and os.path.isfile(stations_master_file_inp):
        siteid9_use = files_rw.read_spotgins_masterfile(stations_master_file_inp)
        siteid9_use = siteid9_use["NAME"]
    elif sptgns_path:
        stations_master_file = os.path.join(
            sptgns_path, "metadata", "stations", "station_master_file.dat"
        )
        siteid9_use = files_rw.read_spotgins_masterfile(stations_master_file)
        siteid9_use = siteid9_use["NAME"]
    else:
        siteid9_use = None

    return dirgen_use, stfi_use, oclo_use, opra_use, siteid9_use


def rm_prov_listing(dir_inp):
    """
    Remove potential PROV folders and related files in the listing folder.

    If a processing is interrupted this temp PROV folder remains,
    occupy a lot of space and can eventually crashs the whole machine/container.

    Thus it must be removed ASAP


    Parameters
    ----------
    dir_inp : str
        Path to the input directeur whose associated PROV folders and files
        need to be removed.

    Functionality
    -------------
    1. Identifies and removes all PROV folders associated with the given directeur.
    2. Identifies and removes specific files (e.g., `.qsub` and `.sh` files)
       related to the directory in the batch listing folder.

    Returns
    -------
    None
    """
    # Get the path to the batch listing folder
    li_batch_fld = os.path.join(gynscmn.get_gin_path(True), "batch", "listing")
    # Extract the base name of the input directory
    dir_basename = os.path.basename(dir_inp)

    # Find and remove all PROV folders associated with the directory
    prov_lis = glob.glob(os.path.join(li_batch_fld, "PROV" + "*" + dir_basename + "*"))
    for prov in prov_lis:
        log.warning("removing PROV folder %s", prov)
        shutil.rmtree(prov)

    # Find and remove specific files (e.g., `.qsub` and `.sh`) in the batch listing folder
    lifil_lis = glob.glob(os.path.join(li_batch_fld, dir_basename) + "*")
    for f in lifil_lis:
        if f.endswith(".qsub") or f.endswith("sh"):
            os.remove(f)

    return None


def archiv_gins_run(dir_inp, archive_folder, verbose=True):
    """
    Archive the GINS run results to the specified archive folder.

    Parameters
    ----------
    dir_inp : str
        Path to the directory containing the GINS run results.
    archive_folder : str
        Path to the archive folder.
    verbose : bool, optional
        If True, enable verbose logging. Default is True.

    Returns
    -------
    None
    """

    if not dir_inp:
        return None

    time.sleep(1)  # wait a proper GINS end

    # do first the PROV folder cleaning
    rm_prov_listing(dir_inp)

    dir_basename = os.path.basename(dir_inp)
    site_id9 = re.search(r"_(....00\w{3})_", dir_inp).group(1)
    arch_fld_site = str(os.path.join(archive_folder, site_id9))
    if not os.path.exists(arch_fld_site):
        os.makedirs(arch_fld_site)

    # get input directories
    dir_batch_fld = os.path.join(gynscmn.get_gin_path(True), "data", "directeur")
    li_batch_fld = os.path.join(gynscmn.get_gin_path(True), "batch", "listing")
    sol_batch_fld = os.path.join(gynscmn.get_gin_path(True), "batch", "solution")
    stat_batch_fld = os.path.join(gynscmn.get_gin_path(True), "batch", "statistiques")

    # get destination
    dir_arch_fld = os.path.join(arch_fld_site, "010_directeurs")
    li_arch_fld = os.path.join(arch_fld_site, "020_listings")
    sol_arch_fld = os.path.join(arch_fld_site, "030_solutions")
    stat_arch_fld = os.path.join(arch_fld_site, "040_statistiques")
    trsh_arch_fld = os.path.join(arch_fld_site, "090_trash")

    # create directories
    for arch_fld in [
        dir_arch_fld,
        li_arch_fld,
        sol_arch_fld,
        stat_arch_fld,
        trsh_arch_fld,
    ]:
        if not os.path.exists(arch_fld):
            os.makedirs(arch_fld)

    log.info(f"archive {dir_basename} run in {arch_fld_site}")
    # move files to their destination
    for batch_fld, arch_fld in zip(
        [dir_batch_fld, li_batch_fld, sol_batch_fld, stat_batch_fld],
        [dir_arch_fld, li_arch_fld, sol_arch_fld, stat_arch_fld],
    ):

        for f in glob.glob(os.path.join(batch_fld, dir_basename) + "*"):
            if f.endswith(".qsub") or f.endswith("sh"):
                os.remove(f)
                continue

                # if verbose:
                #    log.info(f"skip archive of qsub/sh file {f}")
                # continue

            if verbose:
                log.info(f"Archiving {f} to {arch_fld}")

            # check if the file is already in the archive (crash if yes)
            # and move it to the trash folder
            f_bn = os.path.basename(f)
            f_prex = os.path.join(arch_fld, f_bn)
            if os.path.exists(f_prex):
                timstp = utils.get_timestamp()
                f_prex_trsh = os.path.join(trsh_arch_fld, f_bn + "_" + timstp)
                shutil.move(f_prex, f_prex_trsh)
            # move the actual correct file
            shutil.move(f, arch_fld)

    # compress the listings
    for f in glob.glob(os.path.join(li_arch_fld, dir_basename) + "*"):
        if not f.endswith("gz"):
            if verbose:
                log.info(f"Compressing listing {os.path.basename(f)}")
            utils.gzip_compress(f, rm_inp=True)
    return None


def check_arch_sol(rnx_path_inp, archive_folder_inp, verbose=True):
    """
    Check if the solution for the given RINEX file is already archived.

    Parameters
    ----------
    rnx_path_inp : str
        Path to the input RINEX file.
    archive_folder_inp : str
        Path to the archive folder.
    verbose : bool, optional
        If True, enable verbose logging. Default is True.

    Returns
    -------
    bool
        True if the solution is found in the archive, False otherwise.
    """
    epo = conv.rinexname2dt(rnx_path_inp)
    epo_str = epo.strftime("%Y_%j")
    site_id4 = str(os.path.basename(rnx_path_inp)[0:4].upper())

    potential_sol = glob.glob(
        os.path.join(
            archive_folder_inp,
            "*" + site_id4 + "*",
            "030_solutions",
            "*" + site_id4 + "*" + epo_str + "*",
        )
    )

    if len(potential_sol) > 0:
        if verbose:
            log.info(f"Solution(s) found: {potential_sol}")
        return True
    else:
        return False


def const_adapt(const_inp, dir_inp, verbose=True):
    """
    manage
    ERREUR : Galileo integer ambiguity fixing not possible before 7th October 2018.
    solution
    GPS-only before this date
    """

    import yaml

    dir_dic = yaml.load(open(dir_inp), Loader=yaml.FullLoader)
    rnx_date = conv.jjul_cnes2dt(dir_dic["date"]["arc_start"][0])
    if "E" in const_inp and rnx_date < dt.datetime(2018, 10, 7):
        if verbose:
            msg = "Galileo not recommended before 2018-10-07, switch for GPS-only"
            log.warning(msg)
        return "G"
    else:
        return const_inp
