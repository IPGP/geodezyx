#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 28/02/2025 18:52:44

@author: psakic
"""
import re

import geodezyx.operational.gins_runner.gins_common as gynscmn
import geodezyx.operational.gins_runner.gins_dirs_gen as gynsgen
import geodezyx.operational.gins_runner.gins_dirs_run as gynsrun

import multiprocessing as mp
import os
import shutil
import glob

from geodezyx import files_rw

spotgins_wrap = None


def spotgins_run(
    rnxs_path_inp,
    nprocs=8,
    archive_folder_inp=None,
    version="VALIDE_23_2",
    const="G",
    director_generik_path_inp=None,
    director_name_prefix_inp="",
    stations_file_inp=None,
    oceanload_file_inp=None,
    options_prairie_file_inp=None,
    stations_master_file_inp=None,
    force_gen=False,
    force_run=False,
    no_archive=False,
    no_clean_tmp=False,
):

    rnxs_path_use = list(sorted(rnxs_path_inp))
    ##### get the paths of the files needed for SPOTGINS
    dirgen_use, stfi_use, oclo_use, opra_use, siteid9_use = get_spotgins_files(
        director_generik_path_inp,
        stations_file_inp,
        oceanload_file_inp,
        options_prairie_file_inp,
        stations_master_file_inp,
    )

    ##### Multi-processing Wrapper ################
    global spotgins_wrap

    def spotgins_wrap(rnx_mono_path_inp):
        ######## DIRECTORS GENERATION ########
        dirr, tmp_folder = gynsgen.gen_dirs_rnxs(
            rnx_paths_inp=rnx_mono_path_inp,
            director_generik_path=dirgen_use,
            director_name_prefix="SPOTGINS",
            stations_file=stfi_use,
            oceanload_file=oclo_use,
            options_prairie_file=opra_use,
            auto_interval=False,
            force=force_gen,
            sites_id9_series=siteid9_use,
        )

        ######## DIRECTORS RUN ###############
        opt_gins_90_use = "-const " + const
        gynsrun.run_directors(
            dirr,
            opts_gins_90=opt_gins_90_use,
            version=version,
            cmd_mode="exe_gins_dir",
            force=force_run,
        )

        ######## ARCHIVING ####################
        if not no_archive and archive_folder_inp:
            archive_gins_run(dirr, archive_folder_inp)

        ######## CLEANING ####################
        if not no_clean_tmp and os.path.exists(tmp_folder):
            shutil.rmtree(tmp_folder)

    ##### END Multi-processing Wrapper END ################

    pool = mp.Pool(processes=nprocs)
    res_raw = pool.map(spotgins_wrap, rnxs_path_use, chunksize=1)
    pool.close()

    return


def get_spotgins_files(
    director_generik_path_inp,
    stations_file_inp,
    oceanload_file_inp,
    options_prairie_file_inp,
    stations_master_file_inp,
):
    if not gynscmn.get_spotgins_path() and (
        not stations_file_inp or not oceanload_file_inp or not options_prairie_file_inp
    ):
        raise ValueError(
            "SPOTGINS path not set ($SPOTGINS_DIR) and stations_file, oceanload_file or options_prairie_file not provided"
        )
    else:
        sptgns_path = gynscmn.get_spotgins_path()

    if director_generik_path_inp:
        dirgen_use = director_generik_path_inp
    else:
        dirgen_use = os.path.join(
            sptgns_path, "metadata", "directeur", "DIR_SPOTGINS_G20_GE.yml"
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


def archive_gins_run(dir_inp, archive_folder):
    dir_basename = os.path.basename(dir_inp)
    site_id9 = re.match(r"....00\w{3}", dir_basename).group(1)
    arch_fld_site = str(os.path.join(archive_folder, site_id9))
    if not os.path.exists(arch_fld_site):
        os.makedirs(arch_fld_site)

    dir_batch_fld = os.path.join(gynscmn.get_gin_path(True), "data", "directeur")
    li_batch_fld = os.path.join(gynscmn.get_gin_path(True), "batch", "listing")
    sol_batch_fld = os.path.join(gynscmn.get_gin_path(True), "batch", "solution")

    dir_arch_fld = os.path.join(arch_fld_site, "010_directeurs")
    li_arch_fld = os.path.join(arch_fld_site, "020_listings")
    sol_arch_fld = os.path.join(arch_fld_site, "030_solutions")

    for arch_fld in [dir_arch_fld, li_arch_fld, sol_arch_fld]:
        if not os.path.exists(arch_fld):
            os.makedirs(arch_fld)

    for batch_fld, arch_fld in zip([dir_batch_fld, li_batch_fld, sol_batch_fld],
                                   [dir_arch_fld, li_arch_fld, sol_arch_fld]):

        for f in glob.glob(os.path.join(batch_fld, "*" + dir_basename + "*")):
            shutil.move(f, arch_fld)

    return
