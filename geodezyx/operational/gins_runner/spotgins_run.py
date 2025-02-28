#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 28/02/2025 18:52:44

@author: psakic
"""

import geodezyx.operational.gins_runner.gins_common as gynscmn
import geodezyx.operational.gins_runner.gins_dirs_gen as gynsgen
import geodezyx.operational.gins_runner.gins_dirs_run as gynsrun

import multiprocessing as mp
import os

def spotgins_run(
    rnxs_path_inp,
    nprocs=8,
    version="VALIDE_23_1",
    const="G",
    director_name_prefix_inp="",
    director_generik_path_inp=None,
    stations_file_inp=None,
    oceanload_file_inp=None,
    options_prairie_file_inp=None,
    force_gen=False,
    force_run=False,
):

    if not gynscmn.get_spotgins_path() and (
        not stations_file_inp or not oceanload_file_inp or not options_prairie_file_inp
    ):
        raise ValueError(
            "SPOTGINS path not set ($SPOTGINS_DIR) and stations_file, oceanload_file or options_prairie_file not provided"
        )
    else:
        sptgns_path = gynscmn.get_spotgins_path()

    stfi_use = stations_file_inp or os.path.join(
        sptgns_path, "metadata", "directeur", "DIR_SPOTGINS_G20_GE.yml"
    )
    oclo_use = oceanload_file_inp or os.path.join(
        sptgns_path, "metadata", "oceanloading", "load_fes2014b_cf.spotgins"
    )
    opra_use = options_prairie_file_inp or os.path.join(
        sptgns_path, "metadata", "directeur", "options_prairie_static"
    )

    ##### Multi-processing Wrapper ################
    def _spotgins_wrap(rnx_mono_path_inp):

        director_generik_path_use = director_generik_path_inp
        dirr = gynsgen.gen_dirs_rnxs(
            rnx_paths_inp=rnx_mono_path_inp,
            director_generik_path=director_generik_path_use,
            director_name_prefix="_".join(("SPOTGINS", director_name_prefix_inp)),
            stations_file=stfi_use,
            oceanload_file=oclo_use,
            options_prairie_file=opra_use,
            auto_interval=False,
            force=force_gen,
        )

        opt_gins_90_use = "-const " + const
        gynsrun.run_directors(
            dirr,
            opts_gins_90=opt_gins_90_use,
            version=version,
            cmd_mode="exe_gins_dir",
            force=force_run,
        )

    ##### END Multi-processing Wrapper END ################

    pool = mp.Pool(processes=nprocs)
    res_raw = pool.map(_spotgins_wrap, rnxs_path_inp, chunksize=1)
