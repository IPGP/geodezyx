#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26/02/2025 11:06:06

@author: psakic
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
import os
import time

import numpy as np
import yaml

#### geodeZYX modules
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import operational
from geodezyx import utils

#### Import geodezyx GINS submodules
import geodezyx.operational.gins_runner.gins_common as gzgicmn
import geodezyx.operational.gins_runner.gins_prairie as gzgipra
import geodezyx.operational.gins_runner.gins_orbclk as gzgiorb

#### Import the logger
import logging

log = logging.getLogger("geodezyx")


def gen_dirs_from_rnxs(
    rnx_paths_inp,
    director_generik_path,
    director_name_prefix,
    out_director_folder=None,
    temp_data_folder=None,
    stations_file=None,
    oceanload_file=None,
    auto_staocl=False,
    perso_orbclk=False,
    ac="igs",
    repro=2,
    auto_interval=True,
    out_coords="NULL",
    prairie=False,
    prairie_kwargs={"with_historik": 1, "with_wsb": 1},
):
    """
    Generate directors from RINEX files.

    This function can manage both a single RINEX file or a list of RINEX files.
    It produces a director from a RINEX file and writes it to a specified location.

    If the RINEX is not in a gins style folder, it will be copied in the ``temp_data_folder``

    Parameters
    ----------
    rnx_paths_inp : str or list
        Path to a RINEX file or a list of RINEX file paths.
    director_generik_path : str
        Path to the generic director file.
    director_name_prefix : str
        Prefix for the director name.
    out_director_folder : str, optional
        Output folder for the director. Defaults to None.
        is not specified, output director will be
        created in the ``../gin/data/directeur`` folder
    temp_data_folder : str, optional
        Temporary data folder. Defaults to None.
        is not specified RINEXs will be copied in
        an ad hoc folder ``../gin/TEMP_DATA``
    stations_file : str, optional
        Path to the stations file. Defaults to None.
    oceanload_file : str, optional
        Path to the ocean load file. Defaults to None.
    auto_staocl : bool, optional
        Automatically create station and ocean loading files.
        create automatically a station file and a ocean loading file
        with the rinex header
        auto mode is prioritary upon the manu mode
        so, if a path for stat file or ocload file is specified and
        auto_staocl is on, it will be the automatic stat/ocload
        which will be used
        Defaults to False.
    perso_orbclk : bool, optional
        Download and use specific orbits/clocks. Defaults to False.
    ac : str, optional
        Analysis center for perso orbits/clocks. Defaults to "igs".
    repro : int, optional
        Reprocessing version for perso orbits/clocks. Defaults to 2.
    auto_interval : bool, optional
        Automatically find the interval in RINEX and apply it to the director. Defaults to True.
    out_coords : str, optional
        Output coordinates type ('XYZ' or 'FLH'/'PLH'). Defaults to "NULL".
    prairie : bool, optional
        Run prairie externally. Defaults to False.
        if True, prairie_kwargs control the arguments of the function
        prairie_manual (cf above)
    prairie_kwargs : dict, optional
        Arguments for the prairie_manual function. Defaults to {"with_historik": 1, "with_wsb": 1}.

    Returns
    -------
    str or list
        Path of the produced director as a string, or a list of paths if multiple RINEX files are processed.
    """

    # Multi or Single Mode ?
    if type(rnx_paths_inp) is list:
        multimode = True
        rnx_path_lis = rnx_paths_inp
        log.info("******** DIRECTORS GENRATION ********")
    elif type(rnx_paths_inp) is str:
        multimode = False
        rnx_path_lis = [rnx_paths_inp]
    else:
        log.error("check the rinex_paths_in !!!")
        return None

    N = len(rnx_path_lis)
    director_output_path_lis = []
    failed_rinex_date_lis = []

    def _fail_rnx(rnx_path_inp, rnx_dt_inp, message):
        log.error(message + " failed, skip %s", rnx_path_inp)
        failed_rinex_date_lis.append(rnx_dt_inp)
        bool_continue = False
        if multimode:
            bool_continue = True
            return bool_continue
        else:
            raise Exception

    for i, rnx_path in enumerate(rnx_path_lis):
        if multimode:
            log.info(" ======== ", i + 1, "/", N, " ======== ")
            log.info(" === ", os.path.basename(rnx_path))
        rnx_name = os.path.basename(rnx_path)
        rnx_dt = conv.rinexname2dt(rnx_name)
        gins_path = gzgicmn.get_gin_path()

        # be sure there is a TEMP DATA folder
        if not temp_data_folder:
            temp_data_folder = os.path.join(gins_path, "gin", "TEMP_DATA")
        if not os.path.exists(temp_data_folder):
            os.makedirs(temp_data_folder)

        # be sure the RINEX is in gins folder ...
        bool_rnx_in_gin = gzgicmn.check_if_in_gin(rnx_path)
        # ... and copy it otherwise
        if not bool_rnx_in_gin:
            log.info("INFO : will be copied in %s", temp_data_folder)
            rnx_path = gzgicmn.copy_in_gin(rnx_path, temp_data_folder)
        # check if the RINEX is compressed ...
        bool_comp_rnx = operational.check_if_compressed_rinex(rnx_path)
        # ... if not crz2rnx !
        if bool_comp_rnx:
            crinex_path = rnx_path
            rnx_path = operational.crz2rnx(crinex_path, temp_data_folder)
            if not os.path.isfile(rnx_path):
                bool_cntu = _fail_rnx(rnx_path, rnx_dt, "CRZ2RNX")
                if bool_cntu:
                    continue

        # prairie extern
        if prairie:
            log.info("run prairie externally")
            pra_file_path = gzgipra.prairie_manual(
                rnx_path, temp_data_folder, **prairie_kwargs
            )

            if type(pra_file_path) is list and len(pra_file_path) == 0:
                bool_cntu = _fail_rnx(rnx_path, rnx_dt, "manual prairie")
                if bool_cntu:
                    continue

        if not out_director_folder:
            out_director_folder = os.path.join(
                gzgicmn.get_gin_path(), "gin", "data", "directeur"
            )

        stat_lower = rnx_name[0:4]
        stat = stat_lower.upper()

        if perso_orbclk:
            ac_suffix = "_" + ac + str(repro)
        else:
            ac_suffix = ""

        # date
        try:
            srt_epo, end_epo, freq_rnx = operational.rinex_start_end(rnx_path, True)
        except:
            bool_cntu = _fail_rnx(rnx_path, rnx_dt, "get RINEX start/end")
            if bool_cntu:
                continue

        ses = "_" + operational.rinex_session_id(srt_epo, end_epo, full_mode=True)

        # Interval Initalisation
        # interval_line = utils.grep(rnx_path,'INTERVAL',True)
        #  float(interval_line.split()[0])
        if auto_interval:
            freq_rnx_str = "_" + str(int(freq_rnx)).zfill(2) + "s"
        else:
            freq_rnx_str = ""

        # coord type suffix
        if out_coords.upper() in ("FLH", "XYZ"):
            coord_prefix = "_" + out_coords.lower()
        else:
            coord_prefix = ""

        dir_out_fname = utils.join_improved(
            "_",
            director_name_prefix,
            stat_lower,
            str(conv.dt2jjul_cnes(rnx_dt)),
            str(rnx_dt.year),
            conv.dt2doy(rnx_dt),
            freq_rnx_str,
            coord_prefix,
            ses,
            ac_suffix,
            ".yml",
        )

        dir_out_path = os.path.join(out_director_folder, dir_out_fname)
        director_generik_file = open(director_generik_path)
        dir_dic = yaml.load(director_generik_file, Loader=yaml.FullLoader)

        # ADDING RINEX NAME & DATES INTO THE DIR
        # name
        if not prairie:
            rnx_path_gstyl = gzgicmn.make_path_ginsstyle(rnx_path)
            dir_dic["observation"]["interobject_data"][1]["file"] = rnx_path_gstyl
        else:
            pra_path_dir_gstyl = gzgicmn.make_path_ginsstyle(pra_file_path)
            dir_dic["observation"]["interobject_data"][1]["file"] = pra_path_dir_gstyl

        # date
        strt_day, strt_sec = conv.dt2jjul_cnes(srt_epo, False)
        end_day, end_sec = conv.dt2jjul_cnes(end_epo, False)

        dir_dic["date"]["arc_start"][0] = strt_day
        dir_dic["date"]["arc_stop"][0] = end_day
        if "initial_state_vector_date" in list(dir_dic["date"].keys()):
            dir_dic["date"]["initial_state_vector_date"][0] = strt_day

        dir_dic["date"]["arc_start"][1] = strt_sec
        dir_dic["date"]["arc_stop"][1] = end_sec
        if "initial_state_vector_date" in list(dir_dic["date"].keys()):
            dir_dic["date"]["initial_state_vector_date"][1] = strt_sec

        # interval
        if auto_interval:  # dir is modified only if freq_rnx >= 0
            dir_dic["observation"]["removal"]["simulation_stepsize"] = freq_rnx
            # check if its not a static dir
            if (
                dir_dic["parameter"]["adjustment_parameters"]["stations"][
                    "adjustment_frequency"
                ][-1]
                != 0
            ):
                dir_dic["parameter"]["adjustment_parameters"]["stations"][
                    "adjustment_frequency"
                ] = [0, 0, 0, freq_rnx]
        log.info(
            "INFO : auto_interval: %s, simulation_stepsize : %s",
            auto_interval,
            dir_dic["observation"]["removal"]["simulation_stepsize"],
        )
        # output coords
        if out_coords.upper() in ("FLH", "PLH"):
            dir_dic["parameter"]["adjustment_parameters"]["stations"][
                "coordinates"
            ] = "latitude_longitude_height"
        elif out_coords.upper() == "XYZ":
            dir_dic["parameter"]["adjustment_parameters"]["stations"][
                "coordinates"
            ] = "cartesian_xyz"

        # station file and oceanload file , auto or manu mode
        # auto mode is prioritary upon the manu mode

        if not auto_staocl and (not oceanload_file or not stations_file):
            log.warning("auto_staocl is off and no stat/oclo file specified !!!")
        if auto_staocl:
            stations_file, oceanload_file = _dir_auto_staocl(
                stat_lower, rnx_dt, rnx_path, temp_data_folder
            )

        if oceanload_file:
            oceanload_file_ingin = gzgicmn.bring_to_gin(
                oceanload_file, temp_data_folder
            )
            dir_dic["object"]["station"]["ocean_tide_loading"] = (
                gzgicmn.make_path_ginsstyle(oceanload_file_ingin)
            )
        if stations_file:
            stations_file_ingin = gzgicmn.bring_to_gin(
                stations_file, temp_data_folder
            )
            dir_dic["object"]["station"]["station_coordinates"] = (
                gzgicmn.make_path_ginsstyle(stations_file_ingin)
            )

        # =========   ORBITS   =============
        if not perso_orbclk:
            orbpath, horpath = _dir_regular_orbclk(rnx_dt)
        else:
            orbpath, horpath = gzgiorb.download_convert_2_gins_orb_clk(
                rnx_dt, temp_data_folder, ac=ac, repro=repro
            )

        # Exception case where no sp3/clk was found
        if orbpath is None or horpath is None:
            bool_cntu = _fail_rnx(rnx_path, rnx_dt, "orbits/clocks download/conversion")
            if bool_cntu:
                continue

            orbpath = gzgicmn.make_path_ginsstyle(orbpath)
            horpath = gzgicmn.make_path_ginsstyle(horpath)

        dir_dic["model"]["environment"]["gnss_clock"] = horpath
        dir_dic["observation"]["interobject_data"][0]["file"] = orbpath

        # Security CHECKs
        statfile_path_gstyl = dir_dic["object"]["station"]["station_coordinates"]
        statfile_path_full = os.path.join(gins_path, statfile_path_gstyl[6:])

        ocloadfile_path_gstyl = dir_dic["object"]["station"]["ocean_tide_loading"]
        ocloadfile_path_full = os.path.join(gins_path, ocloadfile_path_gstyl[6:])

        gzgicmn.check_stat_in_statfile(stat, statfile_path_full)

        domes = gzgicmn.find_domes_in_statfile(stat, statfile_path_full)
        log.info("INFO : DOMES : %s", domes)
        gzgicmn.check_domes_in_oclo(domes[0], ocloadfile_path_full)

        # Case of kinematic process : need to change keys in user extension
        if "user_extension" in dir_dic:
            log.info("INFO : kinematic process with 'user_extension' fields")
            if "userext_gps__haute_freq" in dir_dic["user_extension"]:
                dir_dic["user_extension"]["userext_gps__haute_freq"] = stat
            if "userext_gps__haute_freq".upper() in dir_dic["user_extension"]:
                dir_dic["user_extension"]["userext_gps__haute_freq".upper()] = stat

            if "userext_addition" in dir_dic["user_extension"]:
                if "gps__haute_freq" in dir_dic["user_extension"]["userext_addition"]:
                    dir_dic["user_extension"]["userext_addition"][
                        "gps__haute_freq"
                    ] = stat
                if (
                    "gps__haute_freq".upper()
                    in dir_dic["user_extension"]["userext_addition"]
                ):
                    dir_dic["user_extension"]["userext_addition"][
                        "gps__haute_freq".upper()
                    ] = stat
        if "userext_addition" in dir_dic["user_extension"]:
            if np.any(
                [
                    "gps__haute_freq".upper() in e
                    for e in dir_dic["user_extension"]["userext_addition"]
                ]
            ):
                log.info("INFO : kinematic process with 'userext_addition' fields")
                index_gps_hf = [
                    "gps__haute_freq".upper() in e
                    for e in dir_dic["user_extension"]["userext_addition"]
                ].index(True)
                dir_dic["user_extension"]["userext_addition"][index_gps_hf] = (
                    "gps__haute_freq ".upper() + stat
                )

        # WRITING THE NEW DIRECTOR
        log.info("writing : %s", dir_out_path)
        with open(dir_out_path, "w+") as outfile:
            outfile.write(yaml.dump(dir_dic, default_flow_style=False))

        # correcting a bug : GINS YAML interpreter can't manage false or true ...
        utils.replace(dir_out_path, "false", "no")
        utils.replace(dir_out_path, "true", "yes")

        director_output_path_lis.append(dir_out_path)

    if multimode:
        if len(failed_rinex_date_lis) > 0:
            print("INFO : ", len(failed_rinex_date_lis), "/", N, " RINEXs failed")
            print("       ", [str(d) for d in failed_rinex_date_lis])
        return director_output_path_lis
    else:
        return dir_out_path


def _dir_regular_orbclk(dt_rinex_inp):

    if dt_rinex_inp <= dt.datetime.now() - dt.timedelta(days=14):
        prod_id = "G20"
    else:
        prod_id = "G20"

    ### ALL THE OTHER ORBITS/CLOCKS ARE NOW OBSOLETE (2025)
    # # GRG or GR2 for clock and orbits
    # if dt.datetime(1998, 1, 4) <= dt_rinex_inp <= dt.datetime(2013, 12, 29):
    #     prod_id = "GR2"
    # #    elif dt_rinex < dt.datetime(1998,1,4):
    # #        prod_id = 'GR1' # obsolete parce que GR1 intégralement inclu dans GR2
    # elif dt_rinex_inp < dt.datetime(1998, 1, 4):
    #     prod_id = "IG2"
    #     log.info("no GRG/GR2 orbit available for day %s", dt_rinex_inp)
    #     log.info("(day before 1998/1/4) => using IG2 orb. instead")
    # else:
    #     prod_id = "GRG"

    horpath = os.path.join("horloges", prod_id, "defaut")
    orbpath = os.path.join("orbites", prod_id, "defaut")

    return orbpath, horpath


def _dir_auto_staocl(stat_lower, dt_rinex, rinex_path, temp_data_folder):
    log.info("* Automatic stat and ocload files generation :")
    randid = "id" + str(np.random.randint(10000, 99999))
    stat_fname = utils.join_improved(
        "_",
        stat_lower,
        str(dt_rinex.year),
        str(conv.dt2doy(dt_rinex)),
        randid,
        ".stat",
    )

    stat_path_out = os.path.join(temp_data_folder, stat_fname)
    files_rw.write_station_file_gins_from_rinex(rinex_path, stat_path_out)
    oclo_fname = utils.join_improved(
        "_",
        stat_lower,
        str(dt_rinex.year),
        str(conv.dt2doy(dt_rinex)),
        randid,
        ".oclo",
    )

    oclo_path_out = os.path.join(temp_data_folder, oclo_fname)
    gzgicmn.write_oclo_file(stat_path_out, oclo_path_out)

    bool_statfil = os.path.isfile(stat_path_out)
    bool_oclofil = os.path.isfile(oclo_path_out)
    iwhile = 0
    while not (bool_statfil and bool_oclofil) or iwhile < 10:
        if iwhile > 6:
            log.warning("waiting for the ocean loading file for a while %i", iwhile)
        iwhile = iwhile + 1
        time.sleep(0.5)
        bool_statfil = os.path.isfile(stat_path_out)
        bool_oclofil = os.path.isfile(oclo_path_out)

    if not os.path.isfile(stat_path_out):
        log.error("no station file genrated !!!")
        stat_path_out = None

    if not os.path.isfile(oclo_path_out):
        log.error("no oceload file genrated !!!")
        oclo_path_out = None

    return stat_path_out, oclo_path_out
