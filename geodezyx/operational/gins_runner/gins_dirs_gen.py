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
import pandas as pd
import yaml

#### geodeZYX modules
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import operational
from geodezyx import utils

#### Import geodezyx GINS submodules
import geodezyx.operational.gins_runner.gins_common as gynscmn
import geodezyx.operational.gins_runner.gins_prairie as gynspra
import geodezyx.operational.gins_runner.legacy.gins_orbclk_convert as gynsorbcvt

#### Import the logger
import logging

log = logging.getLogger("geodezyx")


def gen_dirs_rnxs(
        rnx_paths_inp,
        director_generik_path,
        director_name_prefix,
        out_director_folder=None,
        temp_data_folder=None,
        stations_file=None,
        oceanload_file=None,
        options_prairie_file=None,
        auto_stations_file=False,
        auto_oceanload=False,
        perso_orbclk=False,
        ac="igs",
        repro=2,
        auto_interval=True,
        out_coords="NULL",
        prairie=False,
        prairie_kwargs={"with_historik": 1, "with_wsb": 1},
        force=False,
        verbose=True,
        sites_id9=None,
        add_tropo_sol=True,
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
        an ad hoc folder ``../gin/TMP_GYNS``
    stations_file : str, optional
        Path to the stations file. Defaults to None.
    oceanload_file : str, optional
        Path to the ocean load file. Defaults to None.
    options_prairie_file : str, optional
        Path to the prairie options file. Defaults to None.
    auto_stations_file: bool, optional
    auto_oceanload: bool, optional
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
        Automatically find the interval in RINEX and apply it to the director.
         Defaults to False.
    out_coords : str, optional
        Output coordinates type ('XYZ' or 'FLH'/'PLH'). Defaults to "NULL".
    prairie : bool, optional
        Run prairie externally. Defaults to False.
        if True, prairie_kwargs control the arguments of the function
        prairie_manual (cf above)
    prairie_kwargs : dict, optional
        Arguments for the prairie_manual function.
        Defaults to {"with_historik": 1, "with_wsb": 1}.
    force : bool, optional
        Force the generation of the director even if a solution already exists.
        Defaults to False.
    verbose : bool
        verbose messages. Default is True.
    sites_id9 : iterable, optional
        A list, Pandas series... containing the site_id9 information.
        This is used to extract the site_id9 if the RINEX name is ambiguous.
    add_tropo_sol : bool, optional
        Add tropospheric solution keys to the director and
        then results in the listing and solution. Defaults to True.

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

    gin_path = gynscmn.get_gin_path()  # must remain not extended
    rnx_path_lis = list(sorted(rnx_path_lis))
    n_rnxs = len(rnx_path_lis)
    director_output_path_lis = []
    failed_rinex_date_lis = []
    temp_data_folder_lis = []

    dir_out_path = None
    tmp_fld_use = None

    def _fail_rnx(rnx_path_inp, rnx_dt_inp, message):
        log.error(message + ", skip %s", rnx_path_inp)
        failed_rinex_date_lis.append(rnx_dt_inp)
        bool_continue = False
        if multimode:
            bool_continue = True
        return bool_continue
        # else:
        #    raise Exception

    for i, rnx_path_ori in enumerate(rnx_path_lis):
        if multimode:
            log.info(f"======== {i + 1} / {n_rnxs} ======== ")
        log.info(" *** Generate directeur for: %s ", os.path.basename(rnx_path_ori))
        rnx_name = os.path.basename(rnx_path_ori)
        rnx_dt = conv.rinexname2dt(rnx_name)
        sites_id9_series = pd.Series(sites_id9)
        siteid9, siteid4_up, siteid4_lo = _dir_rnx_site_id(rnx_name, sites_id9_series)

        coord_prefix = (
            f"_{out_coords.lower()}" if out_coords.upper() in ("FLH", "XYZ") else ""
        )

        ac_suffix = f"_{ac}{repro}" if perso_orbclk else ""

        dir_name = utils.join_improved(
            "_",
            director_name_prefix,
            siteid9,
            str(rnx_dt.year),
            conv.dt2doy(rnx_dt),
            str(conv.dt2jjul_cnes(rnx_dt)),
            # freq_rnx_str,
            coord_prefix,
            # ses,
            ac_suffix,
        )

        sols_matching = gynscmn.check_solution(dir_name)

        if len(sols_matching) > 0 and not force:
            _fail_rnx(rnx_path_ori, rnx_dt, "solution already exists")
            continue

        # be sure there is a TEMP DATA folder
        if not temp_data_folder:
            temp_data_folder = gynscmn.get_temp_data_gins_path()
        if not os.path.exists(temp_data_folder):
            os.makedirs(temp_data_folder)

        # go for a sub TEMP folder
        tmp_fld_use = os.path.join(
            temp_data_folder, utils.get_timestamp() + "_" + dir_name
        )
        os.makedirs(tmp_fld_use)

        # be sure the RINEX is in gins folder ...
        rnx_path = gynscmn.bring_to_gin(rnx_path_ori, tmp_fld_use, verbose=verbose)

        # check if the RINEX is compressed ... if not crz2rnx !
        if operational.check_if_compressed_rinex(rnx_path):
            crinex_path = rnx_path
            rnx_path = operational.crz2rnx(crinex_path, tmp_fld_use, verbose=verbose)
            if not os.path.isfile(rnx_path):
                bool_cntu = _fail_rnx(rnx_path, rnx_dt, "CRZ2RNX failed")
                if bool_cntu:
                    continue

        # prairie extern
        pra_file_path = None
        if prairie:
            log.info("run prairie externally")
            pra_file_path = gynspra.prairie_manual(
                rnx_path, tmp_fld_use, **prairie_kwargs
            )

            if type(pra_file_path) is list and len(pra_file_path) == 0:
                bool_cntu = _fail_rnx(rnx_path, rnx_dt, "manual prairie failed")
                if bool_cntu:
                    continue

        if not out_director_folder:
            out_director_folder = os.path.join(gin_path, "gin", "data", "directeur")

        # date

        srt_epo, end_epo, freq_rnx = None, None, None

        try:
            srt_epo, end_epo, freq_rnx = operational.rinex_start_end(rnx_path, True, verbose=verbose)
        except Exception:
            log.warning("Unable to read RINEX, fallback to filename to get its start/end")
            srt_epo = rnx_dt
            end_epo = srt_epo + dt.timedelta(seconds=86360)  # 1 day
            freq_rnx = 30

        if not srt_epo and not end_epo and not freq_rnx:
            srt_epo, end_epo, freq_rnx = None, None, None
            bool_cntu = _fail_rnx(rnx_path, rnx_dt, "get RINEX start/end failed")
            if bool_cntu:
                continue

        dir_out_fname = dir_name + ".yml"
        dir_out_path = os.path.join(out_director_folder, dir_out_fname)
        director_generik_file = open(director_generik_path)
        dir_dic = yaml.load(director_generik_file, Loader=yaml.FullLoader)

        # ADDING RINEX NAME & DATES INTO THE DIR
        # name
        if not prairie:
            rnx_path_gstyl = gynscmn.make_path_ginsstyle(rnx_path)
            dir_dic["observation"]["interobject_data"][1]["file"] = rnx_path_gstyl
        else:
            pra_path_dir_gstyl = gynscmn.make_path_ginsstyle(pra_file_path)
            dir_dic["observation"]["interobject_data"][1]["file"] = pra_path_dir_gstyl

        # date
        strt_day, strt_sec = conv.dt2jjul_cnes(srt_epo, False)
        end_day, end_sec = conv.dt2jjul_cnes(end_epo, False)

        # SPOTGINS compatibility for the full day
        # one more sec will be +/- below
        strt_sec = 19
        end_sec =  86389.0

        dir_dic["date"]["arc_start"][0] = strt_day
        dir_dic["date"]["arc_stop"][0] = end_day
        if "initial_state_vector_date" in list(dir_dic["date"].keys()):
            dir_dic["date"]["initial_state_vector_date"][0] = strt_day

        dir_dic["date"]["arc_start"][1] = (
                strt_sec - 1
        )  # -1 to make it SPOTGINS compatible
        dir_dic["date"]["arc_stop"][1] = (
                end_sec + 1
        )  # +1 to make it SPOTGINS compatible
        if "initial_state_vector_date" in list(dir_dic["date"].keys()):
            dir_dic["date"]["initial_state_vector_date"][1] = strt_sec

        # interval
        if auto_interval:  # dir is modified only if freq_rnx >= 0
            dir_dic = _dir_auto_intrvl(dir_dic, freq_rnx)

        # output coords
        if out_coords.upper() in ("FLH", "PLH"):
            dir_dic["parameter"]["adjustment_parameters"]["stations"][
                "coordinates"
            ] = "latitude_longitude_height"
        elif out_coords.upper() == "XYZ":
            dir_dic["parameter"]["adjustment_parameters"]["stations"][
                "coordinates"
            ] = "cartesian_xyz"

        # ============== STATION FILE & OCEANLOAD FILE ==============
        # station file and oceanload file , auto or manu mode
        # auto mode is prioritary upon the manu mode

        if not auto_stations_file and not stations_file:
            log.warning("auto_stations_file is off and no stations_file specified !!!")

        if not auto_oceanload and not oceanload_file:
            log.warning("auto_oceanload is off and no oceanload_file specified !!!")

        if auto_stations_file:
            stations_file = _dir_auto_stfi(siteid4_lo, rnx_dt, rnx_path, tmp_fld_use)

        if auto_oceanload:
            oceanload_file = _dir_auto_oclo(
                siteid4_lo, rnx_dt, tmp_fld_use, stations_file
            )

        if stations_file:
            stfi_ingin = gynscmn.bring_to_gin(str(stations_file), tmp_fld_use, verbose=verbose)
            dir_dic["object"]["station"]["station_coordinates"] = (
                gynscmn.make_path_ginsstyle(stfi_ingin)
            )
        if oceanload_file:
            oclo_ingin = gynscmn.bring_to_gin(str(oceanload_file), tmp_fld_use, verbose=verbose)
            dir_dic["object"]["station"]["ocean_tide_loading"] = (
                gynscmn.make_path_ginsstyle(oclo_ingin)
            )

        # ============== PRAIRIE OPTIONS ==============
        # NB : prairie options path does not have to be in gin style
        # but in its absolute path
        if options_prairie_file:
            dir_dic["model"]["environment"]["gnss_preprocessing_options"] = (
                os.path.realpath(options_prairie_file)
            )

        # ============== ORBITS/CLOCKS ================
        if perso_orbclk:
            orbpath, horpath = gynsorbcvt.download_convert_2_gins_orb_clk(
                rnx_dt, tmp_fld_use, ac=ac, repro=repro
            )
        else:
            orbpath, horpath = _dir_regular_orbclk(rnx_dt)

        # Exception case where no sp3/clk was found
        if orbpath is None or horpath is None:
            bool_cntu = _fail_rnx(rnx_path, rnx_dt, "orbits/clocks download/conversion")
            if bool_cntu:
                continue

        if perso_orbclk:
            orbpath = gynscmn.make_path_ginsstyle(orbpath)
            horpath = gynscmn.make_path_ginsstyle(horpath)

        dir_dic["model"]["environment"]["gnss_clock"] = horpath
        dir_dic["observation"]["interobject_data"][0]["file"] = orbpath

        # Security CHECKs
        stfi_path_gstyl = dir_dic["object"]["station"]["station_coordinates"]
        stfi_path_full = str(os.path.join(gin_path, stfi_path_gstyl[6:]))

        oclo_path_gstyl = dir_dic["object"]["station"]["ocean_tide_loading"]
        oclo_path_full = str(os.path.join(gin_path, oclo_path_gstyl[6:]))

        gynscmn.check_site_stfl(siteid4_up, stfi_path_full)

        domes = gynscmn.find_domes_in_stfl(siteid4_up, stfi_path_full)
        if verbose:
            log.info("DOMES: %s", domes)
        gynscmn.chek_domes_oclo(domes[0], oclo_path_full)

        # ========= CUSTOM USER EXTENSION (tropo, high freq...) =============
        dir_dic = _dir_userext(dir_dic, siteid9, add_tropo_sol)

        # WRITING THE NEW DIRECTOR
        if verbose:
            log.info("writing: %s", dir_out_path)
        with open(dir_out_path, "w+") as outfile:
            outfile.write(yaml.dump(dir_dic, default_flow_style=False))

        # correcting a bug : GINS YAML interpreter can't manage false or true ...
        utils.replace(dir_out_path, "false", "no")
        utils.replace(dir_out_path, "true", "yes")

        director_output_path_lis.append(dir_out_path)
        temp_data_folder_lis.append(tmp_fld_use)

    if multimode:
        if len(failed_rinex_date_lis) > 0:
            print("INFO : ", len(failed_rinex_date_lis), "/", n_rnxs, " RINEXs failed")
            print("       ", [str(d) for d in failed_rinex_date_lis])
        return director_output_path_lis, temp_data_folder_lis
    else:
        return dir_out_path, tmp_fld_use


def _dir_regular_orbclk(dt_rinex_inp):
    """
    Get the regular orbits and clocks path for a given date
    """

    if dt_rinex_inp <= dt.datetime.now() - dt.timedelta(days=14):
        prod_id = "G20"
    else:
        prod_id = "G20"

    ### ALL THE OTHER ORBITS/CLOCKS ARE NOW OBSOLETE (2025)

    horpath = os.path.join("horloges", prod_id, "defaut")
    orbpath = os.path.join("orbites", prod_id, "defaut")

    return orbpath, horpath


def _dir_auto_intrvl(dir_dic, freq_rnx):
    """
    Automatically set the interval in the director
    """
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
        "auto_interval, simulation_stepsize : %s",
        dir_dic["observation"]["removal"]["simulation_stepsize"],
    )

    return dir_dic


def _dir_auto_stfi(stat_lower, dt_rinex, rinex_path, temp_data_folder):
    """
    Generate an ad hoc station file from a RINEX file
    """
    log.info("* Automatic station file generation")
    randid = "id" + str(np.random.randint(10000, 99999))
    stat_fname = utils.join_improved(
        "_",
        stat_lower,
        str(dt_rinex.year),
        str(conv.dt2doy(dt_rinex)),
        randid,
        ".stat",
    )

    stfi_path_out = os.path.join(temp_data_folder, stat_fname)
    files_rw.write_station_file_gins_from_rinex(rinex_path, stfi_path_out)
    bool_statfil = os.path.isfile(stfi_path_out)
    if not os.path.isfile(stfi_path_out):
        log.error("no station file genrated !!!")
        stfi_path_out = None

    return stfi_path_out


def _dir_auto_oclo(stat_lower, dt_rinex, temp_data_folder, stat_path_inp):
    """
    Generate an ad hoc ocean loading file from a station file
    """
    log.info("* Automatic stat and ocean loading file generation")
    randid = "id" + str(np.random.randint(10000, 99999))

    oclo_fname = utils.join_improved(
        "_",
        stat_lower,
        str(dt_rinex.year),
        str(conv.dt2doy(dt_rinex)),
        randid,
        ".oclo",
    )

    oclo_path_out = os.path.join(temp_data_folder, oclo_fname)
    gynscmn.write_oclo_file(stat_path_inp, oclo_path_out)

    bool_oclofil = os.path.isfile(oclo_path_out)
    iwhile = 0
    while not bool_oclofil or iwhile < 10:
        if iwhile > 6:
            log.warning("waiting for the ocean loading file for a while %i", iwhile)
        iwhile = iwhile + 1
        time.sleep(0.5)
        bool_oclofil = os.path.isfile(oclo_path_out)

    if not os.path.isfile(oclo_path_out):
        log.error("no oceload file genrated !!!")
        oclo_path_out = None

    return oclo_path_out


def _dir_rnx_site_id(rnx_name, sites_id9_series):
    """
    Extract the site_id9 from a Series (from e.g. a SPOTGINS master file)
    based on the RINEX name
    """
    if conv.rinex_regex_search_tester(
            rnx_name, short_name=False, long_name=True
    ):  ### RINEX3
        site_id9 = rnx_name[0:9].upper()
    else:  ### RINEX2
        site_id4 = rnx_name[0:4].upper()
        site_id9 = site_id4 + "00XXX"
        if type(sites_id9_series) is pd.Series:
            ser_bool = sites_id9_series.str[:4].str.match(site_id4)
            if ser_bool.any():
                site_id9 = sites_id9_series.loc[ser_bool].values[0]
            if ser_bool.sum() > 1:
                log.warning("more than one site_id9 found for %s: %s", site_id4,
                            sites_id9_series.loc[ser_bool].values)
                log.warning("taking the first one : %s", site_id9)

    site_id4_upper = site_id9[0:4]
    site_id4_lower = site_id4_upper.lower()

    return site_id9, site_id4_upper, site_id4_lower


def _dir_userext(dir_dic, site_id9, add_tropo_sol=True):
    """
    Customize the user_extension fields with the site_id4_upper
    """
    uext = "user_extension"
    uadd = "userext_addition"

    if uext in dir_dic:
        if add_tropo_sol:
            dir_dic[uext][uadd].append("ADD_TROPO_IN_SOLUTION")
            dir_dic[uext][uadd].append("GPS__AFF_SINEX")

        # GPS__HAUTE_FREQ STAT is not a key/value couple, but is considered a single string
        # we must search this string and concatenate it with the right site code
        if uadd in dir_dic[uext]:
            is_gpshf = ["GPS__HAUTE_FREQ" in e for e in dir_dic[uext][uadd]]
            if np.any(is_gpshf):
                idx_gpshf = is_gpshf.index(True)
                dir_dic[uext][uadd][idx_gpshf] = "GPS__HAUTE_FREQ " + site_id9.upper()

    return dir_dic


def _dir_auto_staocl(stat_lower, dt_rinex, rinex_path, temp_data_folder):
    ##### DEPRECATED, NOW SEPARATED !!!
    log.info("* Automatic stat and ocload files generation")
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
    gynscmn.write_oclo_file(stat_path_out, oclo_path_out)

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
