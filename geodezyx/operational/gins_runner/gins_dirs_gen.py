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
import geodezyx.operational.gins_runner.gins_common  as gzgicmn
import geodezyx.operational.gins_runner.gins_prairie as gzgipra
import geodezyx.operational.gins_runner.gins_orbclk  as gzgiorb

#### Import the logger
import logging
log = logging.getLogger('geodezyx')

def gen_dirs_from_rnxs(
    rinex_paths_in,
    director_generik_path,
    director_name_prefix,
    out_director_folder="",
    temp_data_folder="",
    stations_file="",
    oceanload_file="",
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
    NEW FCT WHICH CAN MANAGE BOTH ONE RINEX OR A LIST OF RINEX

    produce a director from a rinex, and write it in the specific location

    Input:
        rinex_paths_in : can be a rinex path (a string) or a list of paths

    Return :
            the path of the produced director as a string

    If the rinex is not in a gins style folder, it will be copied in the
    ``temp_data_folder``

    If ``temp_data_folder`` is not specified ( == '') rinex will be copied in
    a ad hoc folder ``../gin/TEMP_DATA``

    If ``out_director_folder`` is not specified ( == ''), output director will be
    created in the ``../gin/data/directeur folder``

    auto_staocl : True or False
        create automatically a station file and a ocean loading file
        with the rinex header
        auto mode is prioritary upon the manu mode
        so, if a path for stat file or ocload file is specified and
        auto_staocl is on, it will be the automatic stat/ocload
        which will be used
    perso_orbclk :  True or False
        download and use specifics orbits
        according to ac & repro args
        (they are useless if perso_orbclk aren't activated)
    auto_interval : True or False
        find the interval in RINEX and apply it to the director
    out_coords : 'XYZ' or 'FLH'/'PLH'
        for geocentrical or geographical coords in output.
        any other string leaves the type of the generic director.
    prairie : True or False
        run prairie externally
        if True, prairie_kwargs control the arguments of the function
        prairie_manual (cf above)
    """

    # Multi or Single Mode ?
    if type(rinex_paths_in) is list:
        multimode = True
        rinex_path_lis = rinex_paths_in
        log.info("******** DIRECTORS GENRATION ********")
    elif type(rinex_paths_in) is str:
        multimode = False
        rinex_path_lis = [rinex_paths_in]
    else:
        log.error("check the rinex_paths_in !!!")
        return None

    N = len(rinex_path_lis)
    director_output_path_lis = []
    failed_rinex_date_lis = []

    for i, rinex_path in enumerate(rinex_path_lis):
        if multimode:
            log.info(" ======== ", i + 1, "/", N, " ======== ")
            log.info(" === ", os.path.basename(rinex_path))
        rinex_name = os.path.basename(rinex_path)
        dt_rinex = conv.rinexname2dt(rinex_name)
        gins_path = gzgicmn.get_gins_path()

        # be sure there is a TEMP DATA folder
        if temp_data_folder == "":
            temp_data_folder = os.path.join(gins_path, "gin", "TEMP_DATA")
        if not os.path.exists(temp_data_folder):
            os.makedirs(temp_data_folder)

        # be sure the RINEX is in gins folder ...
        bool_rinex_in_gin = gzgicmn.check_if_file_in_gin_folder(rinex_path)
        # ... and copy it otherwise
        if not bool_rinex_in_gin:
            log.info("INFO : will be copied in %s", temp_data_folder)
            rinex_path = gzgicmn.copy_file_in_gin_folder(rinex_path, temp_data_folder)
        # check if the RINEX is compressed ...
        bool_comp_rnx = operational.check_if_compressed_rinex(rinex_path)
        # ... if not crz2rnx !
        if bool_comp_rnx:
            crinex_path = rinex_path
            rinex_path = operational.crz2rnx(crinex_path, temp_data_folder)
            if not os.path.isfile(rinex_path):
                log.error("something went wrong during CRZ2RNX, skiping the file")
                log.error(crinex_path)

        # prairie extern
        if prairie:
            log.info("run prairie externally")
            pra_file_path = gzgipra.prairie_manual(
                rinex_path, temp_data_folder, **prairie_kwargs
            )

            if type(pra_file_path) is list and len(pra_file_path) == 0:
                log.error(
                    "something went wrong during prairie externally, skip the day"
                )
                log.error(rinex_path)
                failed_rinex_date_lis.append(dt_rinex)
                if multimode:
                    continue
                else:
                    return ""

        if out_director_folder == "":
            out_director_folder = os.path.join(
                gzgicmn.get_gins_path(), "gin", "data", "directeur"
            )

        stat_lower = rinex_name[0:4]
        stat = stat_lower.upper()

        if perso_orbclk:
            calccntr_suffix = "_" + ac + str(repro)
        else:
            calccntr_suffix = ""

        # date
        try:
            strt_epoch, end_epoch, freq_rnx_true = operational.rinex_start_end(
                rinex_path, True
            )
        except:
            log.error(
                "unable to get RINEX start/end, skiping  ... %s", rinex_path
            )

            failed_rinex_date_lis.append(dt_rinex)
            if multimode:
                continue
            else:
                return ""

        ses = "_" + operational.rinex_session_id(strt_epoch, end_epoch, full_mode=1)

        # Interval Initalisation
        # interval_line = utils.grep(rinex_path,'INTERVAL',True)
        #  float(interval_line.split()[0])
        if auto_interval:
            freq_rnx_str = "_" + str(int(freq_rnx_true)).zfill(2) + "s"
        else:
            freq_rnx_str = ""

        # coord type suffix
        if out_coords.upper() in ("FLH", "XYZ"):
            coord_prefix = "_" + out_coords.lower()
        else:
            coord_prefix = ""

        director_output_name = utils.join_improved("_",
                                                   director_name_prefix,
                                                   stat_lower,
                                                   str(conv.dt2jjul_cnes(dt_rinex)),
                                                   str(dt_rinex.year),
                                                   conv.dt2doy(dt_rinex),
                                                   freq_rnx_str,
                                                   coord_prefix,
                                                   ses,
                                                   calccntr_suffix,
                                                   ".yml")

        director_output_path = os.path.join(out_director_folder, director_output_name)
        director_generik_file = open(director_generik_path)
        dir_dic = yaml.load(director_generik_file,Loader=yaml.FullLoader)

        # ADDING RINEX NAME & DATES INTO THE DIR
        # name
        if not prairie:
            rnx_path_gstyl = gzgicmn.make_path_ginsstyle(rinex_path)
            dir_dic["observation"]["interobject_data"][1][
                "file"
            ] = rnx_path_gstyl
        else:
            pra_path_dir_gstyl = gzgicmn.make_path_ginsstyle(pra_file_path)
            dir_dic["observation"]["interobject_data"][1][
                "file"
            ] = pra_path_dir_gstyl

        # date
        strt_day, strt_sec = conv.dt2jjul_cnes(strt_epoch, False)
        end_day, end_sec = conv.dt2jjul_cnes(end_epoch, False)

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
            dir_dic["observation"]["removal"][
                "simulation_stepsize"
            ] = freq_rnx_true
            # check if its not a static dir
            if (
                dir_dic["parameter"]["adjustment_parameters"]["stations"][
                    "adjustment_frequency"
                ][-1]
                != 0
            ):
                dir_dic["parameter"]["adjustment_parameters"]["stations"][
                    "adjustment_frequency"
                ] = [0, 0, 0, freq_rnx_true]
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

        if auto_staocl == False and (oceanload_file == "" or stations_file == ""):
            log.warning(
                "Auto mode is off and no stat or oclo file specified !!!"
            )
        if auto_staocl:
            log.info("* stat and ocload files generation :")
            randid = "id" + str(np.random.randint(10000, 99999))
            sta_fname = utils.join_improved("_", stat_lower,
                                            str(dt_rinex.year),
                                            str(conv.dt2doy(dt_rinex)),
                                            randid, ".stat")

            stations_file = os.path.join(temp_data_folder, sta_fname)
            files_rw.write_station_file_gins_from_rinex(rinex_path, stations_file)
            oclo_fname =  utils.join_improved("_", stat_lower,
                                              str(dt_rinex.year),
                                              str(conv.dt2doy(dt_rinex)),
                                              randid, ".oclo")

            oceanload_file = os.path.join(temp_data_folder, oclo_fname)
            gzgicmn.write_oceanload_file(stations_file, oceanload_file)

            bool_statfil = os.path.isfile(stations_file)
            bool_oclofil = os.path.isfile(oceanload_file)
            iwhile = 0
            while not (bool_statfil and bool_oclofil) or iwhile < 6:
                if iwhile > 6:
                    log.warning("something is going wrong in the loop waiting for the ocean loading file ..., %i",
                        iwhile,
                    )
                iwhile = iwhile + 1
                time.sleep(0.5)
                bool_statfil = os.path.isfile(stations_file)
                bool_oclofil = os.path.isfile(oceanload_file)

            if not os.path.isfile(oceanload_file):
                log.error("no oceload file !!!")

        if oceanload_file != "":
            oceanload_file_ingin = gzgicmn.bring_to_gin_folder(oceanload_file, temp_data_folder)
            dir_dic["object"]["station"]["ocean_tide_loading"] = (
                gzgicmn.make_path_ginsstyle(oceanload_file_ingin)
            )
        if stations_file != "":
            stations_file_ingin = gzgicmn.bring_to_gin_folder(stations_file, temp_data_folder)
            dir_dic["object"]["station"]["station_coordinates"] = (
                gzgicmn.make_path_ginsstyle(stations_file_ingin)
            )

        # =========   ORBITS   =============
        if not perso_orbclk:
            # GRG or GR2 for clock and orbits
            if dt.datetime(1998, 1, 4) <= dt_rinex <= dt.datetime(2013, 12, 29):
                gr_gr2 = "GR2"
            #    elif dt_rinex < dt.datetime(1998,1,4):
            #        gr_gr2 = 'GR1' # obsolete parce que GR1 intégralement inclu dans GR2
            elif dt_rinex < dt.datetime(1998, 1, 4):
                gr_gr2 = "IG2"
                log.info("no GRG/GR2 orbit available for day %s", dt_rinex)
                log.info("(day before 1998/1/4) => using IG2 orb. instead")
            else:
                gr_gr2 = "GRG"
            horpath = os.path.join("horloges", gr_gr2, "defaut")
            orbpath = os.path.join("orbites", gr_gr2, "defaut")
        else:
            orbpath, horpath = gzgiorb.download_convert_2_gins_orb_clk(
                dt_rinex, temp_data_folder, ac=ac, repro=repro
            )

            # Exception case where no sp3/clk was found
            if orbpath is None or horpath is None:
                log.warning(
                    "something went wrong with orbits/clocks download/conversion %s",
                    dt_rinex,
                )
                log.warning("       skiping the day ...")
                failed_rinex_date_lis.append(dt_rinex)
                if multimode:
                    continue
                else:
                    return ""

            orbpath = gzgicmn.make_path_ginsstyle(orbpath)
            horpath = gzgicmn.make_path_ginsstyle(horpath)

        dir_dic["model"]["environment"]["gnss_clock"] = horpath
        dir_dic["observation"]["interobject_data"][0]["file"] = orbpath

        # correcting a pyyaml bug : version is interpreted as a int instead of a string
        #verslis = list(str(dir_dic["version"]))
        # verstr = verslis[0] + verslis[1] + "_" + "_".join(verslis[2:])
        # dir_dic["version"] = verstr
        # comment 202501

        # Security CHECKs
        statfile_path_gstyl = dir_dic["object"]["station"][
            "station_coordinates"
        ]
        statfile_path_full = os.path.join(gins_path, statfile_path_gstyl[6:])

        ocloadfile_path_gstyl = dir_dic["object"]["station"][
            "ocean_tide_loading"
        ]
        ocloadfile_path_full = os.path.join(gins_path, ocloadfile_path_gstyl[6:])

        gzgicmn.check_if_stat_in_stationfile(stat, statfile_path_full)

        domes = gzgicmn.find_DOMESstat_in_stationfile(stat, statfile_path_full)
        log.info("INFO : DOMES : %s", domes)
        gzgicmn.check_if_DOMES_in_oceanloadfile(domes[0], ocloadfile_path_full)

        # Case of kinematic process : need to change keys in user extension
        if "user_extension" in dir_dic:
            log.info("INFO : kinematic process with 'user_extension' fields")
            if "userext_gps__haute_freq" in dir_dic["user_extension"]:
                dir_dic["user_extension"]["userext_gps__haute_freq"] = stat
            if "userext_gps__haute_freq".upper() in dir_dic["user_extension"]:
                dir_dic["user_extension"]["userext_gps__haute_freq".upper()] = stat

            if "userext_addition" in dir_dic["user_extension"]:
                if (
                    "gps__haute_freq"
                    in dir_dic["user_extension"]["userext_addition"]
                ):
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
        log.info("writing : %s", director_output_path)
        with open(director_output_path, "w+") as outfile:
            outfile.write(yaml.dump(dir_dic, default_flow_style=False))

        # correcting a bug : GINS YAML interpreter can't manage false or true ...
        utils.replace(director_output_path, "false", "no")
        utils.replace(director_output_path, "true", "yes")

        director_output_path_lis.append(director_output_path)

    if multimode:
        if len(failed_rinex_date_lis) > 0:
            print("INFO : ", len(failed_rinex_date_lis), "/", N, " RINEXs failed")
            print("       ", [str(d) for d in failed_rinex_date_lis])
        return director_output_path_lis
    else:
        return director_output_path
