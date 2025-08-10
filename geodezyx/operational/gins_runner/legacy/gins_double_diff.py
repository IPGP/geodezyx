#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26/02/2025 11:03:05

@author: psakic
"""

import os
import numpy as np
from geodezyx import conv, utils
import pandas as pd
import yaml
import copy
import shutil
import datetime as dt

#### Import the logger
import logging
log = logging.getLogger('geodezyx')


#################### DOUBLE DIFF ####################


def gen_dirs_from_double_diff(
    dd_files_paths_in,
    director_generik_path,
    director_name_prefix,
    out_director_folder="",
    temp_data_folder="",
    stations_file="",
    oceanload_file="",
    auto_staocl=False,
    perso_orbclk=False,
    gins_style_orb="GRG",
    ac="igs",
    repro=2,
    auto_interval=True,
    out_coords="NULL",
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
    if type(dd_files_paths_in) is list:
        multimode = True
        dd_files_paths_in = dd_files_paths_in
        print("******** DIRECTORS GENRATION ********")
    elif type(dd_files_paths_in) is str:
        multimode = False
        dd_files_paths_in = [dd_files_paths_in]
    else:
        print("ERR : gen_dirs_from_double_diff : check the dd_files_paths_in !!!")
        return None

    N = len(dd_files_paths_in)
    director_output_path_lis = []
    failed_dd_date_lis = []

    for i, dd_path in enumerate(dd_files_paths_in):
        if multimode:
            print(" ======== ", i + 1, "/", N, " ======== ")
            print(" === ", os.path.basename(dd_path))
        # names, dates etc ...
        Table_dd = pd.read_table(dd_path, sep=r"\s+", header=-1)

        dd_name = os.path.basename(dd_path)

        dd_start_jjul = np.min(Table_dd[12])
        dd_start_sec = np.min(Table_dd[13])
        dd_end_jjul = np.max(Table_dd[12])
        dd_end_sec = np.max(Table_dd[13])

        dd_start_epoch = conv.jjul_cnes2dt(dd_start_jjul) + dt.timedelta(
            seconds=dd_start_sec
        )
        dd_end_epoch = conv.jjul_cnes2dt(dd_end_jjul) + dt.timedelta(seconds=dd_end_sec)

        freq_dd_true = np.min(np.abs(np.diff(Table_dd[13])))

        # stat_A_lower = np.unique(Table_dd[10])[0]
        # stat_A = stat_A_lower.upper()
        # stat_B_lower = np.unique(Table_dd[11])[0]
        # stat_B = stat_B_lower.upper()

        stats_lower_list = list(set(list(Table_dd[10]) + list(Table_dd[11])))
        stats_list = [st.upper() for st in stats_lower_list]

        gins_path = get_gins_path()

        # be sure there is a TEMP DATA folder
        if temp_data_folder == "":
            temp_data_folder = os.path.join(gins_path, "gin", "TEMP_DATA")
        if not os.path.exists(temp_data_folder):
            os.makedirs(temp_data_folder)

        # be sure the RINEX is in gins folder ...
        bool_rinex_in_gin = check_if_file_in_gin_folder(dd_path)
        # ... and copy it otherwise
        if not bool_rinex_in_gin:
            print("INFO : will be copied in ", temp_data_folder)
            dd_path = copy_file_in_gin_folder(dd_path, temp_data_folder)

        if out_director_folder == "":
            out_director_folder = os.path.join(
                get_gins_path(), "gin", "data", "directeur"
            )

        if perso_orbclk:
            calccntr_suffix = "_" + ac + str(repro)
        else:
            calccntr_suffix = ""

        # coord type suffix
        if out_coords.upper() in ("FLH", "XYZ"):
            coord_prefix = "_" + out_coords.lower()
        else:
            coord_prefix = ""

        stats_str = "_".join(stats_list)
        director_output_name = (
            director_name_prefix
            + "_"
            + stats_str
            + "_"
            + str(dd_start_jjul)
            + "_"
            + str(dd_start_epoch.year)
            + "_"
            + conv.dt2doy(dd_start_epoch)
            + "_"
            + str(freq_dd_true)
            + coord_prefix
            + calccntr_suffix
            + ".yml"
        )

        director_output_path = os.path.join(out_director_folder, director_output_name)
        director_generik_file = open(director_generik_path)
        director_dic = yaml.load(director_generik_file)
        director_dic_proto = copy.deepcopy(director_dic)

        # ADDING RINEX NAME & DATES INTO THE DIR
        # name
        dd_path_dir_compatible = make_path_ginsstyle(dd_path)
        director_dic["observation"]["interobject_data"][1][
            "file"
        ] = dd_path_dir_compatible

        # date
        strt_day, strt_sec = conv.dt2jjul_cnes(dd_start_epoch, False)
        end_day, end_sec = conv.dt2jjul_cnes(dd_end_epoch, False)

        director_dic["date"]["arc_start"][0] = strt_day
        director_dic["date"]["arc_stop"][0] = end_day
        if "initial_state_vector_date" in list(director_dic["date"].keys()):
            director_dic["date"]["initial_state_vector_date"][0] = strt_day

        director_dic["date"]["arc_start"][1] = float(strt_sec)
        director_dic["date"]["arc_stop"][1] = float(end_sec)
        if "initial_state_vector_date" in list(director_dic["date"].keys()):
            director_dic["date"]["initial_state_vector_date"][1] = float(strt_sec)

        # interval
        if auto_interval:  # dir is modified only if freq_rnx >= 0
            director_dic["observation"]["removal"]["simulation_stepsize"] = int(
                freq_dd_true
            )
            # check if its not a static dir
            if (
                director_dic["parameter"]["adjustment_parameters"]["stations"][
                    "adjustment_frequency"
                ][-1]
                != 0
            ):
                director_dic["parameter"]["adjustment_parameters"]["stations"][
                    "adjustment_frequency"
                ] = [0, 0, 0, freq_dd_true]
        print(
            "INFO : auto_interval",
            auto_interval,
            ", simulation_stepsize : ",
            director_dic["observation"]["removal"]["simulation_stepsize"],
        )
        # output coords
        if out_coords.upper() in ("FLH", "PLH"):
            director_dic["parameter"]["adjustment_parameters"]["stations"][
                "coordinates"
            ] = "latitude_longitude_height"
        elif out_coords.upper() == "XYZ":
            director_dic["parameter"]["adjustment_parameters"]["stations"][
                "coordinates"
            ] = "cartesian_xyz"

        # =========   ORBITS   =============
        if not perso_orbclk:
            # GRG or GR2 for clock and orbits
            if dt.datetime(1998, 1, 4) <= dd_start_epoch <= dt.datetime(2013, 12, 29):
                gr_gr2 = "GR2"
            #    elif dt_rinex < dt.datetime(1998,1,4):
            #        gr_gr2 = 'GR1' # obsolete parce que GR1 intégralement inclu dans GR2
            else:
                gr_gr2 = "GRG"
            horpath = os.path.join("horloges", gr_gr2, "defaut")
            orbpath = os.path.join("orbites", gr_gr2, "defaut")
        else:
            orbpath, horpath = download_convert_2_gins_orb_clk(
                dd_start_epoch, temp_data_folder, ac=ac, repro=repro
            )

            # Exception case where no sp3/clk was found
            if orbpath == None or horpath == None:
                print(
                    "WARN : something went wrong with orbits/clocks download/conversion "  )
                print("       skiping the day ...")
                failed_dd_date_lis.append(dd_start_epoch)
                if multimode:
                    continue
                else:
                    return ""

            orbpath = make_path_ginsstyle(orbpath)
            horpath = make_path_ginsstyle(horpath)

        director_dic["model"]["environment"]["gnss_clock"] = horpath
        director_dic["observation"]["interobject_data"][0]["file"] = orbpath

        # correcting a pyyaml bug : version is interpreted as a int instead of a string
        verslis = list(str(director_dic["version"]))

        verstr = verslis[0] + verslis[1] + "_" + "_".join(verslis[2:])
        director_dic["version"] = verstr

        # Security CHECKs
        statfile_path_ginsstyle = director_dic["object"]["station"][
            "station_coordinates"
        ]
        statfile_path_full = os.path.join(gins_path, statfile_path_ginsstyle[6:])

        ocloadfile_path_ginsstyle = director_dic["object"]["station"][
            "ocean_tide_loading"
        ]
        ocloadfile_path_full = os.path.join(gins_path, ocloadfile_path_ginsstyle[6:])

        # check_stat_in_statfile(stat,statfile_path_full)

        # domes = find_domes_in_statfile(stat,statfile_path_full)
        # print 'INFO : DOMES : ',domes
        # chek_domes_oclo(domes[0],ocloadfile_path_full)

        # Case of kinematic process : need to change keys in user extension
        stat_B = None
        if "user_extension" in director_dic:
            print("INFO : kinematic process with 'user_extension' fields")
            if "userext_gps__haute_freq" in director_dic["user_extension"]:
                director_dic["user_extension"]["userext_gps__haute_freq"] = stat_B
            if "userext_gps__haute_freq".upper() in director_dic["user_extension"]:
                director_dic["user_extension"][
                    "userext_gps__haute_freq".upper()
                ] = stat_B

            if "userext_addition" in director_dic["user_extension"]:
                if (
                    "gps__haute_freq"
                    in director_dic["user_extension"]["userext_addition"]
                ):
                    director_dic["user_extension"]["userext_addition"][
                        "gps__haute_freq"
                    ] = stat_B
                if (
                    "gps__haute_freq".upper()
                    in director_dic["user_extension"]["userext_addition"]
                ):
                    director_dic["user_extension"]["userext_addition"][
                        "gps__haute_freq".upper()
                    ] = stat_B

        if "userext_addition" in director_dic["user_extension"]:
            if np.any(
                [
                    "gps__haute_freq".upper() in e
                    for e in director_dic["user_extension"]["userext_addition"]
                ]
            ):
                print("INFO : kinematic process with 'userext_addition' fields")
                index_gps_hf = [
                    "gps__haute_freq".upper() in e
                    for e in director_dic["user_extension"]["userext_addition"]
                ].index(True)
                director_dic["user_extension"]["userext_addition"][index_gps_hf] = (
                    "gps__haute_freq ".upper() + stat_B
                )

        # return director_dic_proto , director_dic

        # WRITING THE NEW DIRECTOR
        print("INFO : writing : ", director_output_path)
        with open(director_output_path, "w+") as outfile:
            outfile.write(yaml.dump(director_dic, default_flow_style=False))

        # correcting a bug : GINS YAML interpreter can't manage false or true ...
        utils.replace(director_output_path, "false", "no")
        utils.replace(director_output_path, "true", "yes")

        director_output_path_lis.append(director_output_path)

    if multimode:
        if len(failed_dd_date_lis) > 0:
            print("INFO : ", len(failed_dd_date_lis), "/", N, " FILEs failed")
            print("       ", [str(d) for d in failed_dd_date_lis])
        return director_output_path_lis
    else:
        return director_output_path


def double90_runner(inp_prarie_file, temp_data_folder=""):

    # be sure there is a TEMP DATA folder
    if temp_data_folder == "":
        temp_data_folder = os.path.join(get_gins_path(1), "TEMP_DATA")
    if not os.path.exists(temp_data_folder):
        os.makedirs(temp_data_folder)

    curdir = os.getcwd()
    os.chdir(temp_data_folder)

    command = "double90 " + inp_prarie_file
    print("INFO : double90 command lanched :")
    print(command)
    print("")

    p = subprocess.Popen(
        [command], shell=True
    )  # , stdout=subprocess.PIPE , stderr=subprocess.PIPE)
    p.wait()
    outname = os.path.basename(inp_prarie_file).split(".")[0] + ".dd"
    os.rename("./fort.20", outname)
    os.chdir(curdir)

    outpath = os.path.join(temp_data_folder, outname)

    return outpath


def double_diff_binom(
    rinex_path_A, rinex_path_B, temp_data_folder="", final_data_folder=""
):
    """
    DISCONTINUED
    """

    # be sure there is a TEMP DATA folder
    if temp_data_folder == "":
        temp_data_folder = os.path.join(get_gins_path(1), "TEMP_DATA")
    if not os.path.exists(temp_data_folder):
        os.makedirs(temp_data_folder)

    statA = os.path.basename(rinex_path_A)[0:4]
    statB = os.path.basename(rinex_path_B)[0:4]

    strt_epochA, end_epochA, freq_rnx_trueA = operational.rinex_start_end(
        rinex_path_A, True
    )
    strt_epochB, end_epochB, freq_rnx_trueB = operational.rinex_start_end(
        rinex_path_B, True
    )

    year, doy = conv.dt2doy_year(strt_epochA)

    argdict = dict()
    # argdict['-ignore']  = 'EG'

    praA_path = prairie_manual(
        rinex_path_A, temp_data_folder=temp_data_folder, argsdict=argdict
    )
    praB_path = prairie_manual(
        rinex_path_B, temp_data_folder=temp_data_folder, argsdict=argdict
    )

    praABcat_name = utils.join_improved("_", statA, statB, year, doy) + ".pra.cat"
    praABcat_path = os.path.join(temp_data_folder, praABcat_name)

    utils.cat(praABcat_path, praA_path, praB_path)

    dd_path = double90_runner(praABcat_path, temp_data_folder=temp_data_folder)

    if final_data_folder != "":
        shutil.move(dd_path, final_data_folder)
        dd_path = os.path.join(final_data_folder, os.path.basename(dd_path))

    return dd_path


def double_diff_multi(
    rinex_path_list,
    temp_data_folder="",
    final_data_folder="",
    remove=True,
    ignore_glo_gal=False,
):
    """
    rinex_path_list :
        List of RINEXs paths

    temp_data_folder & final_data_folder :
        if empty string (''), files are saved in pygins_runner "TEMP_DATA" folder
        and in the specified folder else

    """

    # be sure there is a TEMP DATA folder
    if temp_data_folder == "":
        temp_data_folder = os.path.join(get_gins_path(1), "TEMP_DATA")
    if not os.path.exists(temp_data_folder):
        os.makedirs(temp_data_folder)

    prarie_path_list = []
    stats_list = []

    for rinex_path in rinex_path_list:

        print("****** conversion to GINS (PRARIE) format of :")
        print(rinex_path)

        bool_comp_rnx = operational.check_if_compressed_rinex(rinex_path)
        # ... if not crz2rnx !
        if bool_comp_rnx:
            crinex_path = rinex_path
            rinex_path = operational.crz2rnx(crinex_path, temp_data_folder)
            if not os.path.isfile(rinex_path):
                print("ERR : something went wrong during CRZ2RNX, skiping the file")
                print(crinex_path)
                continue

        statA = os.path.basename(rinex_path)[0:4]

        strt_epochA, end_epochA, freq_rnx_trueA = operational.rinex_start_end(
            rinex_path, True
        )
        doy, year = conv.dt2doy_year(strt_epochA)

        argdict = dict()

        if ignore_glo_gal:
            argdict["-ignore"] = "EG"

        praA_path = prairie_manual(
            rinex_path, temp_data_folder=temp_data_folder, argsdict=argdict
        )

        prarie_path_list.append(praA_path)
        stats_list.append(statA)

        if bool_comp_rnx:
            os.remove(rinex_path)

    pracat_name = (
        utils.join_improved("_", year, doy, *list(set(stats_list))) + ".pra.cat"
    )
    pracat_path = os.path.join(temp_data_folder, pracat_name)

    print("****** Concatenation of the PRAIRIE files in :")
    print(pracat_path)

    utils.cat(pracat_path, *prarie_path_list)

    print("****** generation of the Double Diff. file :")

    dd_path = double90_runner(pracat_path, temp_data_folder=temp_data_folder)

    print("double diff. file generated in :")
    print(dd_path)

    if final_data_folder != "":
        shutil.move(dd_path, final_data_folder)
        dd_path = os.path.join(final_data_folder, os.path.basename(dd_path))

    if remove:
        os.remove(pracat_path)
        [os.remove(e) for e in prarie_path_list]

    return dd_path


# ==========================================================================
# ==========================================================================


# def check_if_director_path_is_ok(director_path,gins_path=get_gin_path()):
# USELESS FOR THE MOMENT
#    # 2 choices : the director_path is absolute or relative
#    director_name = os.path.basename(director_path)
#    nominal_director_path =  os.path.join(gins_path,'gin/data/directeur',director_name)
#
#    if os.path.isabs(adir):
#        out = True
#
#        if not filecmp.cmp(nominal_director_path , director_path):
#            print 'WARN : nominal_director_path != director_path !!!'
#            print director_path
#            print nominal_director_path
#            out = False
#    else:
#        if os.path.isfile(nominal_director_path):
#            out = True
#        else:
#            print 'WARN : director_path dont exists !!!'
#            out = False
#    return out
