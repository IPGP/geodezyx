# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.operational contains functions to run the 
GNSS processing software GINS. 

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
import collections
import copy
import datetime as dt
import glob
#### Import the logger
import logging
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import yaml
from bs4 import UnicodeDammit  # finding the right encoding for debug

#### geodeZYX modules
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import operational
from geodezyx import time_series
from geodezyx import utils

log = logging.getLogger(__name__)

##########  END IMPORT  ##########


def get_gins_path(extended=False):
    try:
        gs_user = os.environ["GS_USER"]
    except:
        print("ERR : env. var. $GS_USER dont exists !!!")
        return None
        # gs_user =  os.environ['HOME']
    if not extended:
        return gs_user
    else:
        return os.path.join(os.environ["GS_USER"], "gin")


def check_if_stat_in_stationfile(stat, stationfile):
    boolout = utils.check_regex(stationfile, stat)
    if not boolout:
        print("WARN :", stat, "not in", stationfile)
        print("       check your RINEX header and particularly")
        print("       the MARKER NAME field (station 4-char. code)")
    return boolout


def check_if_DOMES_in_oceanloadfile(domes, oclofile):
    boolout = utils.check_regex(oclofile, "^   " + str(domes))
    if not boolout:
        print("WARN :", domes, "not in", oclofile)
    return boolout


def find_DOMESstat_in_stationfile(stat, stationfile):
    fil = open(stationfile)
    for l in fil:
        if stat in l:
            f = l.split()
            return f[0], f[1]
    return 00000, "M000"


def check_if_file_in_gin_folder(a_file_path, gins_path=""):
    if gins_path == "":
        gins_path = os.path.join(get_gins_path(), "gin")

    real_path_file = os.path.realpath(a_file_path)
    real_path_gins = os.path.realpath(gins_path)

    if real_path_gins in real_path_file:
        boolout = True
    else:
        boolout = False

    #    a_file_path_splited = a_file_path.split('/')
    #    try:
    #        i_gin = a_file_path_splited.index('gin')
    #    except:
    #        i_gin = 9999
    #    potential_gin_path_splited = ['/'] + a_file_path_splited[:i_gin+1]
    #    potential_gin_path = os.path.join(*potential_gin_path_splited)
    #
    #    boolout = os.path.samefile(gins_path,potential_gin_path)

    if not boolout:
        print("WARN : ", a_file_path, " not in ", gins_path)
    #        print 'true files:', real_path_file , '&' , real_path_gins
    #        print '(potential gins path : ',potential_gin_path,')'
    return boolout


def check_good_exec_of_GINS(streamin, director_name):
    if "Exécution terminée du fichier" in streamin.read():
        print("INFO : happy end for " + director_name + " :)")
        return True
    else:
        print("WARN : bad end for " + director_name + " :(")
        return False


def make_path_ginsstyle(pathin):
    """input path must be an absolute path with /gin/ inside
    output will be .temp.gin/<rest of the path>"""
    # VERY WEAK SOLUTION ....

    if pathin.startswith(".temp.gin"):
        print("INFO : ", pathin, "already a director compatible path")
        return pathin

    if not "/gin/" in pathin:
        print("ERR : not /gin/ in", pathin)
        return None

    rnx_path_lis = pathin.split("/")
    i_gin = rnx_path_lis.index("gin")
    rnx_path_lis2 = rnx_path_lis[i_gin:]
    rnx_path_lis2[0] = ".temp.gin"
    pathout = os.path.join(*rnx_path_lis2)
    return pathout


def copy_file_in_gin_folder(
    a_file_path, repository_folder="/home/psakicki/gin/TP/GWADA/RINEX"
):
    if repository_folder == "":
        repository_folder = os.path.join(get_gins_path(), "TEMP_DATA")
    if not os.path.exists(repository_folder):
        os.makedirs(repository_folder)

    shutil.copy2(a_file_path, repository_folder)

    out = os.path.join(repository_folder, os.path.basename(a_file_path))

    return out


def write_oceanload_file(station_file, oceanload_out_file, fes_yyyy=2004):
    temp_cmd_fil = os.path.join(os.path.dirname(oceanload_out_file), "oclo.cmd.tmp")
    temp_cmd_filobj = open(temp_cmd_fil, "w")

    temp_cmd_filobj.write(station_file + "\n")
    temp_cmd_filobj.write(oceanload_out_file + "\n")
    if fes_yyyy == 2012:
        exe_loadoce_cmd = "exe_loadoce_fes2012"
    elif fes_yyyy == 2004:
        exe_loadoce_cmd = "exe_loadoce"

    p = subprocess.Popen(
        exe_loadoce_cmd + " < " + temp_cmd_fil,
        shell=True,
        executable="/bin/bash",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    #    os.remove(temp_cmd_fil)
    temp_cmd_filobj.close()

    if True:
        return oceanload_out_file
    else:
        return ""


# def generate_director_from_rinex(rinex_path,director_generik_path,
#                                 director_name_prefix,out_director_folder='',
#                                 temp_data_folder='',stations_file='',
#                                 oceanload_file='',auto_staocl=False,
#                                 perso_orbclk=False,ac='igs',repro=2,
#                                 auto_interval=True):
#
#    """
#    produce a director from a rinex, and write it in the specific location
#
#    Return :
#            the path of the produced director as a string
#
#    If the rinex is not in a gins style folder, it will be copied in the
#    ``temp_data_folder``
#
#    If ``temp_data_folder`` is not specified ( == '') rinex will be copied in
#    a ad hoc folder ``../gin/TEMP_DATA``
#
#    If ``out_director_folder`` is not specified ( == ''), output director will be
#    created in the ``../gin/data/directeur folder``
#
#    automatic mode : (auto_staocl)
#        create automaticaly a station file and a ocean loading file
#        auto mode is prioritary upon the manu mode
#
#        so, if a path for stat file or ocload file is specified and
#        auto_staocl is on, it will be the automatic stat/ocload
#        which will be used
#    perso_orbclk :
#        download and use specifics orbits
#        according to ac & repro args
#        (they are useless if perso_orbclk aren't activated)
#    auto_interval:
#        find the interval in RINEX and apply it to the director
#    """
#
#    rinex_name = os.path.basename(rinex_path)
#    dt_rinex  = geok.rinexname2dt(rinex_name)
#    gins_path = get_gins_path()
#
#    # be sure there is a TEMP DATA folder
#    if temp_data_folder=='':
#        temp_data_folder = os.path.join(gins_path,'gin','TEMP_DATA')
#    if not os.path.exists(temp_data_folder):
#        os.makedirs(temp_data_folder)
#
#    # be sure the RINEX is in gins folder ...
#    bool_rinex_in_gin = check_if_file_in_gin_folder(rinex_path)
#    # ... and copy it otherwise
#    if not bool_rinex_in_gin:
#        print 'INFO : will be copied in ', temp_data_folder
#        rinex_path = copy_file_in_gin_folder(rinex_path,temp_data_folder)
#    # check if the RINEX is compressed ...
#    bool_comp_rnx = check_if_compressed_rinex(rinex_path)
#    # ... if not crz2rnx !
#    if bool_comp_rnx:
#        rinex_path = crz2rnx(rinex_path,temp_data_folder)
#
#    if out_director_folder == '':
#        out_director_folder = os.path.join(get_gins_path(),'gin','data','directeur')
#
#    stat_lower = rinex_name[0:4]
#    stat = stat_lower.upper()
#
#    if perso_orbclk:
#        calccntr_suffix = '_' + ac + str(repro)
#    else:
#        calccntr_suffix = ''
#
#    # date
#    strt_epoch , end_epoch = operational.rinex_start_end(rinex_path)
#    ses = '_' + operational.rinex_session_id(strt_epoch,end_epoch,full_mode=1)
#
#    # Interval Initalisation
#    # interval_line = utils.grep(rinex_path,'INTERVAL',True)
#    #  float(interval_line.split()[0])
#    _ , _ , freq_rnx_true = operational.rinex_start_end(rinex_path,True)
#    if auto_interval:
#        freq_rnx_str = '_' + str(int(freq_rnx_true)).zfill(2) + 's'
#    else:
#        freq_rnx_str = ''
#
#    director_output_name = director_name_prefix + '_' + stat_lower + '_' + str(geok.dt2jjul_cnes(dt_rinex)) + '_' + str(dt_rinex.year) + '_' + geok.dt2doy(dt_rinex) + freq_rnx_str + ses + calccntr_suffix +'.yml'
#    director_output_path = os.path.join(out_director_folder,director_output_name)
#    director_generik_file = open(director_generik_path)
#    director_dic = yaml.load(director_generik_file)
#
#    # ADDING RINEX NAME & DATES INTO THE DIR
#    # name
#    rnx_path_dir_compatible = make_path_ginsstyle(rinex_path)
#    director_dic['observation']['interobject_data'][1]['file'] = rnx_path_dir_compatible
#    # date
#    director_dic['date']['arc_start'][0] = geok.dt2jjul_cnes(strt_epoch)
#    director_dic['date']['arc_stop'][0] = geok.dt2jjul_cnes(end_epoch)
#    director_dic['date']['initial_state_vector_date'][0] = geok.dt2jjul_cnes(strt_epoch)
#
#    director_dic['date']['arc_start'][1] = strt_epoch.hour * 3600 + strt_epoch.second + 19
#    director_dic['date']['arc_stop'][1]  = end_epoch.hour  * 3600 + end_epoch.second + 19
#    director_dic['date']['initial_state_vector_date'][1] = strt_epoch.second + 19
#
#    # inteval
#    if auto_interval: # dir is modified only if freq_rnx >= 0
#        director_dic['observation']['removal']['simulation_stepsize'] = freq_rnx_true
#        # check if its not a static dir
#        if director_dic['parameter']['adjustment_parameters']['stations']['adjustment_frequency'][-1] != 0:
#            director_dic['parameter']['adjustment_parameters']['stations']['adjustment_frequency'] = [0, 0, 0, freq_rnx_true]
#
#
#    #station file and oceanload file , auto or manu mode
#    # auto mode is prioritary upon the manu mode
#
#    if auto_staocl == False and (oceanload_file == '' or stations_file == ''):
#        print "WARN : generate dir : Auto mode is off and no stat or oclo file specified !!!"
#
#    if auto_staocl:
#        print '* stat and ocload files generation :'
#        randid = 'id'+str(np.random.randint(500,999))
#        stations_file_name = stat_lower + '_' + str(dt_rinex.year) + '_' + \
#        str(geok.dt2doy(dt_rinex)) + '_' + randid + '.stat'
#        stations_file = os.path.join(temp_data_folder, stations_file_name)
#        gfc.write_station_file_gins_from_rinex(rinex_path,stations_file)
#        oceanload_file_name = stat_lower + '_' + str(dt_rinex.year) + '_' + \
#        str(geok.dt2doy(dt_rinex)) + '_' + randid + '.oclo'
#        oceanload_file = os.path.join(temp_data_folder , oceanload_file_name)
#        write_oceanload_file(stations_file,oceanload_file)
#        bool_statfil = os.path.isfile(stations_file)
#        bool_oclofil = os.path.isfile(oceanload_file)
#        iwhile = 0
#        while not (bool_statfil and bool_oclofil) or iwhile < 6:
#            iwhile = iwhile + 1
#            time.sleep(.5)
#            bool_statfil = os.path.isfile(stations_file)
#            bool_oclofil = os.path.isfile(oceanload_file)
#
#        print stations_file, bool_oclofil
#        print oceanload_file , bool_oclofil
#        if not os.path.isfile(oceanload_file):
#            print 'ERR : no oceload file !!!'
#
#    if oceanload_file != '':
#        director_dic['object']['station']['ocean_tide_loading'] = make_path_ginsstyle(oceanload_file)
#    if stations_file != '':
#        director_dic['object']['station']['station_coordinates'] = make_path_ginsstyle(stations_file)
#
#    # =========   ORBITS   =============
#    if not perso_orbclk:
#        # GRG or GR2 for clock and orbits
#        if dt_rinex <= dt.datetime(2013,12,29):
#            gr_gr2 = 'GR2'
#        else:
#            gr_gr2 = 'GRG'
#        horpath = os.path.join('horloges',gr_gr2,'defaut')
#        orbpath = os.path.join('orbites',gr_gr2,'defaut')
#    else:
#        orbpath , horpath = download_convert_2_gins_orb_clk(dt_rinex,temp_data_folder,
#                                        ac=ac,repro=repro)
#        orbpath = make_path_ginsstyle(orbpath)
#        horpath = make_path_ginsstyle(horpath)
#
#    director_dic['model']['environment']['gnss_clock'] = horpath
#    director_dic['observation']['interobject_data'][0]['file'] = orbpath
#
#    # correcting a pyyaml bug : version is interpreted as a int instead of a string
#    verslis = list(str(director_dic['version']))
#
#    verstr = verslis[0] + verslis[1] + '_' + '_'.join(verslis[2:])
#    director_dic['version'] = verstr
#
#    # Security CHECKs
#    statfile_path_ginsstyle = director_dic['object']['station']['station_coordinates']
#    statfile_path_full = os.path.join(gins_path,statfile_path_ginsstyle[6:])
#
#    ocloadfile_path_ginsstyle = director_dic['object']['station']['ocean_tide_loading']
#    ocloadfile_path_full = os.path.join(gins_path,ocloadfile_path_ginsstyle[6:])
#
#    check_if_stat_in_stationfile(stat,statfile_path_full)
#
#    domes = find_DOMESstat_in_stationfile(stat,statfile_path_full)
#    print 'INFO : DOMES : ',domes
#    check_if_DOMES_in_oceanloadfile(domes[0],ocloadfile_path_full)
#
#    # Case of kinematic process : need to change keys in user extension
#    if 'user_extension' in director_dic:
#        if 'userext_gps__haute_freq' in director_dic['user_extension']:
#            director_dic['user_extension']['userext_gps__haute_freq'] = stat
#
#    # WRITING THE NEW DIRECTOR
#    print 'INFO : writing : ', director_output_path
#    with open(director_output_path, 'w+') as outfile:
#        outfile.write( yaml.dump(director_dic,default_flow_style=False) )
#
#    # correcting a bug : GINS YAML interpreter can't manage false or true ...
#    utils.replace(director_output_path,'false','no')
#    utils.replace(director_output_path,'true','yes')
#
#    return director_output_path

# def multi_gen_dir_from_rinex_list(rinex_list,director_generik_path,
#    director_name_prefix,out_director_folder='',temp_data_folder='',
#    stations_file='',oceanload_file='',auto_staocl=False,
#    perso_orbclk=False,ac='igs',repro=2,auto_interval=True):
#
#    multi_gen_dir_from_rinex_list.__doc__ = generate_director_from_rinex.__doc__
#
#    N = len(rinex_list)
#    director_output_path = []
#    print "******** DIRECTORS GENRATION ********"
#    for i,rnx in enumerate(rinex_list):
#        print ' ======== ' , i+1 , '/' , N , ' ======== '
#        print ' === ' , os.path.basename(rnx)
#        out_dir_path = generate_director_from_rinex(rnx,director_generik_path,
#        director_name_prefix,out_director_folder,temp_data_folder,
#        stations_file,oceanload_file,auto_staocl=auto_staocl,
#        perso_orbclk=perso_orbclk,ac=ac,repro=repro,
#        auto_interval=auto_interval)
#
#        director_output_path.append(out_dir_path)
#
#    director_output_path = sorted(director_output_path)
#    return director_output_path
#
#


def prairie_manual(
    rinex_paths_in,
    temp_data_folder="",
    force=True,
    with_historik=True,
    with_wsb=True,
    argsdict=dict(),
):
    """
    argsdict :
        a dictionnary so as argsdict[argument] = val
        e.g.
        argsdict['-options'] = "/home/psakicki/THESE/SOFTWARES/GINSv2B/gin/data/prairie/options_GPS.dat"
        argsdict['-wsb']     = "/home/psakicki/THESE/SOFTWARES/GINSv2B/gin/data/prairie/WSB_NEW_1996.dat"

    argument cc2noncc, wsb, options, and out are automatically managed
    """

    if type(rinex_paths_in) is list:
        multimode = True
        rinex_path_lis = rinex_paths_in
    elif type(rinex_paths_in) is str:
        multimode = False
        rinex_path_lis = [rinex_paths_in]

    # be sure there is a TEMP DATA folder
    if temp_data_folder == "":
        temp_data_folder = os.path.join(get_gins_path(1), "TEMP_DATA")
    if not os.path.exists(temp_data_folder):
        os.makedirs(temp_data_folder)

    out_pra_path_lis = []
    for rinex_path in rinex_path_lis:
        #        # be sure the RINEX is in gins folder ...
        #        bool_rinex_in_gin = check_if_file_in_gin_folder(rinex_path)
        #        # ... and copy it otherwise
        #        if not bool_rinex_in_gin:
        expected_rinex_path = os.path.join(
            temp_data_folder, os.path.basename(rinex_path)
        )
        if not os.path.isfile(expected_rinex_path):
            print("INFO : will be copied in ", temp_data_folder)
            rinex_path = copy_file_in_gin_folder(rinex_path, temp_data_folder)
        else:
            print("INFO : ", expected_rinex_path, "already exists")
            rinex_path = expected_rinex_path
        # check if the RINEX is compressed ...
        bool_comp_rnx = operational.check_if_compressed_rinex(rinex_path)
        # ... if not crz2rnx !
        if bool_comp_rnx:
            crinex_path = rinex_path
            rinex_path = operational.crz2rnx(crinex_path, temp_data_folder)
            if not os.path.isfile(rinex_path):
                print("ERR : something went wrong during CRZ2RNX, skiping the file")
                print(crinex_path)

        rinex_name = os.path.basename(rinex_path)
        out_pra_name = rinex_name + ".pra"
        print("INFO : input & output name :", rinex_name, out_pra_name)
        out_pra_path = os.path.join(temp_data_folder, out_pra_name)

        if os.path.isfile(out_pra_path) and not force:
            print(
                "INFO : prairie_manual : ", out_pra_path, "already exists, skiping ..."
            )

        command = "exe_prairie "

        if force:
            forcearg = " -f "
        else:
            forcearg = ""

        command = command + " " + rinex_name + forcearg

        argsdict["-out"] = out_pra_name

        if not "-cc2noncc" in list(argsdict.keys()):
            argsdict["-cc2noncc"] = os.path.join(
                get_gins_path(1), "archives/p1c1bias.2000p"
            )
        if with_historik and not "-historik" in list(argsdict.keys()):
            argsdict["-historik"] = os.path.join(
                get_gins_path(1), "data/constell/historik_glonass"
            )
        if with_wsb and not "-wsb" in list(argsdict.keys()):
            argsdict["-wsb"] = os.path.join(get_gins_path(1), "archives/WSBREF.res.dat")
        if not "-options" in list(argsdict.keys()):
            argsdict["-options"] = os.path.join(
                get_gins_path(1), "data/prairie/options_GPS.dat"
            )

        for k, v in argsdict.items():
            if "-" == k[0]:
                moins = ""
            else:
                moins == "-"

            command = command + " " + moins + k + " " + v

        #        stream = os.popen(command)
        curdir = os.getcwd()
        os.chdir(temp_data_folder)
        p = subprocess.Popen(
            [command], shell=True
        )  # , stdout=subprocess.PIPE , stderr=subprocess.PIPE)
        p.wait()
        os.chdir(curdir)

        print("")
        print(command)
        print("")
        #
        if os.path.isfile(out_pra_path) and os.stat(out_pra_path).st_size > 0:
            print("INFO : prairie_manual : OK for manual prairie :)")
            print("output in ", out_pra_path)
            out_pra_path_lis.append(out_pra_path)
        else:
            print("ERR : prairie_manual : ", out_pra_path, "dont exist or empty :(")

    if len(out_pra_path_lis) == 1:
        return out_pra_path
    else:
        return out_pra_path_lis


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
        print("******** DIRECTORS GENRATION ********")
    elif type(rinex_paths_in) is str:
        multimode = False
        rinex_path_lis = [rinex_paths_in]
    else:
        print("ERR : gen_dirs_from_rnxs : check the rinex_paths_in !!!")
        return None

    N = len(rinex_path_lis)
    director_output_path_lis = []
    failed_rinex_date_lis = []

    for i, rinex_path in enumerate(rinex_path_lis):
        if multimode:
            print(" ======== ", i + 1, "/", N, " ======== ")
            print(" === ", os.path.basename(rinex_path))
        rinex_name = os.path.basename(rinex_path)
        dt_rinex = conv.rinexname2dt(rinex_name)
        gins_path = get_gins_path()

        # be sure there is a TEMP DATA folder
        if temp_data_folder == "":
            temp_data_folder = os.path.join(gins_path, "gin", "TEMP_DATA")
        if not os.path.exists(temp_data_folder):
            os.makedirs(temp_data_folder)

        # be sure the RINEX is in gins folder ...
        bool_rinex_in_gin = check_if_file_in_gin_folder(rinex_path)
        # ... and copy it otherwise
        if not bool_rinex_in_gin:
            print("INFO : will be copied in ", temp_data_folder)
            rinex_path = copy_file_in_gin_folder(rinex_path, temp_data_folder)
        # check if the RINEX is compressed ...
        bool_comp_rnx = operational.check_if_compressed_rinex(rinex_path)
        # ... if not crz2rnx !
        if bool_comp_rnx:
            crinex_path = rinex_path
            rinex_path = operational.crz2rnx(crinex_path, temp_data_folder)
            if not os.path.isfile(rinex_path):
                print("ERR : something went wrong during CRZ2RNX, skiping the file")
                print(crinex_path)

        # prairie extern
        if prairie:
            print("INFO : run prairie externally")
            pra_file_path = prairie_manual(
                rinex_path, temp_data_folder, **prairie_kwargs
            )

            if type(pra_file_path) is list and len(pra_file_path) == 0:
                print(
                    "ERR : something went wrong during prairie externally, skip the day"
                )
                print(rinex_path)
                failed_rinex_date_lis.append(dt_rinex)
                if multimode:
                    continue
                else:
                    return ""

        if out_director_folder == "":
            out_director_folder = os.path.join(
                get_gins_path(), "gin", "data", "directeur"
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
            print(
                "ERR : gen_dirs_from_rnxs : unable to get RINEX start/end, skiping ...",
                end=" ",
            )
            print(rinex_path)
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

        director_output_name = (
            director_name_prefix
            + "_"
            + stat_lower
            + "_"
            + str(geok.dt2jjul_cnes(dt_rinex))
            + "_"
            + str(dt_rinex.year)
            + "_"
            + geok.dt2doy(dt_rinex)
            + freq_rnx_str
            + coord_prefix
            + ses
            + calccntr_suffix
            + ".yml"
        )
        director_output_path = os.path.join(out_director_folder, director_output_name)
        director_generik_file = open(director_generik_path)
        director_dic = yaml.load(director_generik_file)

        # ADDING RINEX NAME & DATES INTO THE DIR
        # name
        if not prairie:
            rnx_path_dir_compatible = make_path_ginsstyle(rinex_path)
            director_dic["observation"]["interobject_data"][1][
                "file"
            ] = rnx_path_dir_compatible
        else:
            pra_path_dir_compatible = make_path_ginsstyle(pra_file_path)
            director_dic["observation"]["interobject_data"][1][
                "file"
            ] = pra_path_dir_compatible

        # date
        strt_day, strt_sec = geok.dt2jjul_cnes(strt_epoch, False)
        end_day, end_sec = geok.dt2jjul_cnes(end_epoch, False)

        director_dic["date"]["arc_start"][0] = strt_day
        director_dic["date"]["arc_stop"][0] = end_day
        if "initial_state_vector_date" in list(director_dic["date"].keys()):
            director_dic["date"]["initial_state_vector_date"][0] = strt_day

        director_dic["date"]["arc_start"][1] = strt_sec
        director_dic["date"]["arc_stop"][1] = end_sec
        if "initial_state_vector_date" in list(director_dic["date"].keys()):
            director_dic["date"]["initial_state_vector_date"][1] = strt_sec

        # interval
        if auto_interval:  # dir is modified only if freq_rnx >= 0
            director_dic["observation"]["removal"][
                "simulation_stepsize"
            ] = freq_rnx_true
            # check if its not a static dir
            if (
                director_dic["parameter"]["adjustment_parameters"]["stations"][
                    "adjustment_frequency"
                ][-1]
                != 0
            ):
                director_dic["parameter"]["adjustment_parameters"]["stations"][
                    "adjustment_frequency"
                ] = [0, 0, 0, freq_rnx_true]
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

        # station file and oceanload file , auto or manu mode
        # auto mode is prioritary upon the manu mode

        if auto_staocl == False and (oceanload_file == "" or stations_file == ""):
            print(
                "WARN : generate dir : Auto mode is off and no stat or oclo file specified !!!"
            )
        if auto_staocl:
            print("* stat and ocload files generation :")
            randid = "id" + str(np.random.randint(10000, 99999))
            stations_file_name = (
                stat_lower
                + "_"
                + str(dt_rinex.year)
                + "_"
                + str(geok.dt2doy(dt_rinex))
                + "_"
                + randid
                + ".stat"
            )
            stations_file = os.path.join(temp_data_folder, stations_file_name)
            gfc.write_station_file_gins_from_rinex(rinex_path, stations_file)
            oceanload_file_name = (
                stat_lower
                + "_"
                + str(dt_rinex.year)
                + "_"
                + str(geok.dt2doy(dt_rinex))
                + "_"
                + randid
                + ".oclo"
            )
            oceanload_file = os.path.join(temp_data_folder, oceanload_file_name)
            write_oceanload_file(stations_file, oceanload_file)

            bool_statfil = os.path.isfile(stations_file)
            bool_oclofil = os.path.isfile(oceanload_file)
            iwhile = 0
            while not (bool_statfil and bool_oclofil) or iwhile < 6:
                if iwhile > 6:
                    print(
                        "WARN : something is going wrong in the loop waiting for the ocean loading file ...",
                        iwhile,
                    )
                iwhile = iwhile + 1
                time.sleep(0.5)
                bool_statfil = os.path.isfile(stations_file)
                bool_oclofil = os.path.isfile(oceanload_file)

            print(stations_file, bool_oclofil)
            print(oceanload_file, bool_oclofil)
            if not os.path.isfile(oceanload_file):
                print("ERR : no oceload file !!!")

        if oceanload_file != "":
            director_dic["object"]["station"]["ocean_tide_loading"] = (
                make_path_ginsstyle(oceanload_file)
            )
        if stations_file != "":
            director_dic["object"]["station"]["station_coordinates"] = (
                make_path_ginsstyle(stations_file)
            )

        print("AAA3")

        # =========   ORBITS   =============
        if not perso_orbclk:
            # GRG or GR2 for clock and orbits
            if dt.datetime(1998, 1, 4) <= dt_rinex <= dt.datetime(2013, 12, 29):
                gr_gr2 = "GR2"
            #    elif dt_rinex < dt.datetime(1998,1,4):
            #        gr_gr2 = 'GR1' # obsolete parce que GR1 intégralement inclu dans GR2
            elif dt_rinex < dt.datetime(1998, 1, 4):
                gr_gr2 = "IG2"
                print("INFO : no GRG/GR2 orbit available for day", dt_rinex)
                print("       (day before 1998/1/4) => using IG2 orb. instead")
            else:
                gr_gr2 = "GRG"
            horpath = os.path.join("horloges", gr_gr2, "defaut")
            orbpath = os.path.join("orbites", gr_gr2, "defaut")
        else:
            orbpath, horpath = download_convert_2_gins_orb_clk(
                dt_rinex, temp_data_folder, ac=ac, repro=repro
            )

            # Exception case where no sp3/clk was found
            if orbpath == None or horpath == None:
                print(
                    "WARN : something went wrong with orbits/clocks download/conversion ",
                    dt_rinex,
                )
                print("       skiping the day ...")
                failed_rinex_date_lis.append(dt_rinex)
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

        check_if_stat_in_stationfile(stat, statfile_path_full)

        domes = find_DOMESstat_in_stationfile(stat, statfile_path_full)
        print("INFO : DOMES : ", domes)
        check_if_DOMES_in_oceanloadfile(domes[0], ocloadfile_path_full)

        # Case of kinematic process : need to change keys in user extension
        if "user_extension" in director_dic:
            print("INFO : kinematic process with 'user_extension' fields")
            if "userext_gps__haute_freq" in director_dic["user_extension"]:
                director_dic["user_extension"]["userext_gps__haute_freq"] = stat
            if "userext_gps__haute_freq".upper() in director_dic["user_extension"]:
                director_dic["user_extension"]["userext_gps__haute_freq".upper()] = stat

            if "userext_addition" in director_dic["user_extension"]:
                if (
                    "gps__haute_freq"
                    in director_dic["user_extension"]["userext_addition"]
                ):
                    director_dic["user_extension"]["userext_addition"][
                        "gps__haute_freq"
                    ] = stat
                if (
                    "gps__haute_freq".upper()
                    in director_dic["user_extension"]["userext_addition"]
                ):
                    director_dic["user_extension"]["userext_addition"][
                        "gps__haute_freq".upper()
                    ] = stat
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
                    "gps__haute_freq ".upper() + stat
                )

        # WRITING THE NEW DIRECTOR
        print("INFO : writing : ", director_output_path)
        with open(director_output_path, "w+") as outfile:
            outfile.write(yaml.dump(director_dic, default_flow_style=False))

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
                    "WARN : something went wrong with orbits/clocks download/conversion ",
                    dt_rinex,
                )
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

        # check_if_stat_in_stationfile(stat,statfile_path_full)

        # domes = find_DOMESstat_in_stationfile(stat,statfile_path_full)
        # print 'INFO : DOMES : ',domes
        # check_if_DOMES_in_oceanloadfile(domes[0],ocloadfile_path_full)

        # Case of kinematic process : need to change keys in user extension
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


def get_temp_data_gins_path():

    gins_path = get_gins_path()

    # be sure there is a TEMP DATA folder
    temp_data_folder = os.path.join(gins_path, "gin", "TEMP_DATA")

    if not os.path.exists(temp_data_folder):
        os.makedirs(temp_data_folder)

    return temp_data_folder


def dirs_copy_generik_2_working(
    director_name_prefix, director_generik_path, temp_data_folder=None
):

    if not temp_data_folder:
        temp_data_folder = get_temp_data_gins_path()

    shutil.copy(director_generik_path, temp_data_folder)
    director_generik_path_tmp = os.path.join(
        temp_data_folder, os.path.basename(director_generik_path)
    )
    director_generik_path_out = os.path.join(
        temp_data_folder, director_name_prefix + ".yml"
    )
    os.rename(director_generik_path_tmp, director_generik_path_out)

    return director_generik_path_out


# def run_director(director_path,opts_gins_pc='',opts_gins_90=''):
#    """ run a simple director """
#    director_name = os.path.basename(director_path)
#    opts_gins_pc = '-F' + opts_gins_pc
#    start = time.time()
#
#    print 'INFO : options ginsPC / gins90 :',opts_gins_pc,'/',opts_gins_90
#
#    if 'IPPP' in opts_gins_90:
#        for grepstr in ('userext_gps__qualiteorb','userext_gps__haute_freq',
#        'userext_gps__hor_interp'):
#            grep_out = utils.grep(director_path,grepstr)
#            if grep_out == '':
#                print "WARN : IPPP mode on, but no",grepstr," in the dir !!!"
#
#    command = "ginspc.bash " + opts_gins_pc + ' ' + director_name + ' ' + opts_gins_90 + ' -v OPERA'
#
#    # l'argument OPERA est super important  !!!
#    # c'est lui qui resoult l'instabilité lors du lancement du subprocess !!!
#    # parce que indirectement exe_gins ne marche pas sans
#    # faire le test avec un exe_gins -v OPERA -fic <fic> et sans
#
#    print 'INFO : submit. command : ' , command
#
#    stream = os.popen(command)
#    stream_str = stream.read()
#    check_good_exec_of_GINS(stream_str,director_name)
#    gins_path = get_gins_path()
#    log_path = os.path.join(gins_path,'python_logs',director_name + ".log")
#
#    if not os.path.exists(os.path.dirname(log_path)):
#        os.makedirs(os.path.dirname(log_path))
#
#    with open(log_path + '.exec', "w") as f:
#        f.write('exec time : ' + str(time.time() - start))
#    with open(log_path , "w") as f:
#        f.write(stream_str)


# def run_director_list(director_path_lis,opts_gins_pc='',opts_gins_90=''):
#    """ run a list of dir in series (on the same 'slot') """
#    N = len(director_path_lis)
#    print "******** DIRECTORS RUNS *******"
#    for i,dirr in enumerate(director_path_lis):
#        print ' ======== ' , i+1 , '/' , N , ' ======== '
#        print 'INFO : launching : ', dirr
#        if dirr[-1] == '~':
#            print 'INFO : geany ~ temp file, skiping this dir. '
#            continue
#        run_director(dirr,opts_gins_pc,opts_gins_90)


def run_directors(
    dir_paths_in, opts_gins_pc="", opts_gins_90="  ", version="OPERA", fic_mode=False
):
    """
    NEW FCT WHICH CAN MANAGE BOTH ONE RINEX OR A LIST OF RINEX, OR A FIC file (170613)
    """
    # Multi or Single Mode ?
    if type(dir_paths_in) is list:
        multimode = True
        director_path_lis = dir_paths_in
        print("******** DIRECTORS RUNS *******")
    elif type(dir_paths_in) is str:
        multimode = False
        director_path_lis = [dir_paths_in]
    else:
        print("ERR : run_directors : check the rinex_paths_in !!!")
        return None

    N = len(director_path_lis)

    for i, director_path in enumerate(director_path_lis):
        if multimode:
            print(" ======== ", i + 1, "/", N, " ======== ")
        print("INFO : launching : ", director_path)
        print("INFO : start at", dt.datetime.now())
        if director_path[-4:] == ".fic":
            fic_mode = True
            print("INFO : input file ends with .fic, fic_mode is activated")
        start = time.time()

        if director_path[-1] == "~":
            print("INFO : geany ~ temp file, skiping this dir. ")
            continue

        director_name = os.path.basename(director_path)
        opts_gins_pc_ope = "-F" + opts_gins_pc

        print("INFO : options ginsPC / gins90 :", opts_gins_pc_ope, "/", opts_gins_90)

        if "IPPP" in opts_gins_90 and not fic_mode:
            for grepstr in (
                "userext_gps__qualiteorb",
                "userext_gps__haute_freq",
                "userext_gps__hor_interp",
                "GPS__QUALITEORB",
                "GPS__HAUTE_FREQ",
                "GPS__HOR_INTERP",
            ):
                grep_out = utils.grep(director_path, grepstr)
                if grep_out == "":
                    print("WARN : IPPP mode on, but no", grepstr, " in the dir !!!")

        if not fic_mode:
            command = (
                "ginspc.bash "
                + opts_gins_pc_ope
                + " "
                + director_name
                + " "
                + opts_gins_90
                + " -v "
                + version
            )
        else:
            command = "exe_gins " + " -fic " + director_name + " -v " + version
        # l'argument OPERA est super important  !!!
        # c'est lui qui resoult l'instabilité lors du lancement du subprocess !!!
        # parce que indirectement exe_gins ne marche pas sans
        # faire le test avec un exe_gins -v OPERA -fic <fic> et sans

        print("INFO : submit. command : ", command)

        gins_path = get_gins_path()
        log_path = utils.create_dir(os.path.join(gins_path, "python_logs"))
        log_path = os.path.join(gins_path, "python_logs", director_name + ".log")

        with open(log_path, "w+") as f:
            process = subprocess.Popen(
                [command],
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                executable="/bin/bash",
            )
            for c in iter(lambda: process.stdout.read(1), ""):
                # https://stackoverflow.com/questions/436220/determine-the-encoding-of-text-in-python4
                # find encoding
                # print(c)
                dammit = UnicodeDammit(c)
                # print(dammit.original_encoding)
                d_deco = c.decode("iso-8859-1", errors="replace")
                sys.stdout.write(d_deco)
                f.write(d_deco)

        check_good_exec_of_GINS(open(log_path), director_name)

        with open(log_path + ".exec", "w") as f:
            f.write("exec time : " + str(time.time() - start))

        print("INFO : end at", dt.datetime.now())
        print("INFO : exec time : ", str(time.time() - start), "seconds")

    # Vieux lancement avec popen mais ca marche plus (170207)
    #        stream = os.popen(command)
    #        stream_str = stream.read()
    #        check_good_exec_of_GINS(stream_str,director_name)
    #        gins_path = get_gins_path()
    #        log_path = os.path.join(gins_path,'python_logs',director_name + ".log")
    #
    #        if not os.path.exists(os.path.dirname(log_path)):
    #            os.makedirs(os.path.dirname(log_path))
    #
    #        with open(log_path + '.exec', "w") as f:
    #            f.write('exec time : ' + str(time.time() - start))
    #        with open(log_path , "w") as f:
    #            f.write(stream_str)
    #        print 'INFO : end at' , dt.datetime.now()
    #        #print 'INFO : exec time : ' , str(time.time() - start))

    return None


def run_director_wrap(intup):
    run_directors(*intup)
    return None


def run_director_list_wrap(tupinp):
    run_directors(*tupinp)


def run_dirs_multislots(
    director_lis,
    slots_lis=["", "U", "L", "R"],
    opts_gins_pc="",
    opts_gins_90="",
    version="OPERA",
    fic_mode=False,
):
    """run a list of dir in parallel (using different 'slots')"""
    """ FRONTEND FUNCTION TO USE """
    if not type(director_lis) is list:
        print("ERR : run_dirs_multislots : director_lis in input is not a list !!!")
        return None
    print("TOTAL NB OF DIRECTORS :", len(director_lis))

    if type(slots_lis) is int:
        slots_lis = [""] * slots_lis

    chunk = utils.chunkIt(director_lis, len(slots_lis))
    pool_size = len(slots_lis)
    pool = mp.Pool(processes=pool_size)
    ZIP = list(
        zip(
            chunk,
            [slot + opts_gins_pc for slot in slots_lis],
            [opts_gins_90] * len(slots_lis),
            [version] * len(slots_lis),
            [fic_mode] * len(slots_lis),
        )
    )
    pool.map(run_director_list_wrap, ZIP)
    return None


def get_director_list(wildcard_dir):
    """with a wildcard (e.g. 'GWADA_MK2*') and return a list of corresponding
    directors found in gin/data/directeur folder"""
    gins_path = get_gins_path()
    di_run_lis = [
        os.path.basename(e)
        for e in glob.glob(os.path.join(gins_path, "data", "directeur", wildcard_dir))
    ]
    return di_run_lis


def get_rinex_list(
    parent_folder,
    specific_stats=[],
    invert=False,
    compressed=True,
    start=dt.datetime(1980, 1, 1),
    end=dt.datetime(2099, 1, 1),
):
    """return :
        all RINEXs found in a parent folder and his subfolders
        (compressed or not)

    parent_folder :
        can be a the path of the partent folder (the rinex archive path)
        but also a RINEX list already found
        (modification of 170524 to gain speed)

    specific_stats :
        MUST BE a list ['STA1','STA2','STA3']
        So, if only one elt in specific_stats, use a syntax as : ['STA1']

    invert :
        False = keeping the specific stats
        True  = removing the specific stats
        or for all stations, leave a empty
        tuple in specific_stats

    NB : the end date is included
    """
    if type(parent_folder) is str:
        wholefilelist = []
        for root, dirs, files in os.walk(parent_folder, topdown=False):
            for name in files:
                wholefilelist.append(os.path.join(root, name))

        wholefilelist = list(set(wholefilelist))

        if compressed:
            rnxregex = operational.rinex_regex()
        else:
            rnxregex = operational.rinex_regex(False)

        rinexfilelist = [fil for fil in wholefilelist if re.search(rnxregex, fil)]
    elif utils.is_iterable(parent_folder):
        # is parent_folder is a list containing already found RINEXs
        rinexfilelist = parent_folder

    specific_stats = [e.lower() for e in specific_stats]

    if specific_stats != []:
        if not invert:
            goodrnxfilelist = [
                fil
                for fil in rinexfilelist
                if os.path.basename(fil)[0:4] in specific_stats
            ]
        else:
            goodrnxfilelist = [
                fil
                for fil in rinexfilelist
                if os.path.basename(fil)[0:4] not in specific_stats
            ]
    else:
        goodrnxfilelist = rinexfilelist

    goodrnxfilelist = sorted(goodrnxfilelist)

    if goodrnxfilelist == []:
        print(
            "WARN : get_rinex_list : no RINEX found, check the path of the parent folder !"
        )
    else:
        print("INFO : get_rinex_list :", len(goodrnxfilelist), "RINEXs found")

    goodrnxfilelist_date = []
    end = end + dt.timedelta(seconds=86399)
    for rnx in goodrnxfilelist:
        dtrnx = conv.rinexname2dt(os.path.basename(rnx))
        if start <= dtrnx <= end:
            goodrnxfilelist_date.append(rnx)
    if len(goodrnxfilelist_date) != len(goodrnxfilelist):
        dif = len(goodrnxfilelist) - len(goodrnxfilelist_date)
        print(
            "INFO : get_rinex_list :",
            dif,
            "RINEXs removed bc. not in the date interval",
        )

    if goodrnxfilelist_date == []:
        print("WARN : get_rinex_list : all RINEXs removed bc. of date interval")

    return goodrnxfilelist_date


def smart_directors_to_run(wildcard_dir="", full_path_out=True):
    """smart runner check directors who worked, and thus give a list of directors
    whitout those who worked

    listing and directeur folders are inspected automatically"""

    gins_path = get_gins_path(True)

    # list of corresponding directors
    dirpath = os.path.join(gins_path, "data", "directeur", wildcard_dir)
    dirlis = glob.glob(dirpath)
    dirlis = [os.path.basename(e) for e in dirlis]

    # list of corresponding already lauched listings
    lipath = os.path.join(gins_path, "batch", "listing", wildcard_dir)
    lilis = glob.glob(lipath)

    lilis = [os.path.basename(e) for e in lilis]
    # must treat new and old independently
    lilis_new = [e for e in lilis if ".yml" in e]
    lilis_old = [e for e in lilis if not ".yml" in e]
    # 1) treating OLD case
    # MAKING LISTS OF TUPLES (PREFIX , <gins or prepars>)
    li_tuple_lis_old = [(e.split(".")[0], e.split(".")[-1]) for e in lilis_old]
    li_finished_lis_old = [e[0] for e in li_tuple_lis_old if "gins" in e[1]]

    # 2) treating NEW case
    # MAKING LISTS OF TUPLES (PREFIX , <gins or prepars>)
    li_tuple_lis_new = [
        (e.split(".")[0] + "." + e.split(".")[1], e.split(".")[-1]) for e in lilis_new
    ]
    li_finished_lis_new = [e[0] for e in li_tuple_lis_new if "gins" in e[1]]

    li_finished_lis = li_finished_lis_old + li_finished_lis_new

    # director to run are those in the dir list minus those which worked
    di_run_lis = list(set(dirlis) - set(li_finished_lis))

    print("for wildcard :", wildcard_dir)
    print("input directors found in the data/dir. folder     :", len(dirlis))
    print("FINISHED listings found in the batch/list. folder :", len(li_finished_lis))
    print("directors who need to be launched                 :", len(di_run_lis))

    di_run_lis = sorted(di_run_lis)

    if full_path_out:
        dirpath = os.path.join(get_gins_path(True), "data", "directeur")
        di_run_lis = [os.path.join(dirpath, d) for d in di_run_lis]

    return di_run_lis


def smart_listing_archive(
    wildcard_dir,
    gins_main_archive,
    gins_anex_archive,
    prepars_archive,
    director_archive,
):
    """for each listing corresponding to the wildcard :
    if it's a prepars => go to the prepars_archive
    if it's a gins without duplicate  => go to the gins_main_archive
    if it's a gins with duplcates => one goes to gins_main_archive
                                     the others in gins_anex_archive"""

    listing_path = os.path.join(
        get_gins_path(), "gin", "batch", "listing", wildcard_dir
    )
    listing_lis = [e for e in glob.glob(listing_path)]

    prepars_lis = [e for e in listing_lis if "prepars" in e]
    gins_lis = [e for e in listing_lis if "gins" in e]

    dirpath = os.path.join(get_gins_path(True), "data", "directeur", wildcard_dir)
    dirlis = [e for e in glob.glob(dirpath)]

    # searching for doublons
    # 2 runs of the same day/stat => keeping the first one
    # Only new case (.yml) is treated ...

    gins_prefix_lis = [os.path.basename(e).split(".")[0] for e in gins_lis]
    #    gins_tuple_lis  = [(e.split('.')[0] , e.split('.')[2] ,e.split('.')[3]) for e in gins_lis]

    count_gins_prefix = collections.Counter(gins_prefix_lis).most_common()
    multi_gins_prefix = [elt for elt, count in count_gins_prefix if count > 1]
    simpl_gins_prefix = [elt for elt, count in count_gins_prefix if count == 1]

    simpl_gins_fullpath = []
    for s in simpl_gins_prefix:
        simpl_gins_fullpath = simpl_gins_fullpath + [e for e in gins_lis if s in e]

    multi_gins_fullpath_main = []
    multi_gins_fullpath_anex = []

    for m in multi_gins_prefix:
        multi_gins_fullpath_temp = [e for e in gins_lis if m in e]
        multi_gins_fullpath_main.append(multi_gins_fullpath_temp[0])
        multi_gins_fullpath_anex = (
            multi_gins_fullpath_anex + multi_gins_fullpath_temp[1:]
        )

    for directory in [
        gins_main_archive,
        gins_anex_archive,
        prepars_archive,
        director_archive,
    ]:
        if not os.path.exists(directory):
            os.makedirs(directory)

    # MOVING
    for p in prepars_lis:
        shutil.move(p, prepars_archive)
    for g in simpl_gins_fullpath:
        shutil.move(g, gins_main_archive)
    for g in multi_gins_fullpath_main:
        shutil.move(g, gins_main_archive)
    for g in multi_gins_fullpath_anex:
        shutil.move(g, gins_anex_archive)
    for d in dirlis:
        shutil.move(d, director_archive)

    return None


def sort_by_stations(archive_path, wildcard, i):
    """i is the indice of the first character of the station name in
    eg : i = 18 for filename KARIB_MK3_FLH_v2__bara_22282_2011_003

    archive path is ABSOLUTE b.c. it can be outside of the gins folder"""

    path = os.path.join(archive_path, wildcard)
    filis = glob.glob(path)

    firstfil = os.path.basename(filis[0])
    print("sure the stat name is good ? ")
    print("e.g.", firstfil[i : i + 4], "for", firstfil)
    print("(you have 5sec to abort if not)")
    time.sleep(5)

    for f in filis:
        stat = os.path.basename(f)[i : i + 4]
        STAT = stat.upper()
        archiv = os.path.join(archive_path, STAT)
        if not os.path.exists(archiv):
            os.makedirs(archiv)
        shutil.move(f, archiv)

    return None


# PERSONAL ORBITS


def sp3_2_gins(sp3_pathin):
    os.chdir(os.path.dirname(sp3_pathin))
    # print("WARN : the path of the GINS conversion tool has been hardcoded !!!")
    kommand_1 = get_gins_path() + "/gins_toolbox/scripts/"
    log.info(kommand_1)
    kommand = kommand_1 + "sp3-gins " + sp3_pathin
    #    stream = os.popen(kommand) # forrtl: severe if we use os.popen
    #    subprocess.call([kommand], shell=True)
    log.info(kommand)
    p = subprocess.Popen(
        [kommand], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = p.communicate()
    return sp3_pathin + ".gin"


def clk_2_gins(clk_pathin):
    os.chdir(os.path.dirname(clk_pathin))
    # print("WARN : the path of the GINS conversion tool has been hardcoded !!!")
    kommand_1 = get_gins_path() + "/gins_toolbox/scripts/"
    log.info(kommand_1)
    kommand = kommand_1 + "clk2gins.sh " + os.path.basename(clk_pathin)
    # kommand = kommand_1 + "clk2gins.sh " + clk_pathin
    #    stream = os.popen(kommand) # forrtl: severe if we use os.popen
    #    subprocess.call([kommand], shell=True)
    log.info(kommand)
    p = subprocess.Popen(
        [kommand], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = p.communicate()
    return clk_pathin + ".gins"


def sort_orbit_gins(pathin, pathout):
    kommand = "sort -n -k1,2 -k3 " + pathin + " > " + pathout
    #    stream = os.popen(kommand) # forrtl: severe if we use os.popen
    #    subprocess.call([kommand], shell=True)
    p = subprocess.Popen(
        [kommand], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = p.communicate()
    return pathout


def bad_sat_finder(orbfilein, egrep_ready=True):
    satdic = dict()
    for l in open(orbfilein):
        f = l.split()
        satid = int(f[0])
        if satid not in satdic:
            satdic[satid] = []
        D = conv.jjul_cnes2dt(int(f[1])) + dt.timedelta(seconds=float(f[2]))
        satdic[satid].append(D)

    epochdic = dict()
    epochlis = []
    for k, v in satdic.items():
        epochdic[k] = len(v)
        epochlis.append(len(v))

    nominal_epoch = utils.most_common(epochlis)

    bad_sat_lis = []

    for k, v in epochdic.items():
        if v != nominal_epoch:
            bad_sat_lis.append(k)

    if egrep_ready:
        bad_sat_lis = ["^  " + str(e) for e in bad_sat_lis]

    return bad_sat_lis


def orbit_cleaner(orbfilein, orbfileout):
    """remove sat of a gins orbit file without a nominal number of epoch"""
    bad_sat_lis = bad_sat_finder(orbfilein)
    goodsatslinlis = utils.grep(orbfilein, bad_sat_lis, invert=1, regex=1)

    filobj = open(orbfileout, "w+")

    for lin in goodsatslinlis:
        filobj.write(lin)

    filobj.close()

    return bad_sat_lis


def download_convert_2_gins_orb_clk(
    centraldate,
    work_folder=None,
    ac="jpl",
    repro=2,
    rm_temp_files=True,
    data_center="cddis",
    force=False,
):

    temp_dir = work_folder

    if not temp_dir:
        temp_dir = os.path.join(get_gins_path(), "gin", "TEMP_DATA")
        log.info("use default forder %s", temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    strtdt = centraldate - dt.timedelta(days=1)
    enddt = centraldate + dt.timedelta(days=1)

    if centraldate > dt.datetime(2013, 12, 29) and repro != 0:
        log.warn(str(centraldate) + " a bit late for repro " + str(repro))

    # defining names & paths
    if repro != 0:
        ac_prefix = ac[0:2] + str(repro)  # os.path.basename(gins_orb_lis[0])[0:3]

    else:
        ac_prefix = ac

    week, dow = conv.dt2gpstime(centraldate)
    week, dow = str(week), str(dow)
    centdate = week + dow + centraldate.strftime("_%Y_%j")
    catoutorb = ac_prefix + centdate + ".3d.orb.gin"
    sortoutorb = ac_prefix + centdate + ".3d.orb.sort.gin"
    cleanoutorb = ac_prefix + centdate + ".3d.orb.clean.gin"

    catoutclk = ac_prefix + centdate + ".3d.clk.gin"

    catoutorb_path = os.path.join(temp_dir, catoutorb)
    sortoutorb_path = os.path.join(temp_dir, sortoutorb)
    cleanoutorb_path = os.path.join(temp_dir, cleanoutorb)
    catoutclk_path = os.path.join(temp_dir, catoutclk)

    # preliminary check, if the orbits already exist
    if (
        os.path.isfile(cleanoutorb_path)
        and os.path.isfile(catoutclk_path)
        and not force
    ):
        print("INFO : converted orbits & clocks already exist, it's nice")
        return cleanoutorb_path, catoutclk_path

    # downloading orbits
    sp3Zlis = operational.multi_downloader_orbs_clks_2(
        temp_dir,
        strtdt,
        enddt,
        prod_types=("sp3",),
        AC_names=(ac,),
        repro=repro,
        archtype="/",
        sorted_mode=0,
        data_center=data_center,
    )
    # sorted_mode must be off !!!
    # elsewhere the cat of the orbits/clock is bad !!!

    # downloading clocks
    # 30sec clock prioritary
    clkZlis = operational.multi_downloader_orbs_clks_2(
        temp_dir,
        strtdt,
        enddt,
        prod_types=("clk",),
        AC_names=(ac,),
        repro=repro,
        archtype="/",
        sorted_mode=0,
        data_center=data_center,
    )
    # sorted_mode must be off !!!
    # elsewhere the cat of the orbits/clock is bad !!!

    # standard clock else
    if clkZlis == []:
        print("INFO : clock download : no 30s clock found, trying std. clocks")
        clkZlis = operational.multi_downloader_orbs_clks(
            temp_dir,
            strtdt,
            enddt,
            sp3clk="clk",
            ac=ac,
            data_center=data_center,
            repro=repro,
            archtype="/",
            sorted_mode=0,
        )
        # sorted_mode must be off !!!
        # elsewhere the cat of the orbits/clock is bad !!!

    # Exception case where no files are found
    if len(sp3Zlis) < 3 or len(clkZlis) < 3:
        print(
            "ERR : download_convert_2_gins_orb_clk : missing sp3/clk for", centraldate
        )
        print("list of sp3 : ", sp3Zlis)
        print("list of clk : ", clkZlis)
        return None, None

    # uncompressing orbits
    sp3lis = []
    for Zfil in sp3Zlis:
        if Zfil.endswith(".Z") or Zfil.endswith(".gz"):
            uncomp_fil = files_rw.unzip_gz_Z(Zfil)
        else:
            uncomp_fil = Zfil
        # some SP3 haven't EOF at the end ...
        if not utils.grep_boolean(uncomp_fil, "EOF"):
            with open(uncomp_fil, "a") as myfile:
                print("INFO : uncompress SP3 : adding a EOF line")
                myfile.write("EOF")
                myfile.close()

        sp3lis.append(uncomp_fil)

    # uncompressing clocks
    clklis = []
    for Zfil in clkZlis:
        if Zfil.endswith(".Z") or Zfil.endswith(".gz"):
            uncomp_fil = files_rw.unzip_gz_Z(Zfil)
        else:
            uncomp_fil = Zfil
        clklis.append(uncomp_fil)

    # converting orbits
    gins_orb_lis = []
    for sp3 in sp3lis:
        gins_orb_lis.append(sp3_2_gins(sp3))

    # converting clocks
    gins_clk_lis = []
    for clk in clklis:
        gins_clk_lis.append(clk_2_gins(clk))

    # cating orbits
    utils.cat(catoutorb_path, *gins_orb_lis)
    sort_orbit_gins(catoutorb_path, sortoutorb_path)

    # test, checking if the cat was OK
    catoutclk_nline = utils.line_count(sortoutorb_path)
    if catoutclk_nline < 7500:
        print("WARN : the 3days orb file (~8000lines) looks small ! ")
        print(catoutclk_nline, "lines found in", sortoutorb_path)

    # cleaning the orbits file
    _ = orbit_cleaner(sortoutorb_path, cleanoutorb_path)
    Nsort = utils.line_count(sortoutorb_path)
    Nclean = utils.line_count(cleanoutorb_path)

    badsatlis = bad_sat_finder(sortoutorb_path, 0)

    if len(badsatlis) != 0:
        print("INFO : ", len(badsatlis), "bad sats.", Nsort - Nclean, "epochs del.")
        print("       bad sats :", badsatlis)

    # cating clocks
    utils.cat(catoutclk_path, *gins_clk_lis)

    if rm_temp_files:
        try:
            os.remove(catoutorb_path)
            os.remove(sortoutorb_path)
            [os.remove(e) for e in sp3lis]
            [os.remove(e) for e in clklis]
            [os.remove(e) for e in gins_orb_lis]
            [os.remove(e) for e in gins_clk_lis]
        except:
            pass
    return cleanoutorb_path, catoutclk_path


def export_results_gins_listing(
    listings_list_in, outpath, static_or_kinematic="kine", outprefix="", coordtype="FLH"
):
    if len(listings_list_in) == 0:
        print("ERR : export_results_gins_listing : listings list is empty ...")
        return None
    if static_or_kinematic == "stat":
        ts = files_rw.read_gins_multi_raw_listings(listings_list_in)
    elif static_or_kinematic == "kine":
        tslist = []
        for lising in listings_list_in:
            tslist.append(files_rw.read_gins(lising))
        ts = time_series.merge_ts(tslist)
    time_series.export_ts(ts, outpath, coordtype, outprefix)
    return outpath


def merge_yaml(yaml1, yaml2, yaml_out=None):

    dic1 = yaml.load(open(yaml1))
    dic2 = yaml.load(open(yaml2))

    def merge(dic1, dic2):
        if isinstance(dic1, dict) and isinstance(dic2, dict):
            for k, v in dic2.items():
                if k not in dic1:
                    dic1[k] = v
                else:
                    dic1[k] = merge(dic1[k], v)
        return dic1

    dic3 = merge(dic1, dic2)

    if yaml_out:
        with open(yaml_out, "w+") as outfile:
            outfile.write(yaml.dump(dic3, default_flow_style=False))

    return dic3


#    return dic1

# yaml1    = '/home/psakicki/GINS/gin/data/EXE_PPP/DIR_REF_PPP.yml'
# yaml2    = '/home/psakicki/GINS/gin/TP/TP_RELAX/DIR_MC0.yml'
# yaml_out = '/home/psakicki/GINS/gin/TP/TP_RELAX/DIR_MC0_perso2.yml'
#
# merge_yaml(yaml1,yaml2,yaml_out)

#################### DOUBLE DIFF ####################


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


# def check_if_director_path_is_ok(director_path,gins_path=get_gins_path()):
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
