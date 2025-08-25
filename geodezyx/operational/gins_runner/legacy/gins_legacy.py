#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 19/04/2025 17:32:02

@author: psakic
"""

import os
import shutil
import subprocess
import time
import yaml
import glob

import geodezyx.operational.gins_runner.gins_common as gynscmn




#   ____  _     _   _                                    __     _
#  / __ \| |   | | | |                                  / _|   | |
# | |  | | | __| | | |     ___  __ _  __ _  ___ _   _  | |_ ___| |_ ___
# | |  | | |/ _` | | |    / _ \/ _` |/ _` |/ __| | | | |  _/ __| __/ __|
# | |__| | | (_| | | |___|  __/ (_| | (_| | (__| |_| | | || (__| |_\__ \
#  \____/|_|\__,_| |______\___|\__, |\__,_|\___|\__, | |_| \___|\__|___/
#                               __/ |            __/ |
#                              |___/            |___/


############### OLD FCTS
def write_oclo_file(station_file, oceanload_out_file, fes_yyyy=2004):
    temp_cmd_fil = os.path.join(os.path.dirname(oceanload_out_file), "oclo.cmd.tmp")
    temp_cmd_filobj = open(temp_cmd_fil, "w")

    temp_cmd_filobj.write(station_file + "\n")
    temp_cmd_filobj.write(oceanload_out_file + "\n")
    if fes_yyyy == 2012:
        exe_loadoce_cmd = "exe_loadoce_fes2012"
    elif fes_yyyy == 2004:
        exe_loadoce_cmd = "exe_loadoce"
    else:
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



def dirs_copy_generik_2_working(
    director_name_prefix, director_generik_path, temp_data_folder=None
):

    if not temp_data_folder:
        temp_data_folder = gynscmn.get_temp_data_gins_path()

    shutil.copy(director_generik_path, temp_data_folder)
    director_generik_path_tmp = os.path.join(
        temp_data_folder, os.path.basename(director_generik_path)
    )
    director_generik_path_out = os.path.join(
        temp_data_folder, director_name_prefix + ".yml"
    )
    os.rename(director_generik_path_tmp, director_generik_path_out)

    return director_generik_path_out


def get_director_list(wildcard_dir):
    """with a wildcard (e.g. 'GWADA_MK2*') and return a list of corresponding
    directors found in gin/data/directeur folder"""
    gins_path = gynscmn.get_gin_path()
    di_run_lis = [
        os.path.basename(e)
        for e in glob.glob(os.path.join(gins_path, "data", "directeur", wildcard_dir))
    ]
    return di_run_lis


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
        stat_up = stat.upper()
        archiv = os.path.join(archive_path, stat_up)
        if not os.path.exists(archiv):
            os.makedirs(archiv)
        shutil.move(f, archiv)

    return None


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


def check_gins_exe(streamin, director_name):
    """
    Check the execution status of a GINS director.

    Parameters
    ----------
    streamin : file-like object
        The input stream to read the execution output from.
    director_name : str
        The name of the director being checked.

    Returns
    -------
    bool
        True if the execution was successful, False otherwise.
    """
    if "Exécution terminée du fichier" in streamin.read():
        print("INFO : happy end for " + director_name + " :)")
        return True
    else:
        print("WARN : bad end for " + director_name + " :(")
        return False

def run_dirs_kwwrap(kwarg):
    run_directors(**kwarg)
    return None


def run_dirs_multi(
    dir_paths_inp,
    nprocs=4,
    opts_gins_pc="",
    opts_gins_90="",
    version="OPERA",
    cmd_mode="exe_gins_dir",
    force=False,
):
    kwargs_lis = []
    for dirr in list(sorted(dir_paths_inp)):
        kwargs = dict()
        kwargs["dir_paths_inp"] = dirr
        kwargs["opts_gins_pc"] = opts_gins_pc
        kwargs["opts_gins_90"] = opts_gins_90
        kwargs["version"] = version
        kwargs["cmd_mode"] = cmd_mode
        kwargs["force"] = force
        kwargs_lis.append(kwargs)

    pool = mp.Pool(processes=nprocs)
    # res_raw = [pool.apply(run_dirs_kwwrap, args=(x,)) for x in kwargs_lis]
    res_raw = pool.map(run_dirs_kwwrap, kwargs_lis, chunksize=1)

    return None



def run_dirs_multislots_custom(
    director_lis,
    slots_lis=["", "U", "L", "R"],
    opts_gins_pc="",
    opts_gins_90="",
    version="OPERA",
    mode="ginspc",
):
    if not type(director_lis) is list:
        log.error("director_lis in input is not a list !!!")
        return None
    log.info("TOTAL NB OF DIRECTORS : %i", len(director_lis))

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
            [mode] * len(slots_lis),
        )
    )
    pool.map(run_director_list_wrap, ZIP)
    return None


def smart_directors_to_run(wildcard_dir="", full_path_out=True):
    """smart runner check directors who worked, and thus give a list of directors
    whitout those who worked

    listing and directeur folders are inspected automatically"""

    gins_path = gynscmn.get_gin_path(True)

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
        dirpath = os.path.join(gynscmn.get_gin_path(True), "data", "directeur")
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
        gynscmn.get_gin_path(), "gin", "batch", "listing", wildcard_dir
    )
    listing_lis = [e for e in glob.glob(listing_path)]

    prepars_lis = [e for e in listing_lis if "prepars" in e]
    gins_lis = [e for e in listing_lis if "gins" in e]

    dirpath = os.path.join(
        gynscmn.get_gin_path(True), "data", "directeur", wildcard_dir
    )
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


def export_results_gins_listing(
    listings_list_in, outpath, static_or_kinematic="kine", outprefix="", coordtype="FLH"
):
    from geodezyx import time_series
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


def run_director_list_wrap(tupinp):
    run_directors(*tupinp)
    return None




#  ______                _   _                                                             _
# |  ____|              | | (_)                                                           | |
# | |__ _   _ _ __   ___| |_ _  ___  _ __     __ _ _ __ __ ___   _____ _   _  __ _ _ __ __| |
# |  __| | | | '_ \ / __| __| |/ _ \| '_ \   / _` | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
# | |  | |_| | | | | (__| |_| | (_) | | | | | (_| | | | (_| |\ v /  __/ |_| | (_| | | | (_| |
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \__, |_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
#                                             __/ |                     __/ |
#                                            |___/                     |___/


#    return dic1

# yaml1    = '/home/psakicki/GINS/gin/data/EXE_PPP/DIR_REF_PPP.yml'
# yaml2    = '/home/psakicki/GINS/gin/TP/TP_RELAX/DIR_MC0.yml'
# yaml_out = '/home/psakicki/GINS/gin/TP/TP_RELAX/DIR_MC0_perso2.yml'
#
# merge_yaml(yaml1,yaml2,yaml_out)


###### FUNCTIONS GRAVEYARD ########

# def get_rinex_list(
#     parent_folder,
#     specific_stats=[],
#     invert=False,
#     compressed=True,
#     start=dt.datetime(1980, 1, 1),
#     end=dt.datetime(2099, 1, 1),
# ):
#     """return :
#         all RINEXs found in a parent folder and his subfolders
#         (compressed or not)
#
#     parent_folder :
#         can be a the path of the partent folder (the rinex archive path)
#         but also a RINEX list already found
#         (modification of 170524 to gain speed)
#
#     specific_stats :
#         MUST BE a list ['STA1','STA2','STA3']
#         So, if only one elt in specific_stats, use a syntax as : ['STA1']
#
#     invert :
#         False = keeping the specific stats
#         True  = removing the specific stats
#         or for all stations, leave a empty
#         tuple in specific_stats
#
#     NB : the end date is included
#     """
#     if type(parent_folder) is str:
#         wholefilelist = []
#         for root, dirs, files in os.walk(parent_folder, topdown=False):
#             for name in files:
#                 wholefilelist.append(os.path.join(root, name))
#
#         wholefilelist = list(set(wholefilelist))
#
#         if compressed:
#             rnxregex = operational.rinex_regex()
#         else:
#             rnxregex = operational.rinex_regex(False)
#
#         rinexfilelist = [fil for fil in wholefilelist if re.search(rnxregex, fil)]
#     elif utils.is_iterable(parent_folder):
#         # is parent_folder is a list containing already found RINEXs
#         rinexfilelist = parent_folder
#
#     specific_stats = [e.lower() for e in specific_stats]
#
#     if specific_stats != []:
#         if not invert:
#             goodrnxfilelist = [
#                 fil
#                 for fil in rinexfilelist
#                 if os.path.basename(fil)[0:4] in specific_stats
#             ]
#         else:
#             goodrnxfilelist = [
#                 fil
#                 for fil in rinexfilelist
#                 if os.path.basename(fil)[0:4] not in specific_stats
#             ]
#     else:
#         goodrnxfilelist = rinexfilelist
#
#     goodrnxfilelist = sorted(goodrnxfilelist)
#
#     if goodrnxfilelist == []:
#         print(
#             "WARN : get_rinex_list : no RINEX found, check the path of the parent folder !"
#         )
#     else:
#         print("INFO : get_rinex_list :", len(goodrnxfilelist), "RINEXs found")
#
#     goodrnxfilelist_date = []
#     end = end + dt.timedelta(seconds=86399)
#     for rnx in goodrnxfilelist:
#         dtrnx = conv.rinexname2dt(os.path.basename(rnx))
#         if start <= dtrnx <= end:
#             goodrnxfilelist_date.append(rnx)
#     if len(goodrnxfilelist_date) != len(goodrnxfilelist):
#         dif = len(goodrnxfilelist) - len(goodrnxfilelist_date)
#         print(
#             "INFO : get_rinex_list :",
#             dif,
#             "RINEXs removed bc. not in the date interval",
#         )
#
#     if goodrnxfilelist_date == []:
#         print("WARN : get_rinex_list : all RINEXs removed bc. of date interval")
#
#     return goodrnxfilelist_date
