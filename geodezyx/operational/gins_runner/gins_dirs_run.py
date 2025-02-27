#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26/02/2025 11:07:24

@author: psakic
"""


#### Import the logger
import shutil
import os
import sys
import time
import datetime as dt
import subprocess
import glob
import collections
import multiprocessing as mp

from bs4 import UnicodeDammit

from geodezyx import files_rw, time_series

#### Import geodezyx GINS submodules
import geodezyx.operational.gins_runner.gins_common as gynscmn


#### geodeZYX modules
from geodezyx import utils

#### Import the logger
import logging

log = logging.getLogger("geodezyx")


def run_directors(
    dir_paths_inp,
    opts_gins_pc="",
    opts_gins_90="  ",
    version="OPERA",
    cmd_mode="exe_gins_dir",
    force=False,
):
    """
    NEW FCT WHICH CAN MANAGE BOTH ONE RINEX OR A LIST OF RINEX, OR A FIC file (170613)

    mode = "ginspc" or "exe_gins_dir" or "exe_gins_fic"
    """
    # Multi or Single Mode ?
    if type(dir_paths_inp) is list:
        multimode = True
        director_path_lis = dir_paths_inp
        print("******** DIRECTORS RUNS *******")
    elif type(dir_paths_inp) is str:
        multimode = False
        director_path_lis = [dir_paths_inp]
    else:
        print("ERR : run_directors : check the rinex_paths_in !!!")
        return None

    ndir = len(director_path_lis)

    for i, director_path in enumerate(director_path_lis):
        if multimode:
            print(" ======== ", i + 1, "/", ndir, " ======== ")
        print("INFO : launching : ", director_path)
        print("INFO : start at", dt.datetime.now())
        if director_path[-4:] == ".fic":
            cmd_mode = "exe_gins_fic"
            print("INFO : input file ends with .fic, fic_mode is activated")
        start = time.time()

        if director_path[-1] == "~":
            print("INFO : geany ~ temp file, skiping this dir. ")
            continue

        sols_exist = gynscmn.check_solution(os.path.basename(director_path))
        if len(sols_exist) > 0 and not force:
            print("INFO : solution", sols_exist, "already exists, skipping ...")
            continue

        director_name = os.path.basename(director_path)
        opts_gins_pc_ope = "-F" + opts_gins_pc

        print("INFO : options ginsPC / gins90 :", opts_gins_pc_ope, "/", opts_gins_90)

        if "IPPP" in opts_gins_90 and cmd_mode != "ginspc":
            _check_dir_keys(director_path)

        if cmd_mode == "ginspc":
            command = " ".join(
                (
                    "ginspc.bash",
                    opts_gins_pc_ope,
                    director_name,
                    opts_gins_90,
                    "-v",
                    version,
                )
            )
        elif cmd_mode == "exe_gins_fic":
            command = " ".join(
                ("exe_gins", "-fic", director_name, "-v", version, opts_gins_90)
            )
        elif cmd_mode == "exe_gins_dir":
            command = " ".join(
                ("exe_gins", "-dir", director_name, "-v", version, opts_gins_90)
            )
        else:
            print("ERR : run_directors : mode not recognized !!!")
            return None

        # l'argument OPERA est super important  !!!
        # c'est lui qui resoult l'instabilité lors du lancement du subprocess !!!
        # parce que indirectement exe_gins ne marche pas sans
        # faire le test avec un exe_gins -v OPERA -fic <fic> et sans

        print("INFO : submit. command : ", command)

        gins_path = gynscmn.get_gin_path()
        log_path = utils.create_dir(os.path.join(gins_path, "python_logs"))
        log_path = os.path.join(gins_path, "python_logs", director_name + ".log")

        # with open(log_path, "w+") as f:
        #
        #     for c in iter(lambda: process.stdout.read(1), ""):
        #         # https://stackoverflow.com/questions/436220/determine-the-encoding-of-text-in-python4
        #         # find encoding
        #         # print(c)
        #         # dammit = UnicodeDammit(c)
        #         # print(dammit.original_encoding)
        #         d_deco = c.decode("iso-8859-1", errors="replace")
        #         # sys.stdout.write(d_deco)
        #         f.write(d_deco)

        process = subprocess.Popen(
            [command],
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            executable="/bin/bash",
        )

        gynscmn.check_gins_exe(open(log_path), director_name)

        with open(log_path + ".exec", "w") as f:
            f.write("exec time : " + str(time.time() - start))

        print("INFO : end at", dt.datetime.now())
        print("INFO : exec time : ", str(time.time() - start), "seconds")

    return None


def run_director_list_wrap(tupinp):
    run_directors(*tupinp)
    return None


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
):

    kwargs_lis = []
    for dirr in list(sorted(dir_paths_inp)):
        kwargs = dict()
        kwargs["dir_paths_inp"] = dirr
        kwargs["opts_gins_pc"] = opts_gins_pc
        kwargs["opts_gins_90"] = opts_gins_90
        kwargs["version"] = version
        kwargs["cmd_mode"] = cmd_mode
        kwargs_lis.append(kwargs)

    pool = mp.Pool(processes=nprocs)
    # res_raw = [pool.apply(run_dirs_kwwrap, args=(x,)) for x in kwargs_lis]
    res_raw = pool.map(run_dirs_kwwrap, kwargs_lis)

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


def _check_dir_keys(director_path_inp):
    for grepstr in (
        "userext_gps__qualiteorb",
        "userext_gps__haute_freq",
        "userext_gps__hor_interp",
        "GPS__QUALITEORB",
        "GPS__HAUTE_FREQ",
        "GPS__HOR_INTERP",
    ):
        grep_out = utils.grep(director_path_inp, grepstr)
        if grep_out == "":
            print("WARN : IPPP mode on, but no", grepstr, " in the dir !!!")
