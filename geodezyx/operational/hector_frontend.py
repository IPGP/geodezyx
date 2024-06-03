#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.operational contains functions to run the 
time series velocities estimation software HECTOR. 

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
import datetime as dt
import glob
#### Import the logger
import logging
import os
import shutil
import subprocess
import time

import matplotlib.pyplot as plt
import numpy as np

from geodezyx import reffram
#### geodeZYX modules
from geodezyx import utils

log = logging.getLogger(__name__)

##########  END IMPORT  ##########
#  _    _ ______ _____ _______ ____  _____
# | |  | |  ____/ ____|__   __/ __ \|  __ \
# | |__| | |__ | |       | | | |  | | |__) |
# |  __  |  __|| |       | | | |  | |  _  /
# | |  | | |___| |____   | | | |__| | | \ \
# |_|  |_|______\_____|  |_|  \____/|_|  \_\
#


def keeping_specific_stats(listoffiles, specific_stats, invert=False):
    """if invert = True : NOT keeping BUT removing specific stats"""
    if not (type(specific_stats) is tuple):
        raise Exception("keeping_specific_stats : specific_stats must be a tuple")
    if not invert:
        outlistoffiles = []
    else:
        outlistoffiles = list(listoffiles)

    if specific_stats != ():
        for fil in listoffiles:
            for stat in specific_stats:
                if stat in fil:
                    if not invert:
                        outlistoffiles.append(fil)
                    else:
                        outlistoffiles.remove(fil)
        return outlistoffiles
    else:
        return listoffiles


def neufile_outlier_removing(
    inp_neufile, generik_conf_file, outdir="", remove_ctl_file=True
):
    """from a NEU file
    => removeoutlier preprocessing
    => 3 MOM files (one per component)"""
    neufile_name = os.path.basename(inp_neufile)
    prefix_inp = neufile_name.replace(".neu", "_")
    inpdir = os.path.dirname(inp_neufile)
    if outdir == "":
        outdir = inpdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    comp_list = ["North", "East", "Up"]

    if utils.grep_boolean(inp_neufile, "Components : YXZ"):
        comp_list_filename = ["Y", "X", "Z"]
    else:
        comp_list_filename = ["N", "E", "U"]

    for comp, comp_fn in zip(comp_list, comp_list_filename):
        work_conf_file = inpdir + "/" + prefix_inp + comp_fn + "_rmoutlier.ctl"
        prefix_out = prefix_inp + comp_fn + "_pre"
        momfile_name = prefix_out + ".mom"
        momfile_path = outdir + "/" + momfile_name
        shutil.copyfile(generik_conf_file, work_conf_file)
        #        work_conf_file_obj = open(work_conf_file,'w+')
        utils.replace(work_conf_file, "DataFile", "DataFile            " + neufile_name)
        utils.replace(work_conf_file, "DataDirectory", "DataDirectory       " + inpdir)
        utils.replace(
            work_conf_file, "OutputFile", "OutputFile          " + momfile_path
        )
        utils.replace(work_conf_file, "component", "component           " + comp)

        p = subprocess.Popen(
            "",
            executable="/bin/bash",
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        command = "removeoutliers " + work_conf_file
        log.info("LAUNCHING : %s", command)
        stdout, stderr = p.communicate(command.encode())
        # logs files
        std_file = open(outdir + "/" + prefix_out + ".std.log", "w")
        std_file.write(stdout.decode("utf-8"))
        std_file.close()
        if stderr:
            log.warning("err.log is not empty, must be checked !")
            log.warning(stderr.decode("utf-8"))
            log.warning("")
            err_file = open(outdir + "/" + prefix_out + ".err.log", "w")
            err_file.write(stderr.decode("utf-8"))
            err_file.close()

        # remove the ctl file
        if remove_ctl_file:
            os.remove(work_conf_file)

    return None


def multi_neufile_outlier_removing(
    inpdir,
    generik_conf_file,
    outdir="",
    extention="neu",
    specific_stats=(),
    invert_specific=False,
    remove_ctl_file=True,
):
    if not os.path.exists(inpdir):
        os.makedirs(inpdir)
    os.chdir(inpdir)

    wildcarded_path = inpdir + "/" + "*." + extention

    listofenufile = glob.glob(wildcarded_path)

    if len(listofenufile) == 0:
        log.warning("no files found in specified directory ...")
        log.warning(wildcarded_path)

    listofenufile = keeping_specific_stats(
        listofenufile, specific_stats, invert=invert_specific
    )

    log.warning("removing outliers of %s timeseries", len(listofenufile))

    for enu in sorted(listofenufile):
        neufile_outlier_removing(
            enu, generik_conf_file, outdir=outdir, remove_ctl_file=remove_ctl_file
        )

    return None


def momfile_trend_processing(
    inp_momfile, generik_conf_file, outdir="", remove_ctl_file=True
):
    momfile_inp_name = os.path.basename(inp_momfile)
    inpdir = os.path.dirname(inp_momfile)
    if outdir == "":
        outdir = inpdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    work_conf_file = inpdir + "/" + momfile_inp_name.replace(".mom", "") + "_trend.ctl"
    prefix_out_name = momfile_inp_name.replace("_pre.mom", "")
    momfile_out_name = momfile_inp_name.replace("pre.mom", "out.mom")
    momfile_out_path = outdir + "/" + momfile_out_name
    shutil.copyfile(generik_conf_file, work_conf_file)
    utils.replace(work_conf_file, "DataFile", "DataFile            " + momfile_inp_name)
    utils.replace(work_conf_file, "DataDirectory", "DataDirectory       " + inpdir)
    utils.replace(
        work_conf_file, "OutputFile", "OutputFile          " + momfile_out_path
    )

    p = subprocess.Popen(
        "",
        executable="/bin/bash",
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    command = "estimatetrend " + work_conf_file
    log.info("LAUNCHING : %s", command)
    stdout, stderr = p.communicate(command.encode())
    # logs files
    sumfile = outdir + "/" + prefix_out_name + ".sum"
    std_file = open(sumfile, "w")
    std_file.write(stdout.decode("utf-8"))
    std_file.close()
    if stderr:
        log.warning("err.log is not empty, must be checked !")
        log.warning(stderr.decode("utf-8"))
        log.warning("")
        err_file = open(outdir + "/" + prefix_out_name + ".err.log", "w")
        err_file.write(stderr.decode("utf-8"))
        err_file.close()

    # plot a graph
    # matplotlib.use('Agg')
    fig = plt.figure()
    try:
        ax = fig.gca()
        ax.set_xlabel("Date")
        ax.set_ylabel("Displacement (m)")
        data = np.loadtxt(momfile_out_path)
        ax.plot(MJD2dt(data[:, 0]), data[:, 1], "b+")
        ax.plot(MJD2dt(data[:, 0]), data[:, 2], "g.")

        # add offsets
        mom_obj = open(inp_momfile, "r+")
        offset_stk = []
        for line in mom_obj:
            if "offset" in line:
                offset_stk.append(float(line.split()[-1]))
        log.info("offset detected for plot : %s", offset_stk)
        for off in offset_stk:
            plt.axvline(MJD2dt(off), color="r")
        mom_obj.close()

        # add trend & bias
        mom_obj = open(sumfile, "r+")
        for line in mom_obj:
            try:
                if "bias :" in line[0:7]:
                    biasline = line
                if "trend:" in line[0:7]:
                    trendline = line
            except:
                continue
        mom_obj.close()

        #    .text(0.1, 0.15, biasline , horizontalalignment='left',verticalalignment='center',transform=ax.transAxes)
        plt.figtext(
            0.1,
            0.05,
            trendline + biasline,
            horizontalalignment="left",
            verticalalignment="center",
            transform=ax.transAxes,
        )
        plt.suptitle(prefix_out_name)

        # save
        fig.set_size_inches(16.53, 11.69 * 0.6667)
        fig.autofmt_xdate()
        fig.savefig(outdir + "/" + prefix_out_name + ".plt.pdf")
        fig.savefig(outdir + "/" + prefix_out_name + ".plt.png")

    except:
        plt.close(fig)

    # remove the ctl file
    if remove_ctl_file:
        os.remove(work_conf_file)

    return None


def MJD2dt(mjd_in):
    # cf http://en.wikipedia.org/wiki/Julian_day
    try:
        return [dt.datetime(1858, 11, 17) + dt.timedelta(m) for m in mjd_in]
    except:
        return dt.datetime(1858, 11, 17) + dt.timedelta(mjd_in)


def multi_momfile_trend_processing(
    inpdir,
    generik_conf_file,
    outdir="",
    extention="pre.mom",
    remove_ctl_file=True,
    specific_stats=(),
    invert_specific=False,
):
    start = time.time()
    if not os.path.exists(inpdir):
        os.makedirs(inpdir)
    os.chdir(inpdir)
    listofmomfile = glob.glob(inpdir + "/" + "*" + extention)
    listofmomfile = keeping_specific_stats(
        listofmomfile, specific_stats, invert=invert_specific
    )

    log.info("processing % files", len(listofmomfile))

    for imom, mom in enumerate(sorted(listofmomfile)):
        log.info("processing file %s/%s", imom + 1, len(listofmomfile))
        momfile_trend_processing(
            mom, generik_conf_file, outdir=outdir, remove_ctl_file=remove_ctl_file
        )

    log.info("multi_momfile_trend_processing exec time : %s", time.time() - start)
    return None


# ===== FCTS 4 multi_sumfiles_trend_extract =====


def sumfiles_to_statdico(inpdir, specific_stats=(), invert_specific=False):
    """this fct search for every sum file in a folder
    a stat dico contains no data
    only the paths to the E,N,U sum files
    statdico[stat] = [path/E.sum,path/N.sum,path/U.sum ]
     for each stat getting the 3 ENU sum files

    Thoses lists will be send in sumfiles_trend_extract
    """
    if not os.path.exists(inpdir):
        os.makedirs(inpdir)
    os.chdir(inpdir)
    listofsumfile = glob.glob(inpdir + "/" + "*" + ".sum")
    listofsumfile = keeping_specific_stats(
        listofsumfile, specific_stats, invert=invert_specific
    )
    listofstat = [
        os.path.basename(sumfil).split(".")[0].split("_")[-2]
        for sumfil in listofsumfile
    ]
    statdico = {}
    for stat in listofstat:
        statdico[stat] = []
        for sumfil in listofsumfile:
            if stat in sumfil:
                statdico[stat].append(sumfil)
    return statdico


def sumfiles_trend_extract(listof3ENUsumfile):
    V, sV = {}, {}
    if len(listof3ENUsumfile) != 3:
        raise Exception("listofENUsumfile != 3")
    for fil in listof3ENUsumfile:
        bnfil = os.path.basename(fil)
        coord = bnfil.split(".")[0][-1]
        for l in open(fil):
            if "trend:" in l:
                f = l.split()
                V[coord] = float(f[1])
                sV[coord] = float(f[3])

    statname = bnfil.split(".")[0].split("_")[-2]
    return statname, V, sV


def get_FLH_from_NEUfile(neufilepath):
    for l in open(neufilepath):
        f = l.split()
        if "Longitude" in l:
            lon = float(f[-1])
        if "Latitude" in l:
            lat = float(f[-1])
        if "Height" in l:
            hau = float(f[-1])
    return lat, lon, hau


def velfile_from_a_list_of_statVsV_tuple(
    listoftup, out_dir, out_prefix, raw_neu_dir="", style="epc"
):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    """
    make a GLOBK style .vel file
    """

    ### HEADER DEFINITION
    if style == "globk":
        outfile = open(out_dir + "/" + out_prefix + ".vel", "w+")
        outfile.write(
            "* Velocity field created with Hector Runner / P. Sakic - La Rochelle Univ. (FRA) \n"
        )
        outfile.write(
            "*  Long         Lat        Evel    Nvel    dEv     dNv    E +-    N +-    Rne      Hvel     dHv    H +-  Site\n"
        )
        outfile.write(
            "*  deg          deg       mm/yr   mm/yr   mm/yr   mm/yr   mm/yr   mm/yr            mm/yr   mm/yr   mm/yr\n"
        )
        k = 1000.0
    elif style == "epc":
        outfile = open(out_dir + "/" + out_prefix + ".vel.neu", "w+")
        k = 1
    elif style == "csv_renag":
        outfile = open(out_dir + "/" + out_prefix + ".vel.csv", "w+")
        outfile.write("Station,V_North,V_East,V_Up,sV_North,sV_East,sV_Up\n")
        k = 1000.0
    elif style == "csv_renag_xyz":
        outfile = open(out_dir + "/" + out_prefix + ".vel.csv", "w+")
        outfile.write("Station,V_X,V_Y,V_Z,sV_X,sV_Y,sV_Z\n")
        k = 1000.0
    elif style == "dataframe":
        column_names = [
            "Station",
            "Latitude",
            "Longitude",
            "V_North",
            "V_East",
            "V_Up",
            "sV_North",
            "sV_East",
            "sV_Up",
        ]
        lines_stk = []
        k = 1

    ### FILLING DEFINITION
    for tup in listoftup:
        stat, V, sV = tup
        try:
            neufile = utils.regex2filelist(raw_neu_dir, stat)[0]
            lat, lon, hau = get_FLH_from_NEUfile(neufile)
        except:
            lat, lon, hau = 0, 0, 0
        if style == "globk":
            line = "{:11.5f}{:11.5f} {:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f} {:6.3f}  {:8.2f}{:8.2f}{:8.2f} {}\n".format(
                lon,
                lat,
                V["E"] * k,
                V["N"] * k,
                0.0,
                0.0,
                sV["E"] * k,
                sV["N"] * k,
                0.0,
                V["U"] * k,
                0.0,
                sV["U"] * k,
                stat + "_GPS",
            )
        elif style == "epc":
            line = "{} {} {} {} {} {} {} {}\n".format(
                stat,
                lat,
                reffram.wrapTo180(lon),
                V["N"] * k,
                V["E"] * k,
                sV["N"] * k,
                sV["E"] * k,
                0,
            )
        elif style == "csv_renag":
            line = "{},{:10.5f},{:10.5f},{:10.5f},{:10.5f},{:10.5f},{:10.5f}\n".format(
                stat,
                V["N"] * k,
                V["E"] * k,
                V["U"] * k,
                sV["N"] * k,
                sV["E"] * k,
                sV["U"] * k,
            )
        elif style == "csv_renag_xyz":
            line = "{},{:10.5f},{:10.5f},{:10.5f},{:10.5f},{:10.5f},{:10.5f}\n".format(
                stat,
                V["X"] * k,
                V["Y"] * k,
                V["Z"] * k,
                sV["X"] * k,
                sV["Y"] * k,
                sV["Z"] * k,
            )
        elif style == "dataframe":
            line = [
                stat,
                lat,
                lon,
                V["N"] * k,
                V["E"] * k,
                V["U"] * k,
                sV["N"] * k,
                sV["E"] * k,
                sV["U"] * k,
            ]
            lines_stk.append(line)

        if style != "dataframe":
            outfile.write(line)

    #### FINALISATION
    if style == "dataframe":
        DF = pd.DataFrame(lines_stk, columns=column_names)
        utils.pickle_saver(DF, out_dir, out_prefix + "_DataFrame")
        return DF
    else:
        outfile.close()
        return None


# def velNEUfile_from_a_list_of_statVsV_tuple(listoftup,out_dir, out_prefix , raw_neu_dir=''):
#    """ make a dirty velocity file compatible with EPC """
#    outfile = open(out_dir +'/' + out_prefix + '.vel.neu','w+')
#    for tup in listoftup:
#        stat , V , sV = tup
#        try:
#            neufile = utils.regex2filelist(raw_neu_dir,stat)[0]
#            lat , lon , hau = get_FLH_from_NEUfile(neufile)
#        except:
#            lat , lon , hau = 0,0,0
#
#        line = '{} {} {} {} {} {} {} {}\n'.format(stat,lat,refframe.wrapTo180(lon),V['N'],V['E'],sV['N'],sV['E'],0)
#        outfile.write(line)
#    return None


def multi_sumfiles_trend_extract(
    inp_dir,
    out_dir,
    out_prefix,
    raw_neu_dir="",
    specific_stats=(),
    invert_specific=False,
    style="epc",
):
    """style = globk OR epc
    make a GLOBK style .vel file or
    make a dirty velocity file compatible with EPC"""

    utils.create_dir(out_dir)

    statdico = sumfiles_to_statdico(
        inp_dir, specific_stats, invert_specific=invert_specific
    )
    listoftup = []
    for stat, listof3ENUsumfile in statdico.items():
        tup = sumfiles_trend_extract(listof3ENUsumfile)
        if tup[0] != stat:
            raise Exception("tup[0] != stat")
        listoftup.append(tup)
    listoftup.sort(key=lambda x: x[0])
    velfile_from_a_list_of_statVsV_tuple(
        listoftup, out_dir, out_prefix, raw_neu_dir, style
    )
