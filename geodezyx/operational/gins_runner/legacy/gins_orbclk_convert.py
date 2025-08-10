#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26/02/2025 11:36:09

@author: psakic
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
#### Import the logger
import logging
import os
import subprocess

#### geodeZYX modules
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import operational
from geodezyx import utils

import geodezyx.operational.gins_runner.gins_common as ginscmn

log = logging.getLogger('geodezyx')

def sp3_2_gins(sp3_pathin):
    os.chdir(os.path.dirname(sp3_pathin))
    # print("WARN : the path of the GINS conversion tool has been hardcoded !!!")
    kommand_1 = ginscmn.get_gin_path() + "/gins_toolbox/scripts/"
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
    kommand_1 = ginscmn.get_gin_path() + "/gins_toolbox/scripts/"
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
        temp_dir = os.path.join(ginscmn.get_gin_path(), "gin", "TEMP_DATA")
        log.info("use default forder %s", temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    strtdt = centraldate - dt.timedelta(days=1)
    enddt = centraldate + dt.timedelta(days=1)

    if centraldate > dt.datetime(2013, 12, 29) and repro != 0:
        log.warning(str(centraldate) + " a bit late for repro " + str(repro))

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
    sp3_zlis = operational.multi_downloader_orbs_clks_2(
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
    clk_zlis = operational.multi_downloader_orbs_clks_2(
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
    if not clk_zlis:
        print("INFO : clock download : no 30s clock found, trying std. clocks")
        clk_zlis = operational.multi_downloader_orbs_clks(
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
    if len(sp3_zlis) < 3 or len(clk_zlis) < 3:
        print(
            "ERR : download_convert_2_gins_orb_clk : missing sp3/clk for", centraldate
        )
        print("list of sp3 : ", sp3_zlis)
        print("list of clk : ", clk_zlis)
        return None, None

    # uncompressing orbits
    sp3lis = []
    for Zfil in sp3_zlis:
        if Zfil.endswith(".Z") or Zfil.endswith(".gz"):
            uncomp_fil = files_rw.unzip_gz_z(Zfil)
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
    for Zfil in clk_zlis:
        if Zfil.endswith(".Z") or Zfil.endswith(".gz"):
            uncomp_fil = files_rw.unzip_gz_z(Zfil)
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