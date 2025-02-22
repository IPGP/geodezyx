#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.operational contains functions to run the 
GNSS processing software TRACK. 

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

#### geodeZYX modules
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import operational
from geodezyx import utils

log = logging.getLogger('geodezyx')

##########  END IMPORT  ##########


def track_runner(
    rnx_rover,
    rnx_base,
    working_dir,
    experience_prefix,
    XYZbase=[],
    XYZrover=[],
    outtype="XYZ",
    mode="short",
    interval=None,
    antmodfile="~/gg/tables/antmod.dat",
    calc_center="igs",
    forced_sp3_path="",
    const="G",
    silent=False,
    rinex_full_path=False,
    run_on_gfz_cluster=False,
    forced_iono_path="",
):

    # paths & files
    working_dir = utils.create_dir(working_dir)
    temp_dir = utils.create_dir(os.path.join(working_dir, "TEMP"))
    out_dir = utils.create_dir(os.path.join(working_dir, "OUTPUT"))

    if operational.check_if_compressed_rinex(rnx_rover):
        rnx_rover = operational.crz2rnx(rnx_rover, temp_dir)
    else:
        shutil.copy(rnx_rover, temp_dir)

    if operational.check_if_compressed_rinex(rnx_base):
        rnx_base = operational.crz2rnx(rnx_base, temp_dir)
    else:
        shutil.copy(rnx_base, temp_dir)

    # RINEX START & END
    rov_srt, rov_end, rov_itv = operational.rinex_start_end(rnx_rover, 1)
    bas_srt, bas_end, bas_itv = operational.rinex_start_end(rnx_base, 1)

    # RINEX NAMES
    rov_name = os.path.basename(rnx_rover)[0:4]
    bas_name = os.path.basename(rnx_base)[0:4]

    rov_name_uper = rov_name.upper()
    bas_name_uper = bas_name.upper()

    srt_str = rov_srt.strftime("%Y_%j")
    exp_full_name = "_".join((experience_prefix, rov_name, bas_name, srt_str))

    out_conf_fil = os.path.join(out_dir, exp_full_name + ".cmd")
    out_result_fil = os.path.join(out_dir, exp_full_name + ".out")

    log.info(out_conf_fil)

    confobj = open(out_conf_fil, "w+")

    # Obs Files
    confobj.write(" obs_file" + "\n")
    ### just the basename, the caracter nb is limited  (20210415)
    if not rinex_full_path:
        confobj.write(
            " ".join((" ", bas_name_uper, os.path.basename(rnx_base), "F")) + "\n"
        )
        confobj.write(
            " ".join((" ", rov_name_uper, os.path.basename(rnx_rover), "K")) + "\n"
        )
    else:
        confobj.write(" ".join((" ", bas_name_uper, rnx_base, "F")) + "\n")
        confobj.write(" ".join((" ", rov_name_uper, rnx_rover, "K")) + "\n")
    confobj.write("\n")

    date = conv.rinexname2dt(os.path.basename(rnx_rover))

    # Nav File
    if forced_sp3_path == "":
        strt_rnd = dt.datetime(*bas_srt.timetuple()[:3])
        end_rnd = dt.datetime(*bas_end.timetuple()[:3])

        orblis = operational.multi_downloader_orbs_clks(
            temp_dir, strt_rnd, end_rnd, archtype="/", calc_center=calc_center
        )

        # sp3Z = orblis[0]
        sp3 = [utils.uncompress(sp3Z) for sp3Z in orblis]
        sp3 = [e if ".sp3" in e[-5:] else e + ".sp3" for e in sp3]
    else:
        if utils.is_iterable(forced_sp3_path):
            sp3 = forced_sp3_path
        else:
            sp3 = [forced_sp3_path]
    for sp3_mono in sp3:
        confobj.write(" ".join((" ", "nav_file", sp3_mono, " sp3")) + "\n")
    confobj.write("\n")

    # Iono file

    if forced_iono_path != "":
        confobj.write(" ionex_file " + forced_iono_path + "\n")

    # Mode
    confobj.write(" mode " + mode + "\n")
    confobj.write("\n")

    # Output
    confobj.write(" pos_root " + exp_full_name + ".pos" + "\n")
    confobj.write(" res_root " + exp_full_name + ".res" + "\n")
    confobj.write(" sum_file " + exp_full_name + ".sum" + "\n")
    confobj.write("\n")

    # Outtype
    confobj.write(" out_type " + outtype + "\n")
    confobj.write("\n")

    # Interval
    if not interval:
        confobj.write(" interval " + str(rov_itv) + "\n")
    else:
        confobj.write(" interval " + str(interval) + "\n")

    confobj.write("\n")

    # Coords
    bool_site_pos = False
    if XYZbase != []:
        if not bool_site_pos:
            confobj.write(" site_pos \n")
            bool_site_pos = True
        XYZbase = [str(e) for e in XYZbase]
        confobj.write(" ".join([" ", bas_name_uper] + XYZbase + ["\n"]))

    if XYZrover != []:
        if not bool_site_pos:
            confobj.write(" site_pos \n")
            bool_site_pos = True
        XYZrover = [str(e) for e in XYZrover]
        confobj.write(" ".join([" ", rov_name_uper] + XYZrover + ["\n"]))

    if bool_site_pos:
        confobj.write("\n")

    # Offsets
    confobj.write(" ante_off \n")

    Antobj_rov, Recobj_rov, Siteobj_rov, Locobj_rov = files_rw.read_rinex_2_dataobjts(
        rnx_rover
    )

    confobj.write(
        " ".join(
            [
                " ",
                rov_name_uper,
                str(Antobj_rov.North_Ecc),
                str(Antobj_rov.East_Ecc),
                str(Antobj_rov.Up_Ecc),
                Antobj_rov.Antenna_Type,
                "\n",
            ]
        )
    )

    Antobj_bas, Recobj_bas, Siteobj_bas, Locobj_bas = files_rw.read_rinex_2_dataobjts(
        rnx_base
    )

    confobj.write(
        " ".join(
            [
                " ",
                bas_name_uper,
                str(Antobj_bas.North_Ecc),
                str(Antobj_bas.East_Ecc),
                str(Antobj_bas.Up_Ecc),
                Antobj_bas.Antenna_Type,
                "\n",
            ]
        )
    )
    confobj.write("\n")

    # Site_stats
    confobj.write(" site_stats \n")
    confobj.write(" " + bas_name_uper + " 0.1 0.1 0.1 0 0 0" + "\n")
    confobj.write(" " + rov_name_uper + " 20 20 20 0.5 0.5 0.5" + "\n")
    confobj.write("\n")

    # constellqtions
    confobj.write(" TR_GNSS " + const + "\n")

    # Misc
    # confobj.write(" USE_GPTGMF"   + '\n')
    confobj.write(" ATM_MODELC GMF 0.5" + "\n")
    confobj.write(" ANTMOD_FILE " + antmodfile + "\n")
    confobj.write(" DCB_FILE " + "~/gg/incremental_updates/tables/dcb.dat.gnss" + "\n")

    confobj.write(" atm_stats" + "\n")
    confobj.write("  all 0.1 0.00030.00023" + "\n")

    confobj.close()
    # END OF FILE WRITING

    dowstring = "".join([str(e) for e in conv.dt2gpstime(date)])
    bigcomand = " ".join(
        ("track -f", out_conf_fil, "-d", conv.dt2doy(date), "-w", dowstring)
    )

    if run_on_gfz_cluster:
        bigcomand = "cjob -c '" + bigcomand + "'"
        executable = "/bin/csh"
    else:
        executable = "/bin/bash"

    log.info("command launched :")
    log.info(bigcomand)

    # START OF PROCESSING
    if not silent:
        os.chdir(temp_dir)
        try:
            subprocess.call(
                [bigcomand], executable=executable, shell=True, timeout=60 * 20
            )
        except subprocess.TimeoutExpired:
            log.warning("command timeout expired, skip")
            pass

        outfiles = []
        outfiles = outfiles + glob.glob(os.path.join(temp_dir, exp_full_name + "*sum*"))
        outfiles = outfiles + glob.glob(os.path.join(temp_dir, exp_full_name + "*pos*"))
        outfiles = outfiles + glob.glob(os.path.join(temp_dir, exp_full_name + "*cmd*"))

        Antobj_rov, Recobj_rov, Siteobj_rov, Locobj_rov = (
            files_rw.read_rinex_2_dataobjts(rnx_rover)
        )

        [shutil.copy(e, out_dir) for e in outfiles]
        [os.remove(e) for e in outfiles]

        log.info("TRACK RUN FINISHED")
        log.info("results available in %s", out_dir)
    else:
        log.info("Silent mode ON: nothing is launched")

    return bigcomand


def run_track(temp_dir, exp_full_name, out_conf_fil, date, rnx_rover):
    dowstring = "".join([str(e) for e in conv.dt2gpstime(date)])
    bigcomand = " ".join(
        ("track -f", out_conf_fil, "-d", conv.dt2doy(date), "-w", dowstring)
    )

    log.info("INFO : command launched :")
    log.info(bigcomand)

    # START OF PROCESSING
    os.chdir(temp_dir)
    subprocess.call([bigcomand], executable="/bin/bash", shell=True)

    outfiles = []
    outfiles = outfiles + glob.glob(os.path.join(temp_dir, exp_full_name + "*sum*"))
    outfiles = outfiles + glob.glob(os.path.join(temp_dir, exp_full_name + "*pos*"))
    outfiles = outfiles + glob.glob(os.path.join(temp_dir, exp_full_name + "*cmd*"))

    Antobj_rov, Recobj_rov, Siteobj_rov, Locobj_rov = files_rw.read_rinex_2_dataobjts(
        rnx_rover
    )

    [shutil.copy(e, out_dir) for e in outfiles]
    [os.remove(e) for e in outfiles]

    log.info("TRACK RUN FINISHED")
    log.info("results available in %s ", out_dir)

    return None
