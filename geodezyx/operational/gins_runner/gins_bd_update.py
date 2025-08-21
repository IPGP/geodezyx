#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27/02/2025 10:37:26

@author: psakic
"""

import argparse
import getpass
import os
import subprocess
import datetime as dt
import re
import sys

from geodezyx import conv, utils
import geodezyx.operational.gins_runner.gins_common as gynscmn

#### Import the logger
import logging

log = logging.getLogger("geodezyx")


def read_credentials(file_path):
    with open(file_path, "r") as file:
        login = file.readline().strip()
        password = file.readline().strip()
    return login, password


def download_rsync(
    file_list,
    remote_user,
    remote_host,
    remote_path,
    local_destination,
    rsync_options=None,
    password=None,
):
    """
    Downloads a list of files using rsync.
    """
    if not rsync_options:
        rsync_options = [
            "--archive",
            "--verbose",
            "--compress",
            "--progress",
            "--relative",
            "--copy-links",
            "--no-perms",
            "--no-owner",
            "--no-group",
            "--omit-dir-times",
        ]  # Default options: archive mode, verbose, and compression

    # Construct the remote source path
    remote_source = f"{remote_user}@{remote_host}:{remote_path}"

    if password:
        rsync_base = ["sshpass", "-p", password, "rsync"] + rsync_options
    else:
        rsync_base = ["rsync"] + rsync_options

    # if exclude_compress:
    #   '--ignore-existing' '--exclude-from=<(echo "$(find /destination -name '*.gz' | sed 's/\.gz$//')"')'

    tmp_rsync_file_lis = f"/tmp/" + utils.get_timestamp() + "_tmp_rsync_file_list.lst"
    utils.write_in_file("\n".join(file_list), tmp_rsync_file_lis)

    # Construct the rsync command
    # /./ : https://askubuntu.com/questions/552120/preserve-directory-tree-while-copying-with-rsync
    rsync_cmd = (
        rsync_base
        + ["--files-from", tmp_rsync_file_lis]
        + [f"{remote_source}/./", f"{local_destination}/"]
    )

    log.debug("Rsync command: %s", " ".join(rsync_cmd))

    # Run the rsync command
    process = subprocess.run(rsync_cmd, stderr=sys.stderr, stdout=sys.stdout, text=True)

    if process.returncode != 0:
        log.error("Rsync failed with return code :( %d", process.returncode)
    else:
        log.info("Rsync completed successfully :)")

    os.remove(tmp_rsync_file_lis)
    return


def bdgins_update(
    date_srt=dt.datetime(2020, 5, 3),
    date_end=None,
    dir_bdgins="",
    login="",
    password="",
    compress=False,
):
    """
    Update the BDGINS repository with the necessary files
    for the given date range.

    Parameters
    ----------
    date_srt : datetime
        The start date for the update.
        The default is 2020-05-03 (origin of G20 products, end of GPS's Selectable Avaiability).
    date_end : datetime
        The end date for the update.
        If not provided, the default is the 15 days ago.
    dir_bdgins : str, optional
        The directory for BDGINS. Defaults to an empty string.
    login : str, optional
        The login for the remote server. Defaults to an empty string.
    password : str, optional
        The password for the remote server. Defaults to an empty string.
    compress : bool, optional
        Whether to gzip-compress the clock files.
        experimental, and not recommended
        Defaults to False.
    
    Returns
    -------
    None
    """
    if not date_end:
        date_end = dt.datetime.now() - dt.timedelta(days=15)
        date_end = conv.round_dt(date_end, "1D", mode="floor", python_dt_out=False)
        date_end = conv.numpy_dt2dt(date_end)

    if not dir_bdgins:
        dir_bdgins = os.path.join(gynscmn.get_gin_path(True), "data")

    if not login:
        home = os.path.expanduser("~")
        ginspcrc = os.path.join(home, ".ginspc")
        if os.path.isfile(ginspcrc):
            login = open(ginspcrc, "r").read().split("=")[1].strip()
        else:
            login = getpass.getuser()

    date = date_srt

    ###### LIST INITIALISATION
    ### misc files
    list_misc = []

    ### full folders
    list_antex = []
    list_constell = []
    list_exe_ppp = []
    list_lunisolaires = []
    list_macromod = []
    list_maree_polaire = []
    list_pole = []
    list_prairie = []

    ### time dependant files
    list_tropo_vmf1 = []
    list_iono_igs = []
    list_orbite_mg3 = []

    list_orbite_g20 = []
    list_orbex_g20 = []
    list_horl_g20 = []

    ###### LIST FILL
    ### full folders
    # (folder's path is added in the rsync command, with subdir destination variable
    ## list_antex.extend(["/"]) ## temporary off because updated continuously
    ## list_constell.extend(["/"])
    ## list_exe_ppp.extend(["/"])
    ## list_lunisolaires.extend(["/"]) ## temporary off because updated continuously
    ## list_macromod.extend(["/"])
    ## list_maree_polaire.extend(["/"])
    ## list_pole.extend(["/"])
    list_prairie.extend(["/"])

    ### misc files: the needed files are considered individually
    # some of them are a redundancy since they must be downloaded in the full folders
    list_misc.extend([f"prairie/igs_satellite_metadata.snx"])
    list_misc.extend([f"pole/nominal_NRO"])
    list_misc.extend([f"ANTEX/igs20.atx"])
    list_misc.extend([f"EXE_PPP/valap_static"])
    list_misc.extend([f"macromod/gnss.xml"])
    list_misc.extend([f"lunisolaires/de440bdlf.ad"])
    list_misc.extend([f"maree_polaire/loading/nominal"])
    l_fil_cons = [
        "constellation_gps.infos",
        "histocom.infos",
        "histogal.infos",
        "historik_glonass",
        "igs_satellite_metadata.snx",
    ]
    list_misc.extend(["constell/" + f for f in l_fil_cons])

    # here is an hybrid tropo/misc file
    list_tropo_vmf1.append(f"orography_ell")

    ### time dependant files
    while date <= date_end + dt.timedelta(days=2): # 2 days later is need for cat orb/clk
        day = str(date.day).zfill(2)
        month = str(date.month).zfill(2)
        year = date.year
        yy = str(year)[2:]

        wk, wkday = conv.dt2gpstime(date)
        doy = str(conv.dt2doy(date)).zfill(3)

        list_tropo_vmf1.extend(
            [
                f"{year}/VMFG_{year}{month}{day}.H00",
                f"{year}/VMFG_{year}{month}{day}.H06",
                f"{year}/VMFG_{year}{month}{day}.H12",
                f"{year}/VMFG_{year}{month}{day}.H18",
            ]
        )
        list_iono_igs.append(f"{year}/igsg{doy}0.{yy}i.Z")
        list_orbite_mg3.append(f"mg3{wk}{wkday}.sp3.Ci9PAU")
        list_orbite_g20.append(f"G20{wk}{wkday}.gin")
        list_orbex_g20.append(f"G20{wk}{wkday}.obx.gz")
        list_horl_g20.append(f"hogps_g20{wk}{wkday}")
        date += dt.timedelta(days=1)

    ###### DESTINATION FOLDERS
    dest_subdir_dic = {
        ### full folders
        #"ANTEX": list_antex,
        #"constell": list_constell,
        #"EXE_PPP": list_exe_ppp,
        #"lunisolaires": list_lunisolaires,
        #"macromod": list_macromod,
        #"maree_polaire": list_maree_polaire,
        "prairie": list_prairie,
        ### misc files
        ".": list_misc,  ## for misc files, destination is in the input path (.)
        ### time dependant files
        "tropo_vmf1": list_tropo_vmf1,
        "ionosphere/igs": list_iono_igs,
        "orbites/SP3/re3": list_orbite_mg3,
        "mesures/gps/orbites/G20": list_orbite_g20,
        "mesures/gps/orbex/G20": list_orbex_g20,
        "mesures/gps/horloges30/G20": list_horl_g20,
    }

    create_dir(dir_bdgins, subdirs=dest_subdir_dic.keys())
    os.chdir(dir_bdgins)

    for subdir, files_list in dest_subdir_dic.items():
        subdir_print = "misc. files" if subdir == "." else subdir
        log.info("Downloading files for %s", subdir_print)

        local_subdir_fullpath = os.path.join(dir_bdgins, subdir)
        remote_subdir_fullpath = os.path.join("/home/gins/MIROIR_STAF", subdir)

        download_rsync(
            files_list,
            login,
            "tite.get.obs-mip.fr",
            remote_subdir_fullpath,
            local_subdir_fullpath,
            password,
        )

    # compress the clock files (experimental, and not recommended)
    if compress:
        for dirr in ["mesures/gps/horloges30/G20", "mesures/gps/orbites/G20"]:
            d = os.path.join(dir_bdgins, dirr)
            cmd = [
                "find",
                d,
                "-name",
                "!",
                ".gz",
                "-exec",
                "gzip",
                "-v",
                "{}",
                r"\;",
            ]
            process = subprocess.run(
                cmd, stderr=sys.stderr, stdout=sys.stdout, text=True
            )

    return None


def create_dir(parent_dir, subdirs):
    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)

    for subdir in subdirs:
        subdir_fullpath = os.path.join(parent_dir, subdir)
        if not os.path.exists(subdir_fullpath):
            os.makedirs(subdir_fullpath)


def main():
    parser = argparse.ArgumentParser(description="Update BDGINS repository.")
    parser.add_argument(
        "-s",
        "--date_srt",
        type=lambda s: conv.date_pattern_2_dt(s),
        required=False,
        default=dt.datetime(2020, 5, 3),
        help=(
            "Start date for the update in various formats. "
            "Default is 2020-05-03 "
            "(origin of G20 products, end of GPS's Selective Availability)."
        ),
    )
    parser.add_argument(
        "-e",
        "--date_end",
        type=lambda s: conv.date_pattern_2_dt(s),
        required=False,
        default=None,
        help=(
            "End date for the update in various formats. "
            "Default is 15 days before the current date."
        ),
    )
    parser.add_argument(
        "-d",
        "--dir_bdgins",
        help=(
            "Directory for the BDGINS repository. "
            "Defaults to a predefined path if not provided."
        ),
        required=False,
    )
    parser.add_argument(
        "-l",
        "--login",
        help=(
            "Login for the remote server. "
            "Defaults to the current user or a value from the `.ginspc` file."
        ),
        required=False,
    )
    parser.add_argument(
        "-p",
        "--password",
        help=(
            "Password for the remote server. "
            "Defaults to an empty string if not provided."
        ),
        required=False,
        default=None,
    )
    parser.add_argument(
        "-c",
        "--compress",
        help=(
            "Flag to enable gzip compression of clock files. "
            "Defaults to `False` and is experimental."
        ),
        action="store_true",
        default=False,
    )
    args = parser.parse_args()

    bdgins_update(
        date_srt=args.date_srt,
        date_end=args.date_end,
        login=args.login,
        password=args.password,
        dir_bdgins=args.dir_bdgins,
        compress=args.compress,
    )


if __name__ == "__main__":
    main()
