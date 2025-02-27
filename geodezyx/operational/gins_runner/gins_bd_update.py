#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27/02/2025 10:37:26

@author: psakic
"""


import argparse
import os
import subprocess
import datetime as dt
import re

from geodezyx import conv

#### Import the logger
import logging
log = logging.getLogger('geodezyx')

def read_credentials(file_path):
    with open(file_path, 'r') as file:
        login = file.readline().strip()
        password = file.readline().strip()
    return login, password

def download_rsync(file_list, remote_user, remote_host, remote_path, local_destination,
                   rsync_options=None, password=None):
    """
    Downloads a list of files using rsync.

    :param file_list: List of file paths to download (relative to the remote source).
    :param remote_user: Remote username (e.g., 'user').
    :param remote_host: Remote host (e.g., 'remote_host').
    :param remote_path: Remote path (e.g., '/path/to/source/').
    :param local_destination: Local destination directory to save the files.
    :param rsync_options: Optional list of additional rsync options (e.g., ['-avz', '--progress']).
    :param password: Optional password for the remote user.
    """
    if not rsync_options:
        rsync_options = ['-avz', '--relative']  # Default options: archive mode, verbose, and compression

    # Construct the remote source path
    remote_source = f"{remote_user}@{remote_host}:{remote_path}"

    if password:
        rsync_base = ['sshpass', '-p', password, 'rsync'] + rsync_options
    else:
        rsync_base = ['rsync'] + rsync_options

    for file in file_list:
        # Construct the rsync command
        # /./ : https://askubuntu.com/questions/552120/preserve-directory-tree-while-copying-with-rsync
        rsync_cmd = rsync_base + [f"{remote_source}/./{file}", f"{local_destination}/"]
        #log.info(rsync_cmd)

        # Run the rsync command
        result = subprocess.run(rsync_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Check if the command was successful
        if result.returncode != 0:
            log.warning(f"Download failed :( : {file}: {result.stderr}")
        else:
            if not re.search(file, result.stdout):
                log.info(f"File unchanged ;) : {file}")
            else:
                log.info(f"Download successful :) :  {file}")


def update_env_files(login, password, dir_bdgins):
    os.chdir(dir_bdgins)
    list_files = subprocess.check_output(
        "grep '/' /root/spotgins/metadata/directeur/DIR_SPOTGINS_G20_GE.yml | grep -v 'clock' | grep -v 'orbite' | grep -v 'ionex' | grep -v 'attitude' | awk -F: '{print $2}'",
        shell=True
    ).decode().splitlines()

    for fic in list_files:
        file = fic.split('/')[-1]
        rep = fic.replace(file, '')
        if not os.path.exists(rep):
            os.makedirs(rep)
        os.chdir(rep)
        commands = f"set xfer:clobber on\nmget /home/gins/MIROIR_STAF/{fic}\nexit"
        download_lftp(login, password, commands)
        os.chdir(dir_bdgins)

    if not os.path.exists('prairie'):
        os.makedirs('prairie')
    os.chdir('prairie')
    commands = "set xfer:clobber on\nmget /home/gins/MIROIR_STAF/prairie/igs_satellite_metadata.snx\nexit"
    download_lftp(login, password, commands)
    os.chdir(dir_bdgins)

def update_bdgins(date_srt, date_end, dir_bdgins, login ="", password =""):

    date = date_srt

    list_trop_file = []
    list_iono_file = []
    list_orbite_g20 = []
    list_orbex_g20 = []
    list_horl_g20 = []
    list_sp3_re3 = []

    while date <= date_end:
        day = str(date.day).zfill(2)
        month = str(date.month).zfill(2)
        year = date.year
        yy = str(year)[2:]

        wk, wkday = conv.dt2gpstime(date)
        doy = str(conv.dt2doy(date)).zfill(3)

        list_trop_file.extend([
            f"{year}/VMFG_{year}{month}{day}.H00",
            f"{year}/VMFG_{year}{month}{day}.H06",
            f"{year}/VMFG_{year}{month}{day}.H12",
            f"{year}/VMFG_{year}{month}{day}.H18"
        ])
        list_iono_file.append(f"{year}/igsg{doy}0.{yy}i.Z")
        list_orbite_g20.append(f"G20{wk}{wkday}.gin")
        list_orbex_g20.append(f"G20{wk}{wkday}.obx.gz")
        list_horl_g20.append(f"hogps_g20{wk}{wkday}")
        #list_sp3_re3.append(f"mg3{wk}{wkday}.sp3.Ci9PAU")
        date += dt.timedelta(days=1)

    dest_subdir_dic = {
        'tropo_vmf1': list_trop_file,
        'ionosphere/igs': list_iono_file,
        'mesures/gps/orbites/G20': list_orbite_g20,
        'mesures/gps/orbex/G20': list_orbex_g20,
        'mesures/gps/horloges30/G20': list_horl_g20,
        # 'orbites/SP3/re3': list_sp3_re3
    }

    create_dir(dir_bdgins, subdirs=dest_subdir_dic.keys())
    os.chdir(dir_bdgins)

    for subdir, files_list in dest_subdir_dic.items():
        log.info("Downloading files for %s", subdir)

        local_subdir_fullpath = os.path.join(dir_bdgins, subdir)
        remote_subdir_fullpath = os.path.join('/home/gins/MIROIR_STAF', subdir)

        download_rsync(files_list,
                       login,
                       'tite.get.obs-mip.fr',
                       remote_subdir_fullpath,
                       local_subdir_fullpath,
                       password)

def create_dir(parent_dir, subdirs):
    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)

    for subdir in subdirs:
        subdir_fullpath = os.path.join(parent_dir, subdir)
        if not os.path.exists(subdir_fullpath):
            os.makedirs(subdir_fullpath)

def download_fes2014b_files(login, password):
    base_dir = '/root/gin/data/OCELOAD/FES2014b'
    os.makedirs(base_dir, exist_ok=True)
    commands = "cd /home/gins/MIROIR_STAF/OCELOAD/FES2014b\nmget *\nquit"
    download_lftp(login, password, commands)

    base_dir = '/root/gin/data/OCELOAD/data'
    os.makedirs(base_dir, exist_ok=True)
    commands = "cd /home/gins/MIROIR_STAF/OCELOAD/FES2014b\nmget *\nquit"
    download_lftp(login, password, commands)






def main():
    parser = argparse.ArgumentParser(description="Update BDGINS repository.")
    parser.add_argument('-s', '--date_srt', type=lambda s: conv.date_pattern_2_dt(s), required=True, help="Start date in various format.")
    parser.add_argument('-e', '--date_end', type=lambda s: conv.date_pattern_2_dt(s), required=True, help="End date in various format.")
    parser.add_argument('-d', '--dir_bdgins',  help="Directory for BDGINS",  required=True)
    parser.add_argument('-l','--login', help="Login for the remote server.",  required=True)
    parser.add_argument('-p','--password', help="Password for the remote server.", required=False, default=None)
    args = parser.parse_args()

    update_bdgins(date_srt=args.date_srt, date_end=args.date_end,
                  login=args.login, password=args.password,
                  dir_bdgins=args.dir_bdgins)


if __name__ == "__main__":
    main()
