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
import sys

from geodezyx import conv, utils

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
    """
    if not rsync_options:
        rsync_options = ['-avz', '--progress', '--relative', '--copy-links']  # Default options: archive mode, verbose, and compression

    # Construct the remote source path
    remote_source = f"{remote_user}@{remote_host}:{remote_path}"

    if password:
        rsync_base = ['sshpass', '-p', password, 'rsync'] + rsync_options
    else:
        rsync_base = ['rsync'] + rsync_options

    tmp_rsync_file_lis = f"/tmp/" + utils.get_timestamp() + "_tmp_rsync_file_list.lst"
    utils.write_in_file("\n".join(file_list), tmp_rsync_file_lis)

    # Construct the rsync command
    # /./ : https://askubuntu.com/questions/552120/preserve-directory-tree-while-copying-with-rsync
    rsync_cmd = rsync_base + ['--files-from', tmp_rsync_file_lis] + [f"{remote_source}/./", f"{local_destination}/"]

    # Run the rsync command
    process = subprocess.run(rsync_cmd, stderr=sys.stderr, stdout=sys.stdout, text=True)

    if process.returncode != 0:
        log.error("Rsync failed with return code :( %d", process.returncode)
    else:
        log.info("Rsync completed successfully :)")

    os.remove(tmp_rsync_file_lis)
    return

def update_bdgins(date_srt, date_end, dir_bdgins, login ="", password =""):

    date = date_srt

    ###### LIST INITIALISATION
    ### misc files do not need initialisation
    list_misc = []

    ### time dependant files
    list_tropo = []
    list_iono = []
    list_orbite_g20 = []
    list_orbex_g20 = []
    list_horl_g20 = []
    list_sp3_re3 = []

    ###### LIST FILL
    ### misc files
    # list_mda_snx_pra = [f"igs_satellite_metadata.snx"]
    # list_pole = [f'nominal_NRO']
    # list_antex = [f'igs20.atx']
    # list_valap_static = [f'valap_static']
    # list_macromod = [f'gnss.xml']
    # list_lunisol = [f'de440bdlf.ad']
    # list_maree_pol = [f'nominal']

    list_misc.extend([f"prairie/igs_satellite_metadata.snx"])
    list_misc.extend([f'pole/nominal_NRO'])
    list_misc.extend([f'ANTEX/igs20.atx'])
    list_misc.extend([f'EXE_PPP/valap_static'])
    list_misc.extend([f'macromod/gnss.xml'])
    list_misc.extend([f'lunisolaires/de440bdlf.ad'])
    list_misc.extend([f'maree_polaire/loading/nominal'])

    ### time dependant files
    while date <= date_end:
        day = str(date.day).zfill(2)
        month = str(date.month).zfill(2)
        year = date.year
        yy = str(year)[2:]

        wk, wkday = conv.dt2gpstime(date)
        doy = str(conv.dt2doy(date)).zfill(3)

        list_tropo.append(f"orography_ell")
        list_tropo.extend([
            f"{year}/VMFG_{year}{month}{day}.H00",
            f"{year}/VMFG_{year}{month}{day}.H06",
            f"{year}/VMFG_{year}{month}{day}.H12",
            f"{year}/VMFG_{year}{month}{day}.H18"
        ])
        list_iono.append(f"{year}/igsg{doy}0.{yy}i.Z")
        list_orbite_g20.append(f"G20{wk}{wkday}.gin")
        list_orbex_g20.append(f"G20{wk}{wkday}.obx.gz")
        list_horl_g20.append(f"hogps_g20{wk}{wkday}")
        #list_sp3_re3.append(f"mg3{wk}{wkday}.sp3.Ci9PAU")
        date += dt.timedelta(days=1)


    ###### DESTINATION FOLDERS

    dest_subdir_dic = {
        ### misc files
        # 'prairie': list_mda_snx_pra,
        # 'pole': list_pole,
        # 'ANTEX': list_antex,
        # 'EXE_PPP': list_valap_static,
        # 'macromod': list_macromod,
        # 'lunisolaires': list_lunisol,
        # 'maree_polaire/loading': list_maree_pol,
        '.': list_misc,

        ### time dependant files
        'tropo_vmf1': list_tropo,
        'ionosphere/igs': list_iono,
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

    # compress the clock files
    for dirr in ['mesures/gps/horloges30/G20', 'mesures/gps/orbites/G20']:
        comp_clk_cmd = ['find', os.path.join(dir_bdgins, dirr), '-name', '!', '.gz', '-exec', 'gzip', '-v', '{}', '\;']
        process = subprocess.run(comp_clk_cmd, stderr=sys.stderr, stdout=sys.stdout, text=True)


    return None



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
