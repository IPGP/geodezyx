#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mansur
This sub-module of geodezyx.operational contains functions to download
gnss data and products from distant IGS servers. 
it can be imported directly with:
from geodezyx import operational
The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License
Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

#### Import the logger
import logging
########## BEGIN IMPORT ##########
#### External modules
import os
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

log = logging.getLogger('geodezyx')

##########  END IMPORT  ##########


## GUS CDDIS FCT FOR REPRO3 ONLY!!!! 2022-10-10
def list_file_in_ftp(
    email="mansur@gfz-potsdam.de",
    host="gdc.cddis.eosdis.nasa.gov",
    directory="gnss/products/repro3/1900",
):
    """
    List files inside one directory in a ftp serve. Note that this uses an anonymous user.

    Parameters
    ----------
    email : login email
    host  : ftp server
    directory  : directory inside the server
    Returns
    -------
    Return a list with the the information about the directory.

    """
    from ftplib import FTP_TLS

    # import sys
    # import ftplib
    ftps = FTP_TLS(host=host)
    ftps.login(user="anonymous", passwd=email)
    ftps.prot_p()
    ftps.cwd(directory)
    dir_list = []
    # print(ftps.dir(dir_list.append))
    ftps.dir(dir_list.append)

    return dir_list


def list_daily_files_ftp(ftp_dir_list=None, fmt="sp3", ac="COD", week=1900, dayw=0):
    """
    List files inside one directory in a ftp serve. Note that this uses an anonymous user.

    Parameters
    ----------
    ftp_dir_list : list with the files in the ftp (comming from the previous fct: list_file_in_ftp)
    fmt : searching format file e.g. sp3
    ac : name of the anlisys center
    week : gps week
    day : day of the week
    Returns
    -------
    Return a list with the the information about the directory.

    """

    ftp_dir_list = ftp_dir_list

    ## List all the files important e.g. SP3, BIA, SNX etc
    files_list = []
    for item in ftp_dir_list:
        splt_item = item.split()
        if splt_item[-1][-2:] not in ["gz"]:
            continue
        else:
            files_list.append(splt_item[-1])
    ##

    ## get all the format needed
    ## sp3
    if fmt in ["SP3", "EPH", "sp3", "eph"]:
        fmt = ["SP3", "EPH", "sp3", "eph"]
    ## clk
    if fmt in ["CLK", "clk"]:
        fmt = ["CLK", "clk"]
    fmt_list = []
    for item in files_list:
        if item[-6:-3] in fmt:
            fmt_list.append(item)

    ##

    ## get the desired AC
    ac = ac
    ac_list = []
    for item in fmt_list:
        if item[0:3] == ac:
            ac_list.append(item)
    ##

    ## Find the daily file
    week = week
    dayw = dayw
    day, year = conv.dt2doy_year(conv.gpstime2dt(week, dayw))
    string_date = str(year) + str(day)
    day_file_lst = []
    counter = 0
    for item in ac_list:
        if string_date in item:
            day_file_lst.append(item)
            counter = counter + 1
    if not day_file_lst:
        return None
    else:
        if counter > 1:
            print("WRN: More than one file found for this AC for this day!!")
            return day_file_lst
        else:
            return day_file_lst[0]
    ##


def download_rp3_cddis(week_ini=None, ac=None, fmt=None, out_path=None):
    """
    Download sp3 weekly files for repro3 products, from cddis (NASA).
    Note: For the moment is only for sp3 and clk.

    Parameters
    ----------
    week_ini : GPS week for download
    fmt : searching format file e.g. sp3
    ac : name of the anlisys center
    out_path : path where the products will be download

    Returns
    -------
    No returns.

    """
    ## arguments from the calling fct
    week_ini = week_ini
    week_final = week_ini + 1
    ac = ac
    fmt = fmt
    ##
    ## Default configuration using the user and cddis from NASA, on this case to repro3
    email = "mansur@gfz-potsdam.de"
    host = "gdc.cddis.eosdis.nasa.gov"
    directory = "gnss/products/repro3/"
    ##
    ## Set the output path. Not that os.chdir do a cd to the path. So
    ## the products are download to the right location
    out_path = out_path
    results_dir = os.path.join(out_path, "Week_" + str(week_ini))
    utils.create_dir(results_dir)
    os.chdir(results_dir)
    ##
    ## creates a list with the files for the week
    dir_wk = directory + str(week_ini)
    ftp_dir_list = list_file_in_ftp(email=email, host=host, directory=dir_wk)
    ##

    ## Do a loop for each day in the week
    for week in range(week_ini, week_final):
        for day in range(0, 7):
            try:
                file = list_daily_files_ftp(
                    ftp_dir_list=ftp_dir_list, fmt=fmt, ac=ac, week=week, dayw=day
                )
                if not file:
                    print("WRN: no file found!")
                    continue
                ## check if file is a list. This is the case when the AC has
                ## more than one solution for the day!! e.g. MIT
                if isinstance(file, list):
                    for fl in file:
                        cmd = (
                            "wget -qr -nH --cut-dirs=9 -t 10 --ftp-user anonymous --ftp-password "
                            + email
                            + " ftps://"
                            + host
                            + "/"
                            + dir_wk
                            + "/"
                            + fl
                        )
                        print(cmd)
                        os.system(cmd)
                else:
                    cmd = (
                        "wget -qr -nH --cut-dirs=9 -t 10 --ftp-user anonymous --ftp-password "
                        + email
                        + " ftps://"
                        + host
                        + "/"
                        + dir_wk
                        + "/"
                        + file
                    )
                    print(cmd)
                    os.system(cmd)

            except Exception as e:
                print(
                    "ERR: THERE IS AN ERROR WITH THIS DATE -- WEEK: "
                    + str(week)
                    + " DAY "
                    + str(day)
                )
                print(e)
                elm_stk = [week, day, e]
                err = pd.DataFrame(elm_stk)
                err = err.T
                err = err.rename(columns={0: "week", 1: "day", 2: "error"})
                pickle_name = "ERROR_" + str(week_ini) + "_" + str(day)
                err.to_pickle(results_dir + "/" + pickle_name + ".pkl")
                continue

    return None
