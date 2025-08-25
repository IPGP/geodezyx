#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic
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

import ftplib

#### Import the logger
import logging
import os
import pathlib
import shutil
import time
import urllib

########## BEGIN IMPORT ##########
#### External modules
from ftplib import FTP, FTP_TLS
from multiprocessing.dummy import Pool as ThreadPool

import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils
from geodezyx.stats import outlier_mad

log = logging.getLogger("geodezyx")


##########  END IMPORT  ##########


def start_end_date_easy(start_year, start_doy, end_year, end_doy):
    """
    generates start/end datetimes from a start/end year/day of year

    Parameters
    ----------
    start_year : int
        start year.
    start_doy : int
        start day of year.
    end_year : int
        end year.
    end_doy : int
        end day of year.

    Returns
    -------
    start : datetime
        converted start datetime.
    end : datetime
        converted end datetime.

    """
    start = conv.doy2dt(start_year, start_doy)
    end = conv.doy2dt(end_year, end_doy)
    return start, end


def effective_save_dir_orbit(
    parent_archive_dir, calc_center, date, archtype="year/doy/"
):
    """
    INTERNAL_FUNCTION

    archtype =
        stat
        stat/year
        stat/year/doy
        year/doy
        year/stat
        week/dow
        wkwwww : use a GFZ's CF-ORB wk<wwww> naming
        OR only '/' for a dirty saving in the parent folder
        ... etc ...
    """
    if archtype == "/":
        return parent_archive_dir

    out_save_dir = parent_archive_dir
    fff = archtype.split("/")
    year = str(date.year)
    doy = conv.dt2doy(date)
    week, dow = conv.dt2gpstime(date)

    for f in fff:
        if "wkwwww" in f:
            f_evaluated = "wk" + str(week).zfill(4)
        else:
            f_evaluated = eval(f)
        out_save_dir = os.path.join(out_save_dir, str(f_evaluated))
    return out_save_dir


#  _    _ _______ _______ _____    _____                      _                 _
# | |  | |__   __|__   __|  __ \  |  __ \                    | |               | |
# | |__| |  | |     | |  | |__) | | |  | | _____      ___ __ | | ___   __ _  __| |
# |  __  |  | |     | |  |  ___/  | |  | |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |
# | |  | |  | |     | |  | |      | |__| | (_) \ v  v /| | | | | (_) | (_| | (_| |
# |_|  |_|  |_|     |_|  |_|      |_____/ \___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|


#### HTTP classic Download


def downloader(
    url, savedir, force=False, check_if_file_already_exists_uncompressed=True
):
    """
    general function to download a file

    can also handle non secure FTP
    """

    if type(url) is tuple:
        need_auth = True
        username = url[1]
        password = url[2]
        url = url[0]
    else:
        need_auth = False
        username = ""
        password = ""

    url_print = str(url)

    rnxname = os.path.basename(url)

    pot_compress_files_list = [os.path.join(savedir, rnxname)]

    if check_if_file_already_exists_uncompressed:
        pot_compress_files_list.append(
            os.path.join(savedir, rnxname.replace(".gz", ""))
        )
        pot_compress_files_list.append(os.path.join(savedir, rnxname.replace(".Z", "")))
        pot_compress_files_list = list(set(pot_compress_files_list))

    for f in pot_compress_files_list:
        if os.path.isfile(f) and (not force):
            log.info(os.path.basename(f) + " already exists locally ;)")
            return None

    ##### LOCAL FILE (particular case for GFZ)
    if os.path.isfile(url):
        log.info("INFO : downloader : the is a local file, a simple copy will be used")
        log.info("       URL : %s", url)
        shutil.copy(url, savedir)

    ##### REMOTE FILE (General case)
    elif ("http" in url) or (("ftp" in url) and not need_auth):
        # managing an authentification
        if need_auth:  # HTTP with Auth
            password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
            top_level_url = url
            password_mgr.add_password(None, top_level_url, username, password)
            handler = urllib.request.HTTPBasicAuthHandler(password_mgr)
            # create "opener" (OpenerDirector instance)
            opener = urllib.request.build_opener(handler)
        else:  # FTP or HTTP without Auth
            opener = urllib.request.build_opener()

        # use the opener to fetch a URL
        try:
            f = opener.open(url)
        except (urllib.error.HTTPError, urllib.error.URLError) as exp:
            log.warning("%s not found :(", rnxname)
            log.warning(url_print)
            log.warning(exp)
            return ""

        log.info("%s found, downloading :)", rnxname)
        data = f.read()
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        outpath = os.path.join(savedir, rnxname)
        with open(outpath, "wb") as code:
            code.write(data)
        return_str = outpath

    elif ("ftp" in url) and need_auth:
        log.critical("MUST BE IMPEMENTED")
        return_str = ""
    else:
        log.error("something goes wrong with the URL")
        log.error(url)
        return_str = ""

    return return_str


def downloader_wrap(intup):
    downloader(*intup)
    return None


#  ______ _______ _____    _____                      _                 _
# |  ____|__   __|  __ \  |  __ \                    | |               | |
# | |__     | |  | |__) | | |  | | _____      ___ __ | | ___   __ _  __| |
# |  __|    | |  |  ___/  | |  | |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |
# | |       | |  | |      | |__| | (_) \ v  v /| | | | | (_) | (_| | (_| |
# |_|       |_|  |_|      |_____/ \___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|


#### FTP DOWNLOAD


class MyFTP_TLS(FTP_TLS):
    """
    This class is a subclass of FTP_TLS from the ftplib module. It is used to create an FTPS client that shares the TLS session.
    This is to avoid the error: ssl.SSLEOFError: EOF occurred in violation of protocol (_ssl.c:2396)
    Source: https://stackoverflow.com/questions/14659154/ftps-with-python-ftplib-session-reuse-required

    Methods
    -------
    ntransfercmd(cmd, rest=None)
        Initiate a data transfer over a new connection.
    """

    def ntransfercmd(self, cmd, rest=None):
        """
        Initiate a data transfer over a new connection.

        Parameters
        ----------
        cmd : str
            The command to send to the server.
        rest : str, optional
            A string that contains a marker representing where the server is to restart the operation's data transfer. Default is None.

        Returns
        -------
        tuple
            The connection and the expected size of the data.

        Notes
        -----
        If the protection level is set to private (i.e., _prot_p is True), the connection is wrapped in an SSL/TLS layer.
        """
        conn, size = FTP.ntransfercmd(self, cmd, rest)
        if self._prot_p:
            conn = self.context.wrap_socket(
                conn, server_hostname=self.host, session=self.sock.session
            )  # this is the fix
        return conn, size


def ftp_dir_list_files(ftp_obj_in):
    """
    Lists the files in the current directory of the FTP object.

    Parameters
    ----------
    ftp_obj_in : FTP object
        The FTP object used to connect to the FTP server.

    Returns
    -------
    list
        A list of filenames in the current directory of the FTP object.

    Notes
    -----
    This function tries to get the list of filenames in the current directory of the FTP object.
    If it encounters a permission error, it checks if the error message is "550 No files found".
    If it is, it logs a warning message and returns an empty list.
    If the error message is different, it raises the exception.
    """
    files = []
    try:
        files = ftp_obj_in.nlst()
    except ftplib.error_perm as resp:
        if str(resp) == "550 No files found":
            log.warning("No files in this directory" + ftp_obj_in.pwd())
        else:
            raise
    return files


def ftp_objt_create(
    secure_ftp_inp=False,
    host="",
    chdir="",
    parallel_download=1,
    user="anonymous",
    passwd="",
    retry_count=3,
):
    """
    This function creates and returns an FTP object and a list of FTP objects for multiple downloads.

    Parameters
    ----------
    secure_ftp_inp : bool
        If True, uses FTPS for secure file transfer. Default is False.
    host : str, optional
        The hostname of the FTP server. Default is an empty string.
    chdir : str, optional
        The directory to change to after connecting to the FTP server. Default is an empty string.
    parallel_download : int, optional
        The number of parallel downloads to be performed. Default is 1.
    user : str, optional
        The username for the FTP server. Default is "anonymous".
    passwd : str, optional
        The password for the FTP server. Default is an empty string.
    retry_count : int, optional
        The number of times to retry creating the FTP object. Default is 3.

    Returns
    -------
    tuple
        The main FTP object for crawling and a list of FTP objects for parallel downloads.

    Notes
    -----
    This function creates an FTP object using the appropriate constructor based on the secure_ftp_inp parameter.
    It then creates a list of FTP objects for multiple downloads.
    If a directory is specified, it changes the current working directory of the main FTP object to that directory.
    """

    # define the right constructor
    if secure_ftp_inp:
        ftp_constuctor = MyFTP_TLS
    else:
        ftp_constuctor = FTP

    # create a list of FTP object for multiple downloads
    ftp_obj_list_out = []
    for i in range(parallel_download):
        for attempt in range(retry_count):
            try:
                current_ftp_obj = ftp_constuctor(host)
                ftp_obj_list_out.append(current_ftp_obj)
                break  # Exit the retry loop if successful
            except Exception as e:
                log.warning("FTP object creation failed on attempt %d", attempt + 1)
                log.warning(e)
                time.sleep(5)
                if attempt == retry_count - 1:
                    log.error("Max retries reached. Could not create FTP object.")

    if secure_ftp_inp:
        [f.login(user, passwd) for f in ftp_obj_list_out]
        [f.prot_p() for f in ftp_obj_list_out]
    else:
        [f.login(user=user, passwd=passwd) for f in ftp_obj_list_out]

    # define the main obj for crawling
    ftp_main = ftp_obj_list_out[0]

    # change the directory of the main ftp obj if we ask for it
    if chdir:
        log.info("Move to: %s", chdir)
        ftp_main.cwd(chdir)

    return ftp_main, ftp_obj_list_out


def ftp_downld_core(ftp_obj, filename, localdir, force=False):
    """
    Performs the FTP download if we are already in the correct FTP folder.
    This is an internal function of ftp_downld.

    Parameters
    ----------
    ftp_obj : FTP object
        The FTP object used to connect to the FTP server.
    filename : str
        The name of the file to be downloaded.
    localdir : str
        The local directory where the downloaded file should be saved.
    force : bool, optional
        If True, forces the download even if the file already exists locally.
        Default is False.

    Returns
    -------
    tuple
        The local path of the downloaded file and a boolean indicating whether the download was successful.

    Notes
    -----
    This function first checks if the local directory exists, if not it creates it.
    Then it checks if the file already exists locally, if it does, it logs a message and returns.
    If the file does not exist, it tries to download the file from the FTP server.
    If the download is successful, it logs a success message and returns the local path and True.
    If the download fails, it logs a failure message and returns the local path and False.
    """

    localpath = os.path.join(localdir, filename)

    if not os.path.isdir(localdir):
        utils.create_dir(localdir)

    dl_go = True
    # Check if the file already exists locally
    bool_dl = False
    if not utils.empty_file_check(localpath):
        if not force:
            log.info(filename + " already exists ;)")
            bool_dl = True
            dl_go = False
        else:
            log.info(filename + " already exists, but re-download forced")
            bool_dl = False
            dl_go = True

    if dl_go:
        try:
            localfile = open(localpath, "wb")
            ftp_obj.retrbinary("RETR " + filename, localfile.write, 1024)
            localfile.close()
            bool_dl = True
            log.info(filename + " downloaded :)")

        except Exception as e:
            log.warning(localpath + " download failed :(")
            log.warning(e)
            bool_dl = False

    return localpath, bool_dl


def ftp_downld_mono(ftp_obj, full_remote_path, localdir, force=False):
    """
    Downloads a file through FTP protocol.

    Parameters
    ----------
    ftp_obj : FTP object
        The FTP object used to connect to the FTP server.
    full_remote_path : str
        The full path of the file on the FTP server.
    localdir : str
        The local directory where the downloaded file should be saved.
    force : bool, optional
        If True, forces the download even if the file already exists locally.
        Default is False.

    Returns
    -------
    tuple
        The output of the ftp_downloader_core function.

    Notes
    -----
    This function changes the current working directory of the FTP object to the directory of the file to be downloaded,
    and then calls the ftp_downloader_core function to download the file.
    """

    filename = os.path.basename(full_remote_path)
    intermed_path = full_remote_path.split("/")[3:]
    intermed_path.remove(filename)
    intermed_path = "/" + "/".join(intermed_path)

    ftp_obj.cwd(intermed_path)

    return ftp_downld_core(ftp_obj, filename, localdir, force=force)


def ftp_downld_wrap(intup):
    """
    This function is a wrapper for the ftp_downld function. It unpacks the input tuple and passes it to the ftp_downld function.

    Parameters
    ----------
    intup : tuple
        A tuple containing the parameters to be passed to the ftp_downld function.

    Returns
    -------
    tuple
        The output of the ftp_downld function.
    """
    outtup = ftp_downld_mono(*intup)
    return outtup


def ftp_downld_front(
    urls,
    savedirs,
    parallel_download=1,
    secure_ftp=False,
    user="anonymous",
    passwd="anonymous@isp.com",
    force=True,
):
    """
    This function is used to download files from FTP servers in parallel.

    Parameters
    ----------
    urls : str or list
        The URL or list of URLs of the files to be downloaded.
    savedirs : str or list
        The directory or list of directories where the downloaded files should be saved.
    parallel_download : int, optional
        The number of parallel downloads to be performed. Default is 1.
    secure_ftp : bool, optional
        If True, uses FTPS for secure file transfer. Default is False.
    user : str, optional
        The username for the FTP server. Default is "anonymous".
    passwd : str, optional
        The password for the FTP server. Default is 'anonymous@isp.com'.
    force : bool, optional
        If True, forces the download even if the file already exists. Default is True.

    Returns
    -------
    out_tup_lis : List of tuples
        Returns a list of tuples containing
        the local path of the downloaded file and
        a boolean indicating whether the download was successful.
        e.g. [(local_path1, True), (local_path2, False), ...]

    Notes
    -----
    This function uses the ThreadPool for parallel downloads.
    """

    # Check if urls and savedirs are iterable, if not convert them to list
    urllist = urls if utils.is_iterable(urls) else [urls]
    savedirlist = savedirs if utils.is_iterable(savedirs) else [savedirs] * len(urllist)
    secure_ftp_use = secure_ftp[0] if utils.is_iterable(secure_ftp) else secure_ftp
    ##### dirty to select secure_ftp 1st elt only....

    # Check if the length of urllist and savedirlist are the same
    if len(urllist) != len(savedirlist):
        log.error(
            "URL & out dirs lists do not have the same length: %s %s",
            len(urllist),
            len(savedirlist),
        )

    # Extract the host from the urls
    urlpathobj = pd.Series(urllist).apply(pathlib.Path)
    host_use = urlpathobj.apply(lambda p: p.parts[1]).unique()[0]

    # Create the FTP object
    ftpobj_main, ftpobj_lis = ftp_objt_create(
        secure_ftp_inp=secure_ftp_use,
        host=host_use,
        parallel_download=parallel_download,
        user=user,
        passwd=passwd,
    )

    # Create a list of FTP objects for parallel downloads
    ftpobj_mp_lis = ftpobj_lis * int(np.ceil(len(urllist) / parallel_download))
    force_lis = [force] * len(urllist)

    # Check if there are less FTP objects than URLs for parallel download
    if len(ftpobj_mp_lis) < len(urllist):
        log.warning(
            "less FTP objects than URL for parallel download, contact the main developper"
        )

    # Create a ThreadPool for parallel downloads
    pool = ThreadPool(parallel_download)

    # Start the parallel downloads
    out_tup_lis = pool.map(
        ftp_downld_wrap, list(zip(ftpobj_mp_lis, urllist, savedirlist, force_lis))
    )

    return out_tup_lis


def ftp_downloader_wo_objects(tupin):
    """
    create the necessary FTP object

    should not be used anymore
    """
    arch_center_main, wwww_dir, filename, localdir = tupin
    ftp_obj_wk = FTP(arch_center_main)
    ftp_obj_wk.login()
    ftp_obj_wk.cwd(wwww_dir)
    localpath, bool_dl = ftp_downld_core(ftp_obj_wk, filename, localdir)
    ftp_obj_wk.close()
    return localpath, bool_dl


def ftp_files_crawler_legacy(urllist, savedirlist, secure_ftp):
    """
    filter urllist,savedirlist generated with download_gnss_rinex with an
    optimized FTP crawl

    """
    ### create a DataFrame based on the urllist and savedirlist lists
    df = pd.concat((pd.DataFrame(urllist), pd.DataFrame(savedirlist)), axis=1)
    df_orig = df.copy()

    ### rename the columns
    if df.shape[1] == 4:
        loginftp = True
        df.columns = ("url", "user", "pass", "savedir")
    else:
        loginftp = False
        df.columns = ("url", "savedir")
        df["user"] = "anonymous"
        df["pass"] = ""

    ### Do the correct split for the URLs
    df = df.sort_values("url")
    df["url"] = df["url"].str.replace("ftp://", "")
    df["dirname"] = df["url"].apply(os.path.dirname)
    df["basename"] = df["url"].apply(os.path.basename)
    df["root"] = [e.split("/")[0] for e in df["dirname"].values]
    df["dir"] = [e1.replace(e2, "")[1:] for (e1, e2) in zip(df["dirname"], df["root"])]
    df["bool"] = False

    #### Initialisation of the 1st variables for the loop
    prev_row_ftpobj = df.iloc[0]
    prev_row_cwd = df.iloc[0]
    ftp_files_list = []
    count_loop = 0  # restablish the connexion after 50 loops (avoid freezing)
    #### Initialisation of the FTP object

    ftpobj, _ = ftp_objt_create(
        secure_ftp_inp=secure_ftp,
        host=prev_row_ftpobj.root,
        user=prev_row_ftpobj.user,
        passwd=prev_row_ftpobj["pass"],
    )

    for irow, row in df.iterrows():
        count_loop += 1

        ####### we recreate a new FTP object if the root URL is not the same
        if row.root != prev_row_ftpobj.root or count_loop > 20:
            ftpobj, _ = ftp_objt_create(
                secure_ftp_inp=secure_ftp,
                host=prev_row_ftpobj.root,
                user=prev_row_ftpobj.user,
                passwd=prev_row_ftpobj["pass"],
            )

            prev_row_ftpobj = row
            count_loop = 0

        ####### we recreate a new file list if the date path is not the same
        if (prev_row_cwd.dir != row.dir) or irow == 0:
            log.info("chdir " + row.dirname)
            ftpobj.cwd("/")

            try:  #### we try to change for the right folder
                ftpobj.cwd(row.dir)
            except:  #### If not possible, then no file in the list
                ftp_files_list = []

            ftp_files_list = ftp_dir_list_files(ftpobj)
            prev_row_cwd = row

            ####### we check if the files is avaiable
        if row.basename in ftp_files_list:
            df.loc[irow, "bool"] = True
            log.info(row.basename + " found on server :)")
        else:
            df.loc[irow, "bool"] = False
            log.warning(row.basename + " not found on server :(")

    df_good = df[df["bool"]].copy()

    df_good["url"] = "ftp://" + df_good["url"]

    ### generate the outputs
    if loginftp:
        urllist_out = list(zip(df_good.url, df_good.user, df_good["pass"]))
    else:
        urllist_out = list(df_good.url)

    savedirlist_out = list(df_good.savedir)

    return urllist_out, savedirlist_out
