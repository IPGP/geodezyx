#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic
This sub-module of geodezyx.operational contains functions to download
gnss data and products from distant IGS servers.
it can be imported directly with:
from geodezyx import operational
The geodezyx toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License
Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/IPGP/geodezyx
"""

import ftplib
import requests
import tqdm

#### Import the logger
import logging
import os
import pathlib
import shutil
import time
import urllib
import urllib.request

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


def download_http(url, output_dir, timeout=120, max_try=4, sleep_time=5):
    """
    Download a file from an HTTP server with retry logic and progress bar.

    Parameters
    ----------
    url : str
        The URL of the file to download.
    output_dir : str
        The directory where the downloaded file will be saved.
    timeout : int, optional
        The timeout for the HTTP connection in seconds. Default is 120 seconds.
    max_try : int, optional
        The maximum number of retry attempts in case of failure. Default is 4.
    sleep_time : int, optional
        The sleep time between retry attempts in seconds. Default is 5 seconds.

    Returns
    -------
    str
        The path to the downloaded file, or an empty string if the download failed.

    Raises
    ------
    AutorinoDownloadError
        If the download fails after the maximum number of retry attempts.
    """

    # Get file size
    log.info("Download file: %s", url)
    response = requests.head(url, timeout=timeout)
    file_size = int(response.headers.get("content-length", 0))

    # Construct output path
    filename = url.split("/")[-1]
    output_path = os.path.join(output_dir, filename)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    dwl = False

    # Check if file already exists
    if os.path.isfile(output_path):
        log.info(f"{filename} already exists locally ;)")
        return (output_path, dwl)

    # Download file with progress bar
    try_count = 0
    while True:
        try:
            response = requests.get(url, stream=True, timeout=timeout)
            with open(output_path, "wb") as f:
                with tqdm.tqdm(
                    total=file_size, unit="B", unit_scale=True, desc=filename
                ) as pbar:
                    for data in response.iter_content(chunk_size=1024):
                        f.write(data)
                        pbar.update(len(data))
            break
        except requests.exceptions.RequestException as e:
            try_count += 1
            if try_count > max_try:
                log.error("download failed after %i attempts", max_try)
                return (url, dwl)

            log.warning("download failed (%s), try %i/%i", str(e), try_count, max_try)
            time.sleep(sleep_time)

    dwl = True
    return (output_path, dwl)


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


######################################################################
######## GENERIC FTP CRAWL INFRASTRUCTURE
######################################################################
# These functions are used by both download_rinex and download_prods
# to crawl FTP directories, check local files, and manage downloads
# in a unified table-based approach.


def regex_match_indir(regex, dir_files_list):
    """
    Match files in a directory against a given regex pattern.

    Parameters
    ----------
    regex : str
        Regex pattern to match filenames.
    dir_files_list : list of str
        List of filenames in the directory.

    Returns
    -------
    str or None
        The first matching filename, or None if no match is found.
    """
    import re

    matches = [file for file in dir_files_list if re.search(regex, file)]
    return matches[0] if matches else None


def regex_match_indir_all(regex, dir_files_list):
    """
    Match all files in a directory against a given regex pattern.

    Parameters
    ----------
    regex : str
        Regex pattern to match filenames.
    dir_files_list : list of str
        List of filenames in the directory.

    Returns
    -------
    list of str
        List of all matching filenames, or empty list if no match is found.
    """
    import re

    matches = [file for file in dir_files_list if re.search(regex, file)]
    return matches


def check_local_file_exists(
    filergx, outdir, local_files_cache, force=False, all_files_mode=False
):
    """
    Check if a file matching a regex exists locally.

    Parameters
    ----------
    filergx : str
        Regex pattern to match filenames.
    outdir : str
        Local output directory path.
    local_files_cache : list of str
        Cached list of local files in the directory.
    force : bool, optional
        Force re-download even if file exists. Default is False.
    all_files_mode : bool, optional
        If True, check for all matching files. Default is False.

    Returns
    -------
    tuple of (bool, bool, str)
        - ok_loc : True if file exists locally and should not be re-downloaded
        - ok_dwl : True if file should be downloaded (not local or force=True)
        - filename : Filename if found locally (semicolon-separated if
          all_files_mode=True), empty string otherwise
    """
    if all_files_mode:
        fillocal_list = regex_match_indir_all(filergx, local_files_cache)

        if fillocal_list:
            valid_files = [f for f in fillocal_list if os.path.getsize(f) > 0]

            if valid_files:
                fil_bns = [os.path.basename(f) for f in valid_files]
                filename = ";".join(fil_bns)

                if not force:
                    log.info("%d file(s) already exist locally ;)", len(valid_files))
                    return True, False, filename
                else:
                    log.info(
                        "%d file(s) already exist locally, but re-download forced",
                        len(valid_files),
                    )
                    return False, True, filename
    else:
        fillocal = regex_match_indir(filergx, local_files_cache)

        if fillocal and os.path.getsize(fillocal) > 0:
            fil_bn = os.path.basename(fillocal)
            if not force:
                log.info("%s already exists locally ;)", fil_bn)
                return True, False, fil_bn
            else:
                log.info("%s already exists locally, but re-download forced", fil_bn)
                return False, True, fil_bn

    return False, False, ""


def get_ftp_connection(
    ftpobj, host, protocol, sftp, user, passwd, prev_host, count_loop, count_nmax
):
    """
    Get or create an FTP connection.

    Creates a new connection if the host has changed, counter exceeds max,
    or this is the first connection. Otherwise returns the existing connection.

    Parameters
    ----------
    ftpobj : FTP or None
        Current FTP connection object.
    host : str
        Target FTP server hostname.
    protocol : str
        Protocol string ("ftp" or "sftp").
    sftp : str or bool
        SFTP mode setting ('auto', True, or False).
    user : str or None
        FTP username.
    passwd : str or None
        FTP password.
    prev_host : str
        Previous host name.
    count_loop : int
        Current operation counter.
    count_nmax : int
        Maximum operations before reconnection.

    Returns
    -------
    tuple of (FTP, str)
        - ftpobj : FTP connection object
        - new_host : str - The host name for connection tracking
    """
    if host != prev_host or count_loop > count_nmax or count_loop == 1:
        if ftpobj:
            ftpobj.close()

        if sftp == "auto":
            sftp_use = protocol == "sftp"
        else:
            sftp_use = bool(sftp)

        ftpobj, _ = ftp_objt_create(
            secure_ftp_inp=sftp_use,
            host=host,
            user=user,
            passwd=passwd,
        )
        return ftpobj, host

    return ftpobj, prev_host


def get_ftp_directory_listing(ftpobj, directory, host, prev_dir):
    """
    Get FTP directory file listing.

    Changes to the specified directory and retrieves the list of files.
    Returns empty list if directory change fails.

    Parameters
    ----------
    ftpobj : FTP
        FTP connection object.
    directory : str
        Remote directory path.
    host : str
        FTP server hostname (for logging/URL generation).
    prev_dir : str
        Previous directory path.

    Returns
    -------
    tuple of (list, list)
        - ftp_files_list : list of filenames in the directory
        - ftp_files_urls : list of complete FTP URLs for all files
        Returns (None, None) if directory is unchanged.
    """
    if prev_dir != directory:
        log.info("chdir " + directory)
        ftpobj.cwd("/")

        try:
            ftpobj.cwd(directory)
            ftp_files_list = ftp_dir_list_files(ftpobj)
            ftp_files_urls = [f"ftp://{host}{directory}/{f}" for f in ftp_files_list]
            return ftp_files_list, ftp_files_urls

        except Exception as e:
            log.warning("unable to chdir to %s, exception %s", directory, e)
            return [], []

    return None, None


def match_files_in_directory(filergx, ftp_files_list, all_files_mode=False):
    """
    Match files in FTP directory using regex pattern.

    Parameters
    ----------
    filergx : str
        Regex pattern to match filenames.
    ftp_files_list : list of str
        List of filenames in the FTP directory.
    all_files_mode : bool, optional
        If True, return all matching files. Default is False.

    Returns
    -------
    str or list of str or None
        Matched filename(s).
    """
    if all_files_mode:
        return regex_match_indir_all(filergx, ftp_files_list)
    else:
        return regex_match_indir(filergx, ftp_files_list)


def update_table_row_with_match(table, irow, file_match, all_files_mode=False):
    """
    Update table row based on file match results.

    Parameters
    ----------
    table : pd.DataFrame
        Table to update.
    irow : int
        Row index to update.
    file_match : str or list of str or None
        Matched filename(s) from FTP directory.
    all_files_mode : bool, optional
        If True, file_match is a list of files. Default is False.
    """
    if file_match:
        table.loc[irow, "ok_dwl"] = True
        if all_files_mode:
            table.loc[irow, "filnam"] = ";".join(file_match)
            log.info(f"{len(file_match)} file(s) found on server :)")
        else:
            table.loc[irow, "filnam"] = file_match
            log.info(file_match + " found on server :)")
    else:
        table.loc[irow, "ok_dwl"] = False
        table.loc[irow, "filnam"] = ""
        log.warning(f"{table.loc[irow, 'filrgx']} not found on server :(")

    table.loc[irow, "crawled"] = True


def generate_download_urls(table, all_files_mode=False):
    """
    Generate complete FTP URLs for files to be downloaded.

    Parameters
    ----------
    table : pd.DataFrame
        Table with crawl results. Must have columns: ok_dwl, ok_loc, host,
        dir, filnam (or rnxnam for backward compat).
    all_files_mode : bool, optional
        If True, handle semicolon-separated filenames. Default is False.
    """
    # support both new generic column name and legacy rnxnam
    nam_col = "filnam" if "filnam" in table.columns else "rnxnam"

    fil_ok_dwl = table["ok_dwl"] & ~table["ok_loc"]

    if all_files_mode:

        def make_urls(row):
            filenames = row[nam_col].split(";")
            # Normalize dir: RINEX dir has no leading '/' (built from pathlib.Path parts),
            # while products dir has a leading '/' (built from os.path.join with absolute
            # basedir).  Strip any leading '/' and always add exactly one.
            dir_clean = row["dir"].lstrip("/")
            urls = [f"ftp://{row['host']}/{dir_clean}/{fname}" for fname in filenames]
            return ";".join(urls)

        table.loc[fil_ok_dwl, "url_true"] = table.loc[fil_ok_dwl].apply(
            make_urls, axis=1
        )
    else:
        # Same normalization for the single-file case
        table.loc[fil_ok_dwl, "url_true"] = table.loc[fil_ok_dwl].apply(
            lambda x: f"ftp://{x['host']}/{x['dir'].lstrip('/')}/{x[nam_col]}", axis=1
        )


def collect_local_files(table):
    """
    Collect local file paths for files that exist locally.

    Parameters
    ----------
    table : pd.DataFrame
        Table with crawl results containing 'ok_loc', 'outdir' columns,
        and either 'filnam' or 'rnxnam' column.

    Returns
    -------
    pd.Series
        Series of local file paths, or empty Series if no local files.
    """
    nam_col = "filnam" if "filnam" in table.columns else "rnxnam"

    ok_loc = table["ok_loc"]
    if ok_loc.sum() > 0:
        return pd.Series(
            table.loc[ok_loc].apply(
                lambda e: os.path.join(e["outdir"], e[nam_col]), axis=1
            )
        )
    else:
        return pd.Series([])


def crawl_ftp_files(
    table,
    sftp="auto",
    user=None,
    passwd=None,
    path_ftp_crawled_files_save=None,
    path_all_ftp_files_save=None,
    force=False,
    all_files_mode=False,
    exclude_patterns=None,
):
    """
    Crawl FTP servers to find available files and update download table.

    This is a generic function used by both RINEX and products downloaders.
    It checks for existing local files, connects to FTP servers, lists remote
    files, and updates the table with availability status and actual file URLs.

    Parameters
    ----------
    table : pd.DataFrame
        Input table containing download metadata with columns:
        - 'host': FTP server hostname
        - 'dir': Remote directory path
        - 'outdir': Local output directory
        - 'filrgx' or 'rnxrgx': Filename regex pattern
        - 'protocol': String indicating protocol ("ftp" or "sftp")
        - 'crawled': Boolean indicating if already crawled
    sftp : str or bool, optional
        SFTP mode setting. Default is 'auto'.
    user : str, optional
        FTP username. Default is None (anonymous).
    passwd : str, optional
        FTP password. Default is None (anonymous).
    path_ftp_crawled_files_save : str, optional
        Path to save the crawled files table as CSV.
    path_all_ftp_files_save : str, optional
        Path to save all discovered FTP files as CSV.
    force : bool, optional
        Force re-download even if files exist locally. Default is False.
    all_files_mode : bool, optional
        If True, download ALL files matching the regex.
        Default is False.
    exclude_patterns : list of str, optional
        List of regex patterns to exclude from matched files.
        Default is None.

    Returns
    -------
    tuple of (pd.DataFrame, pd.Series, pd.Series)
        - table_use : Updated table with crawl results
        - all_ftp_files : All files discovered on FTP servers
        - all_loc_files : Local file paths for files that already exist
    """
    import glob
    import re as re_module

    # Support both generic (filrgx/filnam) and legacy (rnxrgx/rnxnam) column names
    rgx_col = "filrgx" if "filrgx" in table.columns else "rnxrgx"
    nam_col = "filnam" if "filrgx" in table.columns else "rnxnam"

    def _save_crawled_files(table_inp):
        if path_ftp_crawled_files_save:
            table_inp.to_csv(path_ftp_crawled_files_save)

    def _get_and_save_all_ftp_files(all_ftp_files_stk_inp):
        if all_ftp_files_stk_inp:
            all_ftp_files_out = pd.concat(all_ftp_files_stk_inp)
            all_ftp_files_out.reset_index(drop=True, inplace=True)
        else:
            all_ftp_files_out = pd.Series([], dtype=str)

        if path_all_ftp_files_save:
            all_ftp_files_out.to_csv(path_all_ftp_files_save)

        return all_ftp_files_out

    table_use = table.copy()

    # Ensure nam_col exists
    if nam_col not in table_use.columns:
        table_use[nam_col] = ""

    prev_host = ""
    prev_dir = ""
    prev_outdir = ""
    ftp_files_list = []
    all_ftp_files_stk = []
    local_files_lis = []
    count_loop = 0
    count_nmax = 50
    ftpobj = None

    for irow, row in table_use.iterrows():
        if row["crawled"]:
            continue

        # Check local files when output directory changes
        if prev_outdir != row["outdir"]:
            local_files_lis = glob.glob(row["outdir"] + "/*")
            prev_outdir = row["outdir"]

        # Check if file already exists locally
        ok_loc, ok_dwl, filename = check_local_file_exists(
            row[rgx_col], row["outdir"], local_files_lis, force, all_files_mode
        )

        if ok_loc:
            table_use.loc[irow, "ok_loc"] = True
            table_use.loc[irow, "ok_dwl"] = False
            table_use.loc[irow, nam_col] = filename
            continue
        elif ok_dwl and force:
            table_use.loc[irow, "ok_loc"] = False
            table_use.loc[irow, "ok_dwl"] = True
            table_use.loc[irow, nam_col] = filename

        count_loop += 1

        # Get or create FTP connection
        ftpobj, prev_host = get_ftp_connection(
            ftpobj,
            row["host"],
            row["protocol"],
            sftp,
            user,
            passwd,
            prev_host,
            count_loop,
            count_nmax,
        )

        # Save intermediate results and reset counter on reconnection
        if count_loop > count_nmax:
            count_loop = 0
            _save_crawled_files(table_use)
            _get_and_save_all_ftp_files(all_ftp_files_stk)

        # Get file list when directory changes
        ftp_result = get_ftp_directory_listing(
            ftpobj, row["dir"], row["host"], prev_dir
        )
        if ftp_result[0] is not None:
            ftp_files_list, ftp_files_urls = ftp_result
            prev_dir = row["dir"]

            if ftp_files_urls:
                all_ftp_files_stk.append(pd.Series(ftp_files_urls))

        # Match files on server using regex pattern
        file_match = match_files_in_directory(
            row[rgx_col], ftp_files_list, all_files_mode
        )

        # Apply exclude patterns
        if exclude_patterns and file_match:
            if all_files_mode and isinstance(file_match, list):
                for negpatt in exclude_patterns:
                    file_match = [
                        f for f in file_match if not re_module.search(negpatt, f)
                    ]
            elif isinstance(file_match, str):
                for negpatt in exclude_patterns:
                    if re_module.search(negpatt, file_match):
                        file_match = None
                        break

        # Update table based on file availability
        # Use a local wrapper to handle both column naming conventions
        if file_match:
            table_use.loc[irow, "ok_dwl"] = True
            if all_files_mode and isinstance(file_match, list):
                table_use.loc[irow, nam_col] = ";".join(file_match)
                log.info(f"{len(file_match)} file(s) found on server :)")
            else:
                table_use.loc[irow, nam_col] = file_match
                log.info(str(file_match) + " found on server :)")
        else:
            table_use.loc[irow, "ok_dwl"] = False
            table_use.loc[irow, nam_col] = ""
            log.warning(f"{row[rgx_col]} not found on server :(")

        table_use.loc[irow, "crawled"] = True

    # Generate URLs for downloadable files
    generate_download_urls(table_use, all_files_mode)

    # Save final results
    _save_crawled_files(table_use)
    all_ftp_files = _get_and_save_all_ftp_files(all_ftp_files_stk)

    # Clean up FTP connection
    if ftpobj:
        ftpobj.close()

    # Collect local file paths
    all_loc_files = collect_local_files(table_use)

    return table_use, all_ftp_files, all_loc_files


#### GRAVEYARD


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
