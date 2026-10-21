#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download utilities for GNSS data and products.

This sub-module of geodezyx.operational contains functions to download
GNSS data and products from distant IGS servers.

It can be imported directly with::

    from geodezyx import operational

The geodezyx toolbox is a software for simple but useful functions for
Geodesy and Geophysics under the GNU LGPL v3 License.

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)

GitHub repository:
    https://github.com/IPGP/geodezyx
"""

import ftplib
import requests
import tqdm
import socket
import htmllistparse
import re

#### Import the logger
import logging
import os
import pathlib
import time

########## BEGIN IMPORT ##########
#### External modules
from ftplib import FTP, FTP_TLS
from multiprocessing.dummy import Pool as ThreadPool

import numpy as np
import pandas as pd

#### geodezyx modules
from geodezyx import conv
from geodezyx import utils

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
    Generate the effective save directory path for orbit files.

    Internal function used for organizing downloaded orbit files.

    Parameters
    ----------
    parent_archive_dir : str
        Parent archive directory path.
    calc_center : str
        Calculation center identifier.
    date : dt.datetime
        Date for organizing the directory structure.
    archtype : str, optional
        Archive directory structure pattern. Default is "year/doy/".
        Supported patterns:
        - "stat": Station directory
        - "stat/year": Station/Year structure
        - "stat/year/doy": Station/Year/DOY structure
        - "year/doy": Year/DOY structure
        - "year/stat": Year/Station structure
        - "week/dow": GPS week/Day of week structure
        - "wkwwww": GFZ CF-ORB wk<wwww> naming convention
        - "/": Save directly in parent folder

    Returns
    -------
    str
        Effective save directory path.
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


def effective_save_dir(parent_archive_dir, site, date, archtype="stat"):
    """
    Generate the effective save directory path for station files.

    Internal function used for organizing downloaded station files.

    Parameters
    ----------
    parent_archive_dir : str
        Parent archive directory path.
    site : str
        Station site code.
    date : dt.datetime
        Date for organizing the directory structure.
    archtype : str, optional
        Archive directory structure pattern. Default is "stat".
        Supported patterns:
        - "stat": Station directory
        - "stat/year": Station/Year structure
        - "stat/year/doy": Station/Year/DOY structure
        - "year/doy": Year/DOY structure
        - "year/stat": Year/Station structure
        - "week/dow": GPS week/Day of week structure
        - "/": Save directly in parent folder

    Returns
    -------
    str
        Effective save directory path.

    Notes
    -----
    If site is None or empty, it will be replaced with "ALL_STATIONS" in the path.
    """
    if archtype == "/":
        return parent_archive_dir

    if len(archtype) > 0 and archtype.startswith("/"):
        log.warning("The archive type starts with /, remove it to avoid error")

    out_save_dir = parent_archive_dir
    fff = archtype.split("/")
    year = str(date.year)
    doy = conv.dt2doy(date)
    _, _ = year, doy  ## simply to remove the unused linter warning...
    week, dow = conv.dt2gpstime(date)

    # If site is None or empty, replace with a placeholder
    if not site:
        site = "ALL_STATIONS"

    stat = site  # alias: archtype tokens may use either 'site' or 'stat'

    for f in fff:
        out_save_dir = os.path.join(out_save_dir, str(eval(f)))
    return out_save_dir


#  _    _ _______ _______ _____    _____                      _                 _
# | |  | |__   __|__   __|  __ \  |  __ \                    | |               | |
# | |__| |  | |     | |  | |__) | | |  | | _____      ___ __ | | ___   __ _  __| |
# |  __  |  | |     | |  |  ___/  | |  | |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |
# | |  | |  | |     | |  | |      | |__| | (_) \ v  v /| | | | | (_) | (_| | (_| |
# |_|  |_|  |_|     |_|  |_|      |_____/ \___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|


#### HTTP classic Download


class TimeoutHTTPAdapter(requests.adapters.HTTPAdapter):
    """
    Custom HTTP adapter that enforces a timeout on all requests.

    This adapter ensures that socket timeouts are properly applied to
    streaming HTTP requests, preventing indefinite hangs during downloads.
    """

    def __init__(self, timeout, *args, **kwargs):
        self.timeout = timeout
        super().__init__(*args, **kwargs)

    def send(self, request, **kwargs):
        kwargs["timeout"] = self.timeout
        return super().send(request, **kwargs)


def download_http(url, output_dir, timeout=120, max_try=4, sleep_time=5, force=False):
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
        This timeout applies to:
        - Initial connection
        - Headers reading
        - Each chunk reading (per iteration)
    max_try : int, optional
        The maximum number of retry attempts in case of failure. Default is 4.
    sleep_time : int, optional
        The sleep time between retry attempts in seconds. Default is 5 seconds.
    force : bool, optional
        If True, forces the download even if the file already exists locally.
        Default is False.

    Returns
    -------
    output_path : str
        The path to the downloaded file, or the URL string if the download
        failed.
    downloaded : bool
        ``True`` if the file was actually fetched from the server during this
        call, ``False`` if it already existed locally or the download failed.

    Notes
    -----
    The timeout is enforced using socket.settimeout() on the underlying
    HTTP connection, which ensures that the entire download operation
    (including streaming) will timeout if no data is received for the
    specified duration.
    """
    log.info("Download file: %s", url)

    # Construct output path and ensure output directory exists
    filename = url.split("/")[-1]
    output_path = os.path.join(output_dir, filename)
    os.makedirs(output_dir, exist_ok=True)

    dwl = False

    # Check if file already exists – skip network entirely
    if os.path.isfile(output_path) and not force:
        log.info("%s already exists locally ;) – skipping.", filename)
        return output_path, dwl

    # Try to get expected file size for the progress bar (best-effort)
    try:
        head_resp = requests.head(url, timeout=timeout)
        file_size = int(head_resp.headers.get("content-length", 0))
    except requests.exceptions.RequestException:
        file_size = 0

    # Download file with progress bar and retry logic
    try_count = 0
    while try_count < max_try:
        try:
            try_count += 1

            # Create a session with custom timeout handling
            session = requests.Session()
            session.mount("http://", TimeoutHTTPAdapter(timeout=timeout))
            session.mount("https://", TimeoutHTTPAdapter(timeout=timeout))

            response = session.get(url, stream=True, timeout=timeout)
            response.raise_for_status()

            with open(output_path, "wb") as f:
                with tqdm.tqdm(
                    total=file_size, unit="B", unit_scale=True, desc=filename
                ) as pbar:
                    # Read chunks with timeout handling
                    for data in response.iter_content(chunk_size=8192):
                        if not data:
                            break
                        f.write(data)
                        pbar.update(len(data))

            session.close()
            dwl = True
            return output_path, dwl

        except (
            socket.timeout,
            requests.exceptions.Timeout,
            requests.exceptions.ConnectionError,
            requests.exceptions.RequestException,
        ) as e:
            log.warning("download failed (%s), try %i/%i", str(e), try_count, max_try)

            if try_count >= max_try:
                log.error("download failed after %i attempts", max_try)
                # Clean up incomplete file
                if os.path.exists(output_path):
                    try:
                        os.remove(output_path)
                    except Exception:
                        pass
                return url, dwl

            time.sleep(sleep_time)

    return output_path, dwl


#  ______ _______ _____    _____                      _                 _
# |  ____|__   __|  __ \  |  __ \                    | |               | |
# | |__     | |  | |__) | | |  | | _____      ___ __ | | ___   __ _  __| |
# |  __|    | |  |  ___/  | |  | |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |
# | |       | |  | |      | |__| | (_) \ v  v /| | | | | (_) | (_| | (_| |
# |_|       |_|  |_|      |_____/ \___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|


#### FTP DOWNLOAD


class MyFTP_TLS(FTP_TLS):
    """
    Custom FTPS client that reuses TLS session.

    This class is a subclass of FTP_TLS from the ftplib module. It creates
    an FTPS client that shares the TLS session to avoid the error:
    ssl.SSLEOFError: EOF occurred in violation of protocol (_ssl.c:2396)

    Notes
    -----
    Source: https://stackoverflow.com/questions/14659154/ftps-with-python-ftplib-session-reuse-required
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


def create_ftp_obj(
    secure_ftp_inp=False,
    host="",
    chdir="",
    parallel_download=1,
    user="anonymous",
    passwd="",
    retry_count=3,
    timeout=120,
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
    timeout : int, optional
        The timeout for FTP connections in seconds. Default is 120 seconds.

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
                current_ftp_obj = ftp_constuctor(host, timeout=timeout)
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
    if len(ftp_obj_list_out) > 0:
        ftp_main = ftp_obj_list_out[0]
        # change the directory of the main ftp obj if we ask for it
        if chdir:
            log.info("Move to: %s", chdir)
            ftp_main.cwd(chdir)
    else:
        errmsg = "No FTP object(s) were created successfully :("
        log.error(errmsg)
        ftp_main = None

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
    Wrapper function for FTP download operations.

    Unpacks the input tuple and passes it to the ftp_downld_mono function.

    Parameters
    ----------
    intup : tuple
        A tuple containing the parameters to be passed to ftp_downld_mono.

    Returns
    -------
    tuple
        The output of the ftp_downld_mono function.
    """
    outtup = ftp_downld_mono(*intup)
    return outtup


def downld_front(
    urls,
    savedirs,
    protocols=None,
    parallel_download=1,
    secure_ftp=False,
    user="anonymous",
    passwd="anonymous@isp.com",
    force=True,
    timeout=120,
):
    """
    Generic download function that handles both FTP/SFTP and HTTP/HTTPS protocols.

    This function automatically detects the protocol from URLs and calls the
    appropriate download function (FTP or HTTP) for each file.

    Parameters
    ----------
    urls : str or list
        The URL or list of URLs of the files to be downloaded.
    savedirs : str or list
        The directory or list of directories where the downloaded files should be saved.
    protocols : str or list, optional
        Protocol for each URL (ftp, sftp, http, https). If None, detected from URL.
    parallel_download : int, optional
        The number of parallel downloads to be performed. Default is 1.
    secure_ftp : bool or list, optional
        If True, uses FTPS for secure file transfer. Can be per-URL list. Default is False.
    user : str, optional
        The username for the FTP server. Default is "anonymous".
    passwd : str, optional
        The password for the FTP server. Default is 'anonymous@isp.com'.
    force : bool, optional
        If True, forces the download even if the file already exists. Default is True.
    timeout : int, optional
        The timeout for connections in seconds. Default is 120 seconds.

    Returns
    -------
    out_tup_lis : list of tuple
        List of tuples containing the local path of the downloaded file and
        a boolean indicating whether the download was successful.
        e.g. [(local_path1, True), (local_path2, False), ...]

    Notes
    -----
    This function dispatches to ftp_downld_front() for FTP/SFTP URLs and
    http_downld_front() for HTTP/HTTPS URLs.
    """
    import re as re_module

    # Check if urls and savedirs are iterable, if not convert them to list
    urllist = urls if utils.is_iterable(urls) else [urls]
    urllist = list(urllist)
    savedirlist = savedirs if utils.is_iterable(savedirs) else [savedirs] * len(urllist)

    if isinstance(secure_ftp, (list, tuple)):
        secure_ftp_list = list(secure_ftp)
    else:
        secure_ftp_list = [secure_ftp] * len(urllist)

    # Check if the length of urllist and savedirlist are the same
    if len(urllist) != len(savedirlist):
        log.error(
            "URL & out dirs lists do not have the same length: %s %s",
            len(urllist),
            len(savedirlist),
        )
        return []

    # Detect protocol from URLs if not provided
    if protocols is None:
        protocols = []
        for url in urllist:
            if url.startswith(("http://", "https://")):
                protocols.append("http" if url.startswith("http://") else "https")
            elif url.startswith(("ftp://", "sftp://")):
                protocols.append("ftp" if url.startswith("ftp://") else "sftp")
            else:
                # Fallback: try to parse from path structure
                match = re_module.match(r"^([a-z]+)://", url)
                if match:
                    protocols.append(match.group(1))
                else:
                    protocols.append("ftp")  # Default to FTP

    # Separate URLs by protocol
    ftp_indices = []
    http_indices = []

    for idx, proto in enumerate(protocols):
        if proto in ("ftp", "sftp"):
            ftp_indices.append(idx)
        elif proto in ("http", "https"):
            http_indices.append(idx)

    out_tup_lis = []

    # Download FTP/SFTP files
    if ftp_indices:
        ftp_urls = [urllist[i] for i in ftp_indices]
        ftp_savedirs = [savedirlist[i] for i in ftp_indices]
        ftp_secure = [secure_ftp_list[i] for i in ftp_indices]

        log.info("Downloading %d FTP/SFTP file(s)", len(ftp_urls))
        ftp_results = ftp_downld_front(
            ftp_urls,
            ftp_savedirs,
            parallel_download=parallel_download,
            secure_ftp=ftp_secure,
            user=user,
            passwd=passwd,
            force=force,
            timeout=timeout,
        )
        out_tup_lis.extend(ftp_results)

    # Download HTTP/HTTPS files
    if http_indices:
        http_urls = [urllist[i] for i in http_indices]
        http_savedirs = [savedirlist[i] for i in http_indices]

        log.info("Downloading %d HTTP/HTTPS file(s)", len(http_urls))
        http_results = http_downld_front(
            http_urls,
            http_savedirs,
            parallel_download=parallel_download,
            force=force,
            timeout=timeout,
        )
        out_tup_lis.extend(http_results)

    return out_tup_lis


def http_downld_wrap(intup):
    """
    Wrapper function for HTTP download operations.

    Parameters
    ----------
    intup : tuple
        A tuple containing (url, output_dir, force, timeout).

    Returns
    -------
    tuple
        Result of the download_http function.
    """
    url, output_dir, force, timeout = intup
    return download_http(url, output_dir, timeout=timeout, force=force)


def http_downld_front(
    urls,
    savedirs,
    parallel_download=1,
    force=False,
    timeout=120,
):
    """
    Download files from HTTP/HTTPS servers in parallel.

    Parameters
    ----------
    urls : list
        List of URLs of the files to be downloaded.
    savedirs : list
        List of directories where the downloaded files should be saved.
    parallel_download : int, optional
        The number of parallel downloads to be performed. Default is 1.
    force : bool, optional
        If True, forces the download even if the file already exists. Default is False.
    timeout : int, optional
        The timeout for HTTP connections in seconds. Default is 120 seconds.

    Returns
    -------
    out_tup_lis : list of tuple
        List of tuples containing the local path of the downloaded file and
        a boolean indicating whether the download was successful.
    """
    # Ensure inputs are lists
    urllist = list(urls) if utils.is_iterable(urls) else [urls]
    savedirlist = (
        list(savedirs)
        if utils.is_iterable(savedirs)
        else [savedirs] * len(urllist)
    )

    if len(urllist) != len(savedirlist):
        log.error(
            "URL & out dirs lists do not have the same length: %s %s",
            len(urllist),
            len(savedirlist),
        )
        return []

    # Create list of arguments for parallel downloads
    force_lis = [force] * len(urllist)
    timeout_lis = [timeout] * len(urllist)

    # Create a ThreadPool for parallel downloads
    pool = ThreadPool(parallel_download)

    # Start the parallel downloads
    out_tup_lis = pool.map(
        http_downld_wrap, list(zip(urllist, savedirlist, force_lis, timeout_lis))
    )
    pool.close()
    pool.join()

    return out_tup_lis


def ftp_downld_front(
    urls,
    savedirs,
    parallel_download=1,
    secure_ftp=False,
    user="anonymous",
    passwd="anonymous@isp.com",
    force=True,
    timeout=120,
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
    timeout : int, optional
        The timeout for FTP connections in seconds. Default is 120 seconds.

    Returns
    -------
    out_tup_lis : list of tuple
        List of tuples containing the local path of the downloaded file and
        a boolean indicating whether the download was successful.
        e.g. [(local_path1, True), (local_path2, False), ...]

    Notes
    -----
    This function uses the ThreadPool for parallel downloads.
    """

    # Check if urls and savedirs are iterable, if not convert them to list
    urllist = urls if utils.is_iterable(urls) else [urls]
    urllist = list(urllist)
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
    urlpathobj = pd.Series([pathlib.Path(e) for e in urllist])
    host_use = urlpathobj.apply(lambda p: p.parts[1]).unique()[0]

    # Create the FTP object
    ftpobj_main, ftpobj_lis = create_ftp_obj(
        secure_ftp_inp=secure_ftp_use,
        host=host_use,
        parallel_download=parallel_download,
        user=user,
        passwd=passwd,
        timeout=timeout,
    )

    if not ftpobj_main:
        warnmsg = (
            "No FTP object could be created successfully. FTP download is aborted."
        )
        log.warning(warnmsg)
        return []

    # Create a list of FTP objects for parallel downloads
    ftpobj_mp_lis = ftpobj_lis * int(np.ceil(len(urllist) / parallel_download))
    force_lis = [force] * len(urllist)

    # Check if there are less FTP objects than URLs for parallel download
    if len(ftpobj_mp_lis) < len(urllist):
        warnmsg = "less FTP objects than URL for parallel download, contact the main developper"
        log.warning(warnmsg)

    # Create a ThreadPool for parallel downloads
    pool = ThreadPool(parallel_download)

    # Start the parallel downloads
    out_tup_lis = pool.map(
        ftp_downld_wrap, list(zip(ftpobj_mp_lis, urllist, savedirlist, force_lis))
    )

    return out_tup_lis


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


def check_loc_file(
    filergx, outdir, local_files_cache, force=False, all_files_mode=False
) -> tuple[bool, bool, str]:
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
    ok_loc : bool
        True if file exists locally and should not be re-downloaded.
    ok_dwl : bool
        True if file should be downloaded (not local or force=True).
    filename : str
        Filename if found locally (semicolon-separated if
        all_files_mode=True), empty string otherwise.
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


def renew_ftp_obj(
    ftpobj,
    host,
    protocol,
    sftp,
    user,
    passwd,
    prev_host,
    count_loop,
    count_nmax,
    timeout=120,
) -> tuple[FTP | None, str]:
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
    timeout : int, optional
        The timeout for FTP connections in seconds. Default is 120 seconds.

    Returns
    -------
    ftpobj : FTP or None
        FTP connection object.
    host_out : str
        The host name for connection tracking.
    """
    if host != prev_host or count_loop > count_nmax or count_loop == 1:
        if ftpobj:
            ftpobj.close()

        if sftp == "auto":
            sftp_use = protocol == "sftp"
        else:
            sftp_use = bool(sftp)

        ftpobj, _ = create_ftp_obj(
            secure_ftp_inp=sftp_use,
            host=host,
            user=user,
            passwd=passwd,
            timeout=timeout,
        )

        if not ftpobj:
            warnmsg = "No FTP object could be created successfully."
            log.warning(warnmsg)

        host_out = host
        return ftpobj, host_out

    host_out = prev_host
    return ftpobj, host_out


def list_ftp_dir_core(ftpobj):
    """
    Lists the files in the current directory of the FTP object.

    Parameters
    ----------
    ftpobj : FTP object
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
        files = ftpobj.nlst()
    except ftplib.error_perm as resp:
        if str(resp) == "550 No files found":
            log.warning("No files in this directory" + ftpobj.pwd())
        else:
            raise
    return files


def list_ftp_dir(
    ftpobj, dir_inp, host, prev_dir
) -> tuple[list[str] | None, list[str] | None]:
    """
    Get FTP directory file listing.

    Changes to the specified directory and retrieves the list of files.
    Returns empty list if directory change fails.

    Parameters
    ----------
    ftpobj : FTP
        FTP connection object.
    dir_inp : str
        Remote directory path.
    host : str
        FTP server hostname (for logging/URL generation).
    prev_dir : str
        Previous directory path.

    Returns
    -------
    ftp_files_list : list of str
        Filenames in the remote directory. None if directory is unchanged.
    ftp_files_urls : list of str
        Complete FTP URLs for all files in the directory. None if directory
        is unchanged. Empty list if the directory change failed.
    """
    if prev_dir != dir_inp:
        log.info("chdir " + dir_inp)
        ftpobj.cwd("/")

        try:
            ftpobj.cwd(dir_inp)
            ftp_files_list = list_ftp_dir_core(ftpobj)
            ftp_files_urls = [
                f"ftp://{host}/{dir_inp.lstrip('/')}/{f}" for f in ftp_files_list
            ]
            return ftp_files_list, ftp_files_urls

        except Exception as e:
            log.warning("unable to chdir to %s, exception %s", dir_inp, e)
            return [], []

    return None, None


def match_files_in_dir(filergx, ftp_files_list, all_files_mode=False):
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

    Returns
    -------
    None
        Updates the input table in-place with crawl results.
    """
    if file_match:
        table.loc[irow, "ok_dwl"] = True
        if all_files_mode:
            table.loc[irow, "filename"] = ";".join(file_match)
            log.info(f"{len(file_match)} file(s) found on server :)")
        else:
            table.loc[irow, "filename"] = file_match
            log.info(file_match + " found on server :)")
    else:
        table.loc[irow, "ok_dwl"] = False
        table.loc[irow, "filename"] = ""
        log.warning(f"{table.loc[irow, 'regex']} not found on server :(")

    table.loc[irow, "crawled"] = True


def generate_download_urls(table, all_files_mode=False):
    """
    Generate complete URLs for files to be downloaded using the correct protocol.

    Parameters
    ----------
    table : pd.DataFrame
        Table with crawl results. Must have columns: ok_dwl, ok_loc, host,
        dir, filename, protocol.
    all_files_mode : bool, optional
        If True, handle semicolon-separated filenames. Default is False.

    Returns
    -------
    None
        Updates the 'url_true' column in the input table in-place.
    """
    nam_col = "filename"

    fil_ok_dwl = table["ok_dwl"] & ~table["ok_loc"]

    if all_files_mode:

        def make_urls(row):
            filenames = row[nam_col].split(";")
            # Get protocol, default to ftp if not specified
            protocol = row.get("protocol", "ftp").lower()
            # Normalize dir: RINEX dir has no leading '/' (built from pathlib.Path parts),
            # while products dir has a leading '/' (built from os.path.join with absolute
            # basedir).  Strip any leading '/' and always add exactly one.
            dir_clean = row["dir"].lstrip("/")
            urls = [f"{protocol}://{row['host']}/{dir_clean}/{fname}" for fname in filenames]
            return ";".join(urls)

        table.loc[fil_ok_dwl, "url_true"] = table.loc[fil_ok_dwl].apply(
            make_urls, axis=1
        )
    else:
        # Same normalization for the single-file case
        def make_url(row):
            protocol = row.get("protocol", "ftp").lower()
            return f"{protocol}://{row['host']}/{row['dir'].lstrip('/')}/{row[nam_col]}"

        table.loc[fil_ok_dwl, "url_true"] = table.loc[fil_ok_dwl].apply(
            make_url, axis=1
        )


def collect_local_files(table):
    """
    Collect local file paths for files that exist locally.

    Parameters
    ----------
    table : pd.DataFrame
        Table with crawl results containing 'ok_loc', 'outdir' columns,
        and 'filename' column.

    Returns
    -------
    pd.Series
        Series of local file paths, or empty Series if no local files.
    """
    nam_col = "filename"

    ok_loc = table["ok_loc"]
    if ok_loc.sum() > 0:
        return pd.Series(
            table.loc[ok_loc].apply(
                lambda e: os.path.join(e["outdir"], e[nam_col]), axis=1
            )
        )
    else:
        return pd.Series([])


def list_http_files(url_base, date_start=None, date_end=None):
    """
    Generic function to list files available on an HTTP server with year/day structure.

    This function works with any HTTP server that organizes files in the following
    directory structure: {url_base}/{year}/{doy}/ where doy is the day of year.

    If url_base already ends with /year/doy/, it will list files directly from that directory.

    Parameters
    ----------
    url_base : str
        The base URL of the repository (e.g., 'http://example.com/rinex/').
        Can end with a trailing slash or include year/doy structure.
    date_start : str or datetime.date, optional
        The starting date for filtering files. Can be in formats: 'YYYY-MM-DD',
        'YYYYMMDD', day-of-year 'YYYY-DDD', or as a datetime object.
        If None, no lower bound is applied. Default is None.
    date_end : str or datetime.date, optional
        The ending date for filtering files. Can be in formats: 'YYYY-MM-DD',
        'YYYYMMDD', day-of-year 'YYYY-DDD', or as a datetime object.
        If None, no upper bound is applied. Default is None.

    Returns
    -------
    list of str
        A list of URLs pointing to the files that match the specified date criteria.

    Notes
    -----
    - The function fetches directory listings from the server.
    - Dates are converted to day-of-year (doy) for fetching from the server.
    - Logs the number of files found for each day.
    - Works with any HTTP server with year/day-of-year directory structure.

    Examples
    --------
    >>> files = list_http_files(
    ...     "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/",
    ...     date_start="2020-01-01",
    ...     date_end="2020-01-31"
    ... )
    >>> len(files) > 0
    True
    """
    # Parse date inputs using geodezyx's date converter
    if date_start is not None:
        date_start = conv.date_pattern2dt(date_start)
    if date_end is not None:
        date_end = conv.date_pattern2dt(date_end)

    rnx_path_lst = []

    # Detect if url_base already contains year/doy structure using regex
    # Check if last two path components match year/doy pattern (YYYY/DDD)
    url_clean = url_base.rstrip("/")
    has_year_doy = bool(re.search(r'/\d{4}/\d{3}/?$', url_clean))

    # Build list of URLs to fetch and their date labels for logging
    urls_to_fetch = []

    # Case 1: URL already contains year/doy structure - use it directly
    if has_year_doy:
        url_day = url_base if url_base.endswith("/") else url_base + "/"
        urls_to_fetch.append((url_day, None))
    # Case 2: Date range provided - iterate over dates
    elif date_start is not None and date_end is not None:
        for current_date in conv.dt_range(date_start, date_end):
            year = current_date.year
            doy = current_date.timetuple().tm_yday
            url_day = url_base.rstrip("/") + "/" + f"{year}/{doy:03d}/"
            urls_to_fetch.append((url_day, f"{year}-{doy:03d}"))
    # Case 3: Fallback - fetch year/day structure from server listing
    else:
        log.warning(
            "No date range specified. Will attempt to list all available dates from server."
        )
        try:
            cwd, listing = htmllistparse.fetch_listing(url_base, timeout=30)
            for item in listing:
                url_year = url_base.rstrip("/") + "/" + item.name
                year = item.name.strip("/")
                if not re.match(r"^\d{4}/?$", item.name):
                    continue
                year = int(year)

                cwd, listing = htmllistparse.fetch_listing(url_year, timeout=30)
                for subitem in listing:
                    url_day = url_year.rstrip("/") + "/" + subitem.name
                    day = subitem.name.strip("/")
                    urls_to_fetch.append((url_day, f"{year}-{day}"))
        except Exception as e:
            log.warning(f"Failed to fetch listing from {url_base}: {e}")

    # Common processing: fetch and collect files for all URLs
    for url_day, date_label in urls_to_fetch:
        try:
            cwd, listing = htmllistparse.fetch_listing(url_day, timeout=30)
            for file in listing:
                rnx_path = url_day.rstrip("/") + "/" + file.name
                rnx_path_lst.append(rnx_path)
            # Log with appropriate label
            if date_label:
                log.info(f"{len(listing)} files found for day {date_label}")
            else:
                log.info(f"{len(listing)} files found in {url_day}")
        except Exception as e:
            log.warning(f"Failed to fetch listing from {url_day}: {e}")

    return rnx_path_lst


def crawl_files(
    table,
    sftp="auto",
    user=None,
    passwd=None,
    path_ftp_crawled_files_save=None,
    path_all_ftp_files_save=None,
    force=False,
    all_files_mode=False,
    exclude_patterns=None,
    timeout=120,
) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    """
    Crawl FTP/SFTP/HTTP/HTTPS servers to find available files and update download table.

    This is a generic function used by both RINEX and products downloaders.
    It checks for existing local files, connects to remote servers, lists remote
    files, and updates the table with availability status and actual file URLs.

    Supports both FTP-based protocols (ftp, sftp) and HTTP-based protocols (http, https).

    Parameters
    ----------
    table : pd.DataFrame
        Input table containing download metadata with columns:
        - 'host': Server hostname
        - 'dir': Remote directory path
        - 'outdir': Local output directory
        - 'regex': Filename regex pattern
        - 'protocol': String indicating protocol ("ftp", "sftp", "http", or "https")
        - 'url_theo': Full theoretical URL (for HTTP)
        - 'crawled': Boolean indicating if already crawled
    sftp : str or bool, optional
        SFTP mode setting. Default is 'auto'.
    user : str, optional
        Server username. Default is None (anonymous).
    passwd : str, optional
        Server password. Default is None (anonymous).
    path_ftp_crawled_files_save : str, optional
        Path to save the crawled files table as CSV.
    path_all_ftp_files_save : str, optional
        Path to save all discovered files as CSV.
    force : bool, optional
        Force re-download even if files exist locally. Default is False.
    all_files_mode : bool, optional
        If True, download ALL files matching the regex. Default is False.
    exclude_patterns : list of str, optional
        List of regex patterns to exclude from matched files. Default is None.
    timeout : int, optional
        Connection timeout in seconds. Default is 120 seconds.

    Returns
    -------
    table_use : pd.DataFrame
        Updated table with crawl results.
    all_files : pd.Series
        All files discovered on servers.
    all_loc_files : pd.Series
        Local file paths for files that already exist.
    """
    import glob
    import re as re_module

    rgx_col = "regex"
    nam_col = "filename"

    def _save_crawled_files(table_inp):
        if path_ftp_crawled_files_save:
            table_inp.to_csv(path_ftp_crawled_files_save)

    def _save_all_files(all_files_stk_inp):
        if all_files_stk_inp:
            all_files_out = pd.concat(all_files_stk_inp)
            all_files_out.reset_index(drop=True, inplace=True)
        else:
            all_files_out = pd.Series([], dtype=str)

        if path_all_ftp_files_save:
            all_files_out.to_csv(path_all_ftp_files_save, header=False)

        return all_files_out

    table_use = table.copy()

    # Ensure nam_col exists
    if nam_col not in table_use.columns:
        table_use[nam_col] = ""

    prev_host = ""
    prev_dir = ""
    prev_outdir = ""
    ftp_files_list = []
    all_files_stk = []
    local_files_lis = []
    count_loop = 0
    count_nmax = 50
    ftpobj = None

    for irow, row in table_use.iterrows():
        if row["crawled"]:
            continue

        protocol = row.get("protocol", "ftp").lower()
        is_http = protocol in ("http", "https")

        # Check local files when output directory changes
        if prev_outdir != row["outdir"]:
            local_files_lis = glob.glob(row["outdir"] + "/*")
            prev_outdir = row["outdir"]

        # Check if file already exists locally
        outchk = check_loc_file(
            row[rgx_col], row["outdir"], local_files_lis, force, all_files_mode
        )
        ok_loc, ok_dwl, filename = outchk

        if ok_loc:
            table_use.loc[irow, "ok_loc"] = True
            table_use.loc[irow, "ok_dwl"] = False
            table_use.loc[irow, nam_col] = filename
            table_use.loc[irow, "crawled"] = True
            continue
        elif ok_dwl and force:
            table_use.loc[irow, "ok_loc"] = False
            table_use.loc[irow, "ok_dwl"] = True
            table_use.loc[irow, nam_col] = filename

        count_loop += 1

        # Handle HTTP/HTTPS protocol
        if is_http:
            url_theo = row.get("url_theo", "")
            try:
                # Build URL from structured fields: protocol, host, dir
                protocol = row.get("protocol", "http").lower()
                host = row.get("host", "")
                dir_path = row.get("dir", "").lstrip("/")
                url_base = f"{protocol}://{host}/{dir_path}/" if dir_path else f"{protocol}://{host}/"

                # Parse date from row if available
                date_val = row.get("date", None)

                # Use list_http_files to fetch available files
                http_files_list = list_http_files(url_base, date_start=date_val, date_end=date_val)

                if http_files_list:
                    all_files_stk.append(pd.Series(http_files_list))
                    # Try to match with regex
                    http_filenames = [str(os.path.basename(url)) for url in http_files_list]
                    file_match = match_files_in_dir(row[rgx_col], http_filenames, all_files_mode)

                    if file_match:
                        # Find the full URL(s) for the matched file(s)
                        if all_files_mode and isinstance(file_match, list):
                            matching_urls = [
                                url for url in http_files_list
                                if os.path.basename(url) in file_match
                            ]
                            if matching_urls:
                                table_use.loc[irow, "url_true"] = ";".join(matching_urls)
                                table_use.loc[irow, "ok_dwl"] = True
                                table_use.loc[irow, nam_col] = ";".join(file_match)
                                log.info(f"{len(file_match)} file(s) found on server :)")
                            else:
                                table_use.loc[irow, "ok_dwl"] = False
                        else:
                            matching_url = [
                                url for url in http_files_list
                                if url.endswith(file_match)
                            ]
                            if matching_url:
                                table_use.loc[irow, "url_true"] = matching_url[0]
                                table_use.loc[irow, "ok_dwl"] = True
                                table_use.loc[irow, nam_col] = file_match
                                log.info(f"{file_match} found on server :)")
                            else:
                                table_use.loc[irow, "ok_dwl"] = False
                    else:
                        table_use.loc[irow, "ok_dwl"] = False
                        log.warning(f"{row[rgx_col]} not found on server :(")
                else:
                    table_use.loc[irow, "ok_dwl"] = False
                    log.warning(f"No files found at {url_theo}")

            except Exception as e:
                table_use.loc[irow, "ok_dwl"] = False
                log.warning(f"Failed to crawl {url_theo}: {e}")

        # Handle FTP/SFTP protocol
        else:
            # Get or create FTP connection
            ftpobj, prev_host = renew_ftp_obj(
                ftpobj,
                row["host"],
                protocol,
                sftp,
                user,
                passwd,
                prev_host,
                count_loop,
                count_nmax,
                timeout=timeout,
            )

            # Save intermediate results and reset counter on reconnection
            if count_loop > count_nmax:
                count_loop = 0
                _save_crawled_files(table_use)
                _save_all_files(all_files_stk)

            # Get file list when directory changes
            if not ftpobj:
                log.warning("FTP object is None, skipping row %d", irow)
                table_use.loc[irow, "crawled"] = True
                continue
            else:
                ftp_result = list_ftp_dir(ftpobj, row["dir"], row["host"], prev_dir)

            file_match = None
            if ftp_result[0] is not None:
                ftp_files_list, ftp_files_urls = ftp_result
                prev_dir = row["dir"]

                if ftp_files_urls:
                    all_files_stk.append(pd.Series(ftp_files_urls))
                # Match files on server using regex pattern
                if ftp_files_list:
                    file_match = match_files_in_dir(row[rgx_col], ftp_files_list, all_files_mode)

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
    all_files = _save_all_files(all_files_stk)

    # Clean up FTP connection
    if ftpobj:
        ftpobj.close()

    # Collect local file paths
    all_loc_files = collect_local_files(table_use)

    return table_use, all_files, all_loc_files


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
    timeout=120,
) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    """
    Crawl FTP servers to find available files and update download table.

    This is a wrapper around the generic crawl_files function.
    It checks for existing local files, connects to FTP servers, lists remote
    files, and updates the table with availability status and actual file URLs.

    Now also supports HTTP/HTTPS servers through crawl_files.

    Parameters
    ----------
    table : pd.DataFrame
        Input table containing download metadata with columns:
        - 'host': FTP server hostname
        - 'dir': Remote directory path
        - 'outdir': Local output directory
        - 'regex': Filename regex pattern
        - 'protocol': String indicating protocol ("ftp", "sftp", "http", or "https")
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
        Path to save all discovered files as CSV.
    force : bool, optional
        Force re-download even if files exist locally. Default is False.
    all_files_mode : bool, optional
        If True, download ALL files matching the regex.
        Default is False.
    exclude_patterns : list of str, optional
        List of regex patterns to exclude from matched files.
        Default is None.
    timeout : int, optional
        The timeout for connections in seconds. Default is 120 seconds.

    Returns
    -------
    table_use : pd.DataFrame
        Updated table with crawl results.
    all_files : pd.Series
        All files discovered on servers.
    all_loc_files : pd.Series
        Local file paths for files that already exist.
    """
    # Call the generic crawl_files function which handles both FTP and HTTP
    return crawl_files(
        table,
        sftp=sftp,
        user=user,
        passwd=passwd,
        path_ftp_crawled_files_save=path_ftp_crawled_files_save,
        path_all_ftp_files_save=path_all_ftp_files_save,
        force=force,
        all_files_mode=all_files_mode,
        exclude_patterns=exclude_patterns,
        timeout=timeout,
    )
