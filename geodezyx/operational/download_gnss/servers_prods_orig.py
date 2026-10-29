#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Server definitions and dispatcher for GNSS product downloads.

Contains ``_server_select_products``, the central dispatcher used
exclusively by :mod:`download_prods`.
"""

import logging

log = logging.getLogger("geodezyx")


def _server_select_prods(archive_center, mgex=False, repro=0):
    """
    Resolve archive center name to FTP host, base directory, and protocol.

    Parameters
    ----------
    archive_center : str
        Name of the IGS archive/data center.
    mgex : bool, optional
        Get MGEX products. Default is False.
    repro : int, optional
        Reprocessing campaign number. Default is 0 (operational).

    Returns
    -------
    host : str
        FTP server hostname
    basedir : str
        Base directory on the server
    protocol : str
        "ftp" or "sftp"
    secure_ftp : str
        whether SFTP is needed
    """
    mgex_str = "mgex/" if mgex else ""
    protocol = "ftp"
    secure_ftp = False

    if archive_center == "cddis":
        host = "gdc.cddis.eosdis.nasa.gov"
        basedir = "/pub/gps/products/" + mgex_str
        protocol = "sftp"
        secure_ftp = True

    elif archive_center == "cddis_glonass":
        host = "cddis.gsfc.nasa.gov"
        basedir = "/pub/glonass/products/" + mgex_str

    elif archive_center == "esa":
        host = "gssc.esa.int"
        basedir = "/gnss/products/" + mgex_str

    elif archive_center == "ign":
        host = "igs.ign.fr"
        basedir = "/pub/igs/products/" + mgex_str

    elif archive_center == "ign_iono":
        host = "igs-rf.ign.fr"
        basedir = "/pub/"

    elif archive_center == "ensg":
        host = "igs.ensg.ign.fr"
        basedir = "/pub/igs/products/" + mgex_str

    elif archive_center == "whu":
        host = "igs.gnsswhu.cn"
        basedir = "/pub/gps/products/" + mgex_str

    elif archive_center == "ign_rf":
        host = "igs-rf.ign.fr"
        basedir = "/pub/" + mgex_str

    elif archive_center == "ensg_rf":
        host = "igs-rf.ensg.ign.fr"
        basedir = "/pub/" + mgex_str

    else:
        log.error("Unknown archive center: %s", archive_center)
        return None, None, None, None

    return host, basedir, protocol, secure_ftp
