#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.utils contains functions for Shell-like
 operations in Python.

it can be imported directly with:
from geodezyx import utils

The geodezyx toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/IPGP/geodezyx
"""


########## BEGIN IMPORT ##########
#### External modules
import fnmatch
import glob
import gzip
#### Import the logger
import logging
import os
import re
import shutil
import subprocess
import zipfile

#### geodeZYX modules
from geodezyx import utils

log = logging.getLogger('geodezyx')

##########  END IMPORT  ##########


def gzip_compress(inp_path, out_dir=None, out_fname=None, rm_inp=False):
    """
    Compress a file using gzip.

    Parameters
    ----------
    inp_path : str
        Path to the input file to compress.
    out_dir : str, optional
        Output directory. Default is None (same as input file).
    out_fname : str, optional
        Output filename. Default is None (input filename + ".gz").
    rm_inp : bool, optional
        If True, remove the input file after compression. Default is False.

    Returns
    -------
    str
        Path to the compressed output file.
    """
    if not out_dir:
        out_dir = os.path.dirname(inp_path)
    if not out_fname:
        out_fname = os.path.basename(inp_path) + ".gz"
    out_path = os.path.join(out_dir, out_fname)
    with open(inp_path, 'rb') as f_in:
        with gzip.open(out_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    if rm_inp and os.path.isfile(inp_path):
        os.remove(inp_path)

    return out_path

def extract_zip(zip_path: str, extract_dir: str) -> list:
    """
    Extract a zip archive and return the list of extracted file paths.

    Parameters
    ----------
    zip_path : str
        Path to the ``.zip`` file.
    extract_dir : str
        Directory into which the contents will be extracted.

    Returns
    -------
    list of str
        Absolute paths of all extracted files.
    """
    os.makedirs(extract_dir, exist_ok=True)
    extracted = []
    with zipfile.ZipFile(zip_path, "r") as zf:
        for member in zf.namelist():
            zf.extract(member, extract_dir)
            extracted.append(os.path.join(extract_dir, member))
    log.info("Extracted %d file(s) from %s", len(extracted), os.path.basename(zip_path))
    return extracted


def uncompress(pathin,dirout = '', opts='-f'):
    """
    Uncompress a file using the uncompress command.

    .. deprecated::
        Use :func:`geodezyx.files_rw.unzip_gz_z` instead.

    Parameters
    ----------
    pathin : str
        Path to the file to uncompress.
    dirout : str, optional
        Output directory. Default is '' (current directory).
    opts : str, optional
        Options for the uncompress command. Default is '-f'.

    Returns
    -------
    str or None
        Path to the uncompressed file, or None if input file does not exist.
    """
    log.warning("function discontinued, use files_rw.unzip_gz_z() instead")
    if not os.path.isfile(pathin):
        log.error('uncompress : %s doesnt exist !!!', pathin)
        return None
    komand = 'uncompress ' + opts + ' ' + pathin
    subprocess.call([komand], shell=True)
    pathout_temp = '.'.join(pathin.split('.')[:-1])

    if dirout == '':
        pathout = pathout_temp
    else:
        pathout = os.path.join(dirout,os.path.basename(pathout_temp))
        shutil.move(pathout_temp,pathout)
    return pathout