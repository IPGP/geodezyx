#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27/08/2025 22:34:53

@author: psakic
"""

import datetime as dt
import os
import re
from io import BytesIO
import numpy as np
import pandas as pd
from geodezyx import files_rw

# Shared column names for the navigation data DataFrame
_NAV_COLS = [
    "SVclockBias",
    "SVclockDrift",
    "SVclockDriftRate",
    "IODE",
    "Crs",
    "DeltaN",
    "M0",
    "Cuc",
    "Eccentricity",
    "Cus",
    "sqrtA",
    "TimeEph",
    "Cic",
    "OMEGA",
    "CIS",
    "Io",
    "Crc",
    "omega",
    "OMEGA DOT",
    "IDOT",
    "CodesL2",
    "GPSWeek",
    "L2Pflag",
    "SVacc",
    "SVhealth",
    "TGD",
    "IODC",
    "TransTime",
    "FitIntvl",
    "spare1",
    "spare2",
]


def _read_rinex_nav_v2(f, fn):
    """
    Parse the body of a RINEX 2 navigation file (file object already positioned
    after the END OF HEADER line).

    Parameters
    ----------
    f : file object
        Open file positioned just after the END OF HEADER line.
    fn : str
        File path (used to determine satellite system from the extension).

    Returns
    -------
    nav : pd.DataFrame
    """
    startcol = 3   # column where numerical data starts on continuation lines
    nfloat = 19    # fixed-width field width
    nline = 7      # continuation lines per record (excluding header line)

    sv = []
    epoch = []
    raws = ""

    while True:
        headln = f.readline()
        if not headln:
            break

        sv.append(headln[:2])
        year = int(headln[2:5])
        if 80 <= year <= 99:
            year += 1900
        elif year < 80:  # good till year 2180
            year += 2000

        epoch.append(
            dt.datetime(
                year=year,
                month=int(headln[5:8]),
                day=int(headln[8:11]),
                hour=int(headln[11:14]),
                minute=int(headln[14:17]),
                second=int(headln[17:20]),
                microsecond=int(headln[21]) * 100000,
            )
        )

        raw = (
            headln[22:].rstrip()
            + "".join(f.readline()[startcol:].rstrip() for _ in range(nline - 1))
            + f.readline()[startcol:79].rstrip()
        )
        raws += raw + "\n"

    raws = raws.replace("D", "E")
    strio = BytesIO(raws.encode())
    darr = np.genfromtxt(strio, delimiter=nfloat)

    nav = pd.DataFrame(darr, epoch, _NAV_COLS)

    rinexnav_type = os.path.basename(fn)[-1]
    if rinexnav_type == "n":
        nav["sys"] = pd.Series(["G"] * len(nav.index), index=nav.index)
    else:
        nav["sys"] = pd.Series(
            [rinexnav_type.upper()] * len(nav.index), index=nav.index
        )
    nav["svni"] = pd.Series(np.array(sv), index=nav.index)

    return nav


def _read_rinex_nav_v3(f):
    """
    Parse the body of a RINEX 3 mixed navigation file (file object already
    positioned after the END OF HEADER line).

    RINEX 3 record lengths vary by constellation:
      - R (GLONASS) and S (SBAS): 3 data lines  -> 15 values total
      - G, C, E, J, I and others: 7 data lines  -> 31 values total

    Shorter records are padded with NaN to reach 31 columns.

    Parameters
    ----------
    f : file object
        Open file positioned just after the END OF HEADER line.

    Returns
    -------
    nav : pd.DataFrame
    """
    # Number of data lines per constellation (default 7)
    _NLINE = {"R": 3, "S": 3}
    _NLINE_DEFAULT = 7
    _NCOLS = 31  # max columns (fixed-width, 19 chars each)

    def _parse_fw(line, col_start):
        """Extract 19-char fixed-width float fields starting at col_start."""
        part = line[col_start:].rstrip("\n").rstrip("\r")
        return [
            part[i : i + 19].replace("D", "E").strip()
            for i in range(0, len(part), 19)
            if part[i : i + 19].strip()
        ]

    sv = []
    epoch = []
    rows = []

    while True:
        headln = f.readline()
        if not headln:
            break

        fields = headln[:23].split(" ")
        sv.append(fields[0])

        epoch.append(
            dt.datetime(
                year=int(fields[1]),
                month=int(fields[2]),
                day=int(fields[3]),
                hour=int(fields[4]),
                minute=int(fields[5]),
                second=int(fields[6]),
            )
        )

        # Determine number of data lines for this constellation
        nline_rec = _NLINE.get(fields[0][0], _NLINE_DEFAULT)

        # Fixed-width parsing: header values at col 23, data line values at col 4
        str_vals = _parse_fw(headln, 23)
        for _ in range(nline_rec):
            str_vals += _parse_fw(f.readline(), 4)

        vals = np.array([float(v) if v else np.nan for v in str_vals], dtype=float)
        # Pad with NaN to reach _NCOLS columns
        if len(vals) < _NCOLS:
            vals = np.concatenate([vals, np.full(_NCOLS - len(vals), np.nan)])
        rows.append(vals)

    darr = np.vstack(rows)
    nav = pd.DataFrame(darr, epoch, _NAV_COLS)

    const = [e[0] for e in sv]
    svni = [int(e[1:]) for e in sv]
    nav["sys"] = pd.Series(np.array(const), index=nav.index)
    nav["svni"] = pd.Series(np.array(svni), index=nav.index)

    return nav


def read_rinex_nav_v3_header(f):
    """
    Parse IONOSPHERIC CORR and TIME SYSTEM CORR records from a RINEX 3
    navigation file header.

    Compressed files (.gz / .Z) are decompressed automatically.

    Parameters
    ----------
    f : str
        Path to the RINEX 3 navigation file (plain or .gz / .Z compressed).

    Returns
    -------
    iono_corr_dic : dict
        Keys are the 4-character correction type (e.g. ``'GPSA'``, ``'GPSB'``,
        ``'GAL'``, ``'BDSA'`` …).  Values are numpy arrays of the numeric
        parameters found on that line.

    time_sys_corr_dic : dict
        Keys are the 4-character correction type (e.g. ``'GAUT'``, ``'GPUT'``,
        ``'GLUT'`` …).  Values are numpy arrays of the numeric parameters
        (2 floats + 2 ints) found on that line.

    Examples
    --------
    >>> iono, tsys = read_rinex_nav_v3_header('BRDC00GOP_R_20191760000_01D_MN.rnx')
    >>> iono['GPSA']
    array([ 4.6566e-09,  1.4901e-08, -5.9605e-08, -1.1921e-07])
    >>> tsys['GAUT']
    array([ 1.9557774067e-08, -1.154631946e-14,  1.72800e+05,  2.059e+03])
    """
    if f.endswith(".gz") or f.endswith(".Z") or f.endswith(".z"):
        fn_use = files_rw.unzip_gz_z(f)
    else:
        fn_use = f

    iono_corr_dic = {}
    time_sys_corr_dic = {}

    with open(fn_use, "r") as f:
        for line in f:
            if "END OF HEADER" in line:
                break

            if "IONOSPHERIC CORR" in line:
                # cols 0-4: key (4 chars), cols 5-59: numeric values (free format)
                key = line[:4].strip()
                vals = np.array(line[5:60].split(), dtype=float)
                iono_corr_dic[key] = vals

            elif "TIME SYSTEM CORR" in line:
                # cols 0-4: key (4 chars), cols 5-59: numeric values (free format)
                # values can be concatenated (e.g. "1.23e-08-4.56e-14"), so use
                # a regex to split on boundaries between numbers
                key = line[:4].strip()
                raw = line[5:60].strip()
                # split on whitespace OR on sign boundary between floats/ints
                tokens = re.findall(r'[+-]?[\d.]+(?:[eE][+-]?\d+)?', raw)
                vals = np.array(tokens, dtype=float)
                time_sys_corr_dic[key] = vals

    return iono_corr_dic, time_sys_corr_dic


def read_rinex_nav(fn, writeh5=None):
    """
    Read a RINEX 2 or RINEX 3 navigation file into a DataFrame.

    Based on Michael Hirsch readRinexNav function.
    http://gage14.upc.es/gLAB/HTML/GPS_Navigation_Rinex_v2.11.html

    The RINEX version is auto-detected from the file header.
    Compressed files (.gz / .Z) are decompressed automatically.

    RINEX3 reading is optimized for GOP files.

    Parameters
    ----------
    fn : str
        Path to the RINEX navigation file (plain or .gz / .Z compressed).
    writeh5 : ignored
        Kept for backward compatibility.

    Returns
    -------
    nav : pd.DataFrame
        Navigation data with columns defined in _NAV_COLS plus 'sys' and 'svni'.
    """
    # Decompress file if needed
    if fn.endswith(".gz") or fn.endswith(".Z") or fn.endswith(".z"):
        fn = files_rw.unzip_gz_z(fn)

    version = None

    with open(fn, "r") as f:
        # Detect version and skip to end of header
        while True:
            line = f.readline()
            if "RINEX VERSION" in line:
                version = 2 if float(line.split()[0]) < 3 else 3
            if "END OF HEADER" in line:
                break

        if version == 2:
            return _read_rinex_nav_v2(f, fn)
        else:
            return _read_rinex_nav_v3(f)
