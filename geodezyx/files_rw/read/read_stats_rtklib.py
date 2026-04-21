#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.files_rw contains reading functions to
import RTKLIB solution status (.stat) files.

it can be imported directly with:
from geodezyx import files_rw

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox

RTKLIB stat file format reference: rtkpos.c (RTKLIB source)

Record types
------------
$POS,week,tow,stat,posx,posy,posz,posxf,posyf,poszf
$VELACC,week,tow,stat,vele,veln,velu,acce,accn,accu,velef,velnf,veluf,accef,accnf,accuf
$CLK,week,tow,stat,clk1,clk2,clk3,clk4,clk5,clk6
$ION,week,tow,stat,sat,az,el,ion,ion_fixed
$TROP,week,tow,stat,rcv,ztd,ztdf
$HWBIAS,week,tow,stat,frq,bias,biasf
$SAT,week,tow,sat,frq,az,el,resp,resc,vsat,snr,fix,slip,lock,outc,slipc,rejc,icbias,bias,bias_var,lambda
"""

########## BEGIN IMPORT ##########
import gzip
import logging
import os
from collections import defaultdict

import numpy as np
import pandas as pd

from geodezyx import conv

log = logging.getLogger("geodezyx")
##########  END IMPORT  ##########


# ---------------------------------------------------------------------------
# Column definitions for each record type (without the leading $TYPE field)
# ---------------------------------------------------------------------------
_STAT_COLS = {
    "POS": [
        "week", "tow", "stat",
        "posx", "posy", "posz",
        "posxf", "posyf", "poszf",
    ],
    "VELACC": [
        "week", "tow", "stat",
        "vele", "veln", "velu",
        "acce", "accn", "accu",
        "velef", "velnf", "veluf",
        "accef", "accnf", "accuf",
    ],
    "CLK": [
        "week", "tow", "stat",
        "clk_idx",
        "clk_gps", "clk_glo", "clk_gal", "clk_bds", "clk_irn", "clk_qzs",
    ],
    "ION": [
        "week", "tow", "stat",
        "sat", "az", "el",
        "ion", "ion_fixed",
    ],
    "TROP": [
        "week", "tow", "stat",
        "rcv", "ztd", "ztdf",
    ],
    "HWBIAS": [
        "week", "tow", "stat",
        "frq", "bias", "biasf",
    ],
    "SAT": [
        "week", "tow",
        "sat", "frq",
        "az", "el",
        "resp", "resc",
        "vsat", "snr",
        "fix", "slip", "lock",
        "outc", "slipc", "rejc",
        "icbias", "bias", "bias_var", "lambda",
    ],
}

# CLK record: the actual file has clk_idx + 6 clk values  (7 numeric fields
# after stat).  RTKLIB writes: stat,clk1,clk2,...,clk6  but the example file
# shows an extra integer field right after stat (likely the receiver index).
# We keep it as clk_idx to be safe.



def read_rtklib_stat(filein, keys=None):
    """Read a RTKLIB solution status (.stat) file.

    Parameters
    ----------
    filein : str
        Path to the .stat file (plain text or gzip-compressed ``.stat.gz``).
    keys : str or list of str, optional
        Record type(s) to parse, e.g. ``'SAT'`` or ``['POS', 'CLK']``.
        If ``None`` (default), all record types present in the file are returned.

    Returns
    -------
    dict
        Dictionary with record-type keys (``'POS'``, ``'VELACC'``,
        ``'CLK'``, ``'SAT'``, ``'ION'``, ``'TROP'``, ``'HWBIAS'``) and
        pandas DataFrame values.  Each DataFrame contains an ``epoch``
        column (Python :class:`datetime.datetime`) computed from the
        ``week`` / ``tow`` GPS time fields, plus all other fields defined
        in the RTKLIB source.  Only record types actually present in the
        file are included in the output dictionary.
    """
    raw = defaultdict(list)

    if keys is not None:
        keys = {keys} if isinstance(keys, str) else set(keys)

    opener = gzip.open if str(filein).endswith(".gz") else open
    with opener(filein, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] != "$":
                continue
            parts = line.split(",")
            rectype = parts[0][1:]  # strip leading '$'
            if keys is not None and rectype not in keys:
                continue
            raw[rectype].append(parts[1:])  # drop the $TYPE field

    result = {}

    for rectype, rows in raw.items():
        cols = _STAT_COLS.get(rectype)

        if cols is None:
            # Unknown record type – store as raw strings with generic names
            log.warning("Unknown RTKLIB stat record type: $%s – storing raw", rectype)
            max_fields = max(len(r) for r in rows)
            cols = [f"field{i}" for i in range(max_fields)]

        # Pad / trim rows to the expected number of columns
        n = len(cols)
        padded = []
        for r in rows:
            if len(r) < n:
                r = r + [""] * (n - len(r))
            padded.append(r[:n])

        df = pd.DataFrame(padded, columns=cols)

        # Cast numeric columns (everything except string columns like 'sat')
        str_cols = {"sat"} if rectype in ("SAT", "ION") else set()
        for col in df.columns:
            if col in str_cols:
                continue
            df[col] = pd.to_numeric(df[col], errors="coerce")

        # Add epoch column
        if "week" in df.columns and "tow" in df.columns:
            df["epoch"] = df.apply(
                lambda row: conv.gpstime2dt(int(row["week"]), float(row["tow"]), dow_input=False),
                axis=1,
            )
            # Move epoch to front
            cols_ordered = ["epoch"] + [c for c in df.columns if c != "epoch"]
            df = df[cols_ordered]

        result[rectype] = df.reset_index(drop=True)

    return result




