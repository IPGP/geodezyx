#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27/08/2025 22:34:53

@author: psakic
"""

import datetime as dt
import os
from io import BytesIO
import numpy as np
import pandas as pd

def read_rinex_nav(fn, writeh5=None, version=2):
    """
    Based on Michael Hirsch readRinexNav function
    Manage RINEX3
    http://gage14.upc.es/gLAB/HTML/GPS_Navigation_Rinex_v2.11.html
    """

    with open(fn, "r") as f:
        # find end of header, which has non-constant length
        while True:

            line = f.readline()

            if "RINEX VERSION" in line:
                if float(line.split()[0]) < 3:
                    version = 2
                else:
                    version = 3

            if "END OF HEADER" in line:
                break

            if version == 2:
                startcol = 3  # column where numerical data starts
                nfloat = 19  # number of text elements per float data number
                nline = 7  # number of lines per record
            elif version == 3:
                startcol = 5  # column where numerical data starts
                nfloat = 19  # number of text elements per float data number
                nline = 7  # number of lines per record

        # handle frame by frame
        sv = []
        epoch = []
        raws = ""

        while True:
            headln = f.readline()
            if not headln:
                break
            # handle the header
            if version == 2:
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
                    + "".join(
                        f.readline()[startcol:].rstrip() for _ in range(nline - 1)
                    )
                    + f.readline()[startcol:79].rstrip()
                )

            if version == 3:
                fields = headln[:23].split(" ")
                sv.append(fields[0])
                year = int(fields[1])

                epoch.append(
                    dt.datetime(
                        year=year,
                        month=int(fields[2]),
                        day=int(fields[3]),
                        hour=int(fields[4]),
                        minute=int(fields[5]),
                        second=int(fields[6]),
                    )
                )

                """
                now get the data.
                Use rstrip() to chomp newlines consistently on Windows and Python 2.7/3.4
                Specifically [:-1] doesn't work consistently as .rstrip() does here.
                """
                raw = (
                    headln[23:].rstrip()
                    + " "
                    + "".join(
                        f.readline()[startcol:].rstrip() + " " for _ in range(nline - 1)
                    )
                    + f.readline()[startcol:80].rstrip()
                )

            raws += raw + "\n"

    raws = raws.replace("D", "E")

    strio = BytesIO(raws.encode())
    darr = np.genfromtxt(strio, delimiter=nfloat)

    nav = pd.DataFrame(
        darr,
        epoch,
        [
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
        ],
    )

    if version == 3:
        const = [e[0] for e in sv]
        svni = [int(e[1:]) for e in sv]
        nav["sys"] = pd.Series(np.array(const), index=nav.index)
        nav["svni"] = pd.Series(np.array(svni), index=nav.index)
    elif version == 2:
        rinexnav_type = os.path.basename(fn)[-1]
        if rinexnav_type == "n":
            nav["sys"] = pd.Series(["G"] * len(nav.index), index=nav.index)
        else:
            nav["sys"] = pd.Series(
                [rinexnav_type.upper()] * len(nav.index), index=nav.index
            )
        nav["svni"] = pd.Series(np.array(sv), index=nav.index)

    if False:
        h5fn = fn.with_suffix(".h5")
        print(("saving NAV data to {}".format(h5fn)))
        nav.to_hdf(h5fn, key="NAV", mode="a", complevel=6, append=False)

    return nav
