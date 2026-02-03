#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 09/12/2025 15:36:24

@author: psakic
"""

# Add missing imports required by the module
import datetime as dt
import logging
import os
import struct
import numpy as np
import warnings

# Module logger
log = logging.getLogger('geodezyx')

#  _                         _____                          _
# | |                       / ____|                        | |
# | |     ___  __ _ _ __   | (___   ___  ___ ___  _ __   __| |___
# | |    / _ \/ _` | '_ \   \___ \ / _ \/ __/ _ \| '_ \ / _` / __|
# | |___|  __/ (_| | |_) |  ____) |  __/ (_| (_) | | | | (_| \__ \
# |______\___|\__,_| .__/  |_____/ \___|\___\___/|_| |_|\__,_|___/
#                  | |
#                  |_|

#### LEAP SECONDS MANAGEMENT
def leapseconds_harcoded_list():
    """
    Returns a hardcoded list of leap seconds.

    This function provides a list of tuples where each tuple contains a datetime object
    representing the date and time when a leap second was added, and an integer representing
    the total number of leap seconds that have been added up to that date.

    Returns
    -------
    list of tuple
        A list of tuples where each tuple contains:
        - datetime: The date and time when a leap second was added.
        - int: The total number of leap seconds added up to that date.

    Note
    ----
    The list is hardcoded and may be outdated if the IERS's bulletin C has been updated.
    The last known update is included in the warning message.

    Universal Time (UT1) and International Atomic Time (TAI) were defined as equal
    in 1958. When UTC was introduced in 1972, UT1 had shifted by around 10 seconds
    in relation to TAI.
    We therefore chose an initial offset of 10 seconds between UTC and TAI.

    the initial 10 sec are added by default in frontend find_leapsecond()


    """
    hardcoded_leapsec_lis = [
        (dt.datetime(1972, 7, 1, 0, 0), 1),
        (dt.datetime(1973, 1, 1, 0, 0), 2),
        (dt.datetime(1974, 1, 1, 0, 0), 3),
        (dt.datetime(1975, 1, 1, 0, 0), 4),
        (dt.datetime(1976, 1, 1, 0, 0), 5),
        (dt.datetime(1977, 1, 1, 0, 0), 6),
        (dt.datetime(1978, 1, 1, 0, 0), 7),
        (dt.datetime(1979, 1, 1, 0, 0), 8),
        (dt.datetime(1980, 1, 1, 0, 0), 9),
        (dt.datetime(1981, 7, 1, 0, 0), 10),
        (dt.datetime(1982, 7, 1, 0, 0), 11),
        (dt.datetime(1983, 7, 1, 0, 0), 12),
        (dt.datetime(1985, 7, 1, 0, 0), 13),
        (dt.datetime(1988, 1, 1, 0, 0), 14),
        (dt.datetime(1990, 1, 1, 0, 0), 15),
        (dt.datetime(1991, 1, 1, 0, 0), 16),
        (dt.datetime(1992, 7, 1, 0, 0), 17),
        (dt.datetime(1993, 7, 1, 0, 0), 18),
        (dt.datetime(1994, 7, 1, 0, 0), 19),
        (dt.datetime(1996, 1, 1, 0, 0), 20),
        (dt.datetime(1997, 7, 1, 0, 0), 21),
        (dt.datetime(1999, 1, 1, 0, 0), 22),
        (dt.datetime(2006, 1, 1, 0, 0), 23),
        (dt.datetime(2009, 1, 1, 0, 0), 24),
        (dt.datetime(2012, 7, 1, 0, 0), 25),
        (dt.datetime(2015, 7, 1, 0, 0), 26),
        (dt.datetime(2017, 1, 1, 0, 0), 27)
    ]
    log.warning(
        "Harcoded leap second list loaded\nIt might be wrong if IERS's bulletin C has been updated!\nlast known update: %s",
        hardcoded_leapsec_lis[-1]
    )

    return hardcoded_leapsec_lis


def leapseconds_parse_post2404(leapsec_file_path='/usr/share/zoneinfo/leap-seconds.list'):
    """
    Parses the leap seconds file and returns a list of leap seconds.
    This function is compatible with post Ubuntu 24.04 versions.

    This function reads a file containing leap second information and extracts
    the timestamps and the corresponding leap second values. The timestamps are
    converted to datetime objects.

    Parameters
    ----------
    leapsec_file_path : str, optional
        The path to the leap seconds file. The default is '/usr/share/zoneinfo/leap-seconds.list'.

    Returns
    -------
    list of tuple
        A list of tuples where each tuple contains:
        - datetime: The date and time when a leap second was added.
        - int: The total number of leap seconds added up to that date.

    Note
    ----
    The list is hardcoded and may be outdated if the IERS's bulletin C has been updated.
    The last known update is included in the warning message.

    Universal Time (UT1) and International Atomic Time (TAI) were defined as equal
    in 1958. When UTC was introduced in 1972, UT1 had shifted by around 10 seconds
    in relation to TAI.
    We therefore chose an initial offset of 10 seconds between UTC and TAI.

    the initial 10 sec are added by default in frontend find_leapsecond()
    thus, they are substracted here.
    """
    leapsec_lis = []

    with open(leapsec_file_path, 'r') as file:
        for line in file:
            # Skip header lines
            if line.startswith('#'):
                continue

            # Split the line into components
            parts = line.split()
            if len(parts) < 2:
                continue

            # Extract the timestamp and the leap second value
            timestamp = int(parts[0])
            leapsec = int(parts[1])
            leapsec = leapsec - 10 # 10 sec offset between TAI and UTC, re-added later

            # Convert the timestamp to a datetime object
            # Avoid dependency on external ntp2dt function; convert directly here
            date_leapsec = dt.datetime(1900, 1, 1) + dt.timedelta(seconds=timestamp)

            # Append the parsed data to the list
            leapsec_lis.append((date_leapsec, leapsec))

    return leapsec_lis

def leapseconds_parse_pre2404(leapsec_file_path='/usr/share/zoneinfo/right/UTC'):
    """
    Parses the leap seconds file and returns a list of leap seconds.
    This function is compatible with pre Ubuntu 24.04 versions.

    This function reads a file containing leap second information and extracts
    the timestamps and the corresponding leap second values. The timestamps are
    converted to datetime objects.

    Parameters
    ----------
    leapsec_file_path : str, optional
        The path to the leap seconds file.
        The default is '/usr/share/zoneinfo/right/UTC'.

    Returns
    -------
    list of tuple
        A list of tuples where each tuple contains:
        - datetime: The date and time when a leap second was added.
        - int: The total number of leap seconds added up to that date.

    Notes
    -----
    The list is hardcoded and may be outdated if the IERS's bulletin C has been updated.
    The last known update is included in the warning message.

    Universal Time (UT1) and International Atomic Time (TAI) were defined as equal
    in 1958. When UTC was introduced in 1972, UT1 had shifted by around 10 seconds
    in relation to TAI.
    We therefore chose an initial offset of 10 seconds between UTC and TAI.

    the initial 10 sec are added by default in frontend find_leapsecond()

    Source
    ------
    http://stackoverflow.com/questions/19332902/extract-historic-leap-seconds-from-tzdata

    """

    def _print_leaps(leap_lst):
        """
        INTERNAL_FUNCTION
        """
        # leap_lst is tuples: (timestamp, num_leap_seconds)
        outlist = []
        for ts, num_secs in leap_lst:
            dtime = (dt.datetime.fromtimestamp(ts - num_secs + 1))
            outlist.append((dtime, num_secs))
        return outlist

    f = open(leapsec_file_path, 'rb')

    tzfile_magic = 'TZif'.encode('US-ASCII')

    fmt = ">4s c 15x 6l"
    size = struct.calcsize(fmt)
    (tzfile_magic, tzfile_format, ttisgmtcnt, ttisstdcnt, leapcnt, timecnt,
     typecnt, charcnt) = struct.unpack(fmt, f.read(size))
    # print("DEBUG: tzfile_magic: {} tzfile_format: {} ttisgmtcnt: {} ttisstdcnt: {} leapcnt: {} timecnt: {} typecnt: {} charcnt: {}".format(tzfile_magic, tzfile_format, ttisgmtcnt, ttisstdcnt, leapcnt, timecnt, typecnt, charcnt))

    # Make sure it is a tzfile(5) file
    assert tzfile_magic == tzfile_magic, (
        "Not a tzfile; file magic was: '{}'".format(tzfile_magic))

    # comments below show struct codes such as "l" for 32-bit long integer
    offset = (timecnt * 4  # transition times, each "l"
              + timecnt * 1  # indices tying transition time to ttinfo values, each "B"
              + typecnt * 6  # ttinfo structs, each stored as "lBB"
              + charcnt * 1)  # timezone abbreviation chars, each "c"

    f.seek(offset, 1)  # seek offset bytes from current position

    fmt = '>{}l'.format(leapcnt * 2)
    # print("DEBUG: leapcnt: {}  fmt: '{}'".format(leapcnt, fmt))
    size = struct.calcsize(fmt)
    data = struct.unpack(fmt, f.read(size))

    lst = [(data[i], data[i + 1]) for i in range(0, len(data), 2)]
    assert all(lst[i][0] < lst[i + 1][0] for i in range(len(lst) - 1))
    assert all(lst[i][1] == lst[i + 1][1] - 1 for i in range(len(lst) - 1))

    final_leap_lis = _print_leaps(lst)

    return final_leap_lis


def extract_leapseconds_from_system():
    """
    Extract a leap seconds list from the system (Ubuntu-based systems)
    If the leap seconds file is not found, a hardcoded list is returned

    Note
    ----
    Universal Time (UT1) and International Atomic Time (TAI) were defined as equal
    in 1958. When UTC was introduced in 1972, UT1 had shifted by around 10 seconds
    in relation to TAI.
    We therefore chose an initial offset of 10 seconds between UTC and TAI.

    the initial 10 sec are added by default in frontend find_leapsecond()
    """

    leapsec_fname_post2404 = "/usr/share/zoneinfo/leap-seconds.list"
    leapsec_fname_pre2404 = '/usr/share/zoneinfo/right/UTC'

    leapsec_lis = []
    success_parsing = False

    if os.path.isfile(leapsec_fname_pre2404) and not success_parsing:
        try:
            leapsec_lis = leapseconds_parse_pre2404(leapsec_fname_pre2404)
            success_parsing = True
        except Exception as e:
            log.warning("Error while parsing leap seconds file %s", leapsec_fname_pre2404)
            pass

    if os.path.isfile(leapsec_fname_post2404) and not success_parsing:
        try:
            leapsec_lis = leapseconds_parse_post2404(leapsec_fname_post2404)
            success_parsing = True
        except Exception:
            log.warning("Error while parsing leap seconds file %s", leapsec_fname_post2404)
            pass

    if not success_parsing:
        leapsec_lis = leapseconds_harcoded_list()

    return leapsec_lis

LEAP_SEC_LIS = extract_leapseconds_from_system()

def find_leapsecond(dtin, leapsec_lis_inp=LEAP_SEC_LIS,
                    apply_initial_delta=True):
    """
    Find the TAI-UTC leap second for a given datetime

    Parameters
    ----------
    dtin : datetime.datetime
        Epoch for which the leap second is researched

    leapsec_lis_inp : list, optional
        A list of leap second, provided by extract_leapseconds_from_system()
        automatically determined if given list is empty

    apply_initial_delta : bool, optional
        apply the initial delta of 10 seconds between UTC and TAI
        See note below
        default is True

    Returns
    -------
    leapsec_out : int
        The leap second for the given epoch

    Note
    ----
    Universal Time (UT1) and International Atomic Time (TAI) were defined as equal
    in 1958. When UTC was introduced in 1972, UT1 had shifted by around 10 seconds
    in relation to TAI.
    We therefore chose an initial offset of 10 seconds between UTC and TAI.

    In 1972, the leap-second system was introduced so that the broadcast UTC
    seconds could be made exactly equal to the standard SI second, while still
    maintaining the UTC time of day and changes of UTC date synchronized with
    those of UT1 (the solar time standard that superseded GMT).[11] By then,
    the UTC clock was already 10 seconds behind TAI, which had been
    synchronized with UT1 in 1958, but had been counting true SI seconds
    since then. After 1972, both clocks have been ticking in SI seconds,
    so the difference between their readouts at any time is 10 seconds plus
    the total number of leap seconds that have been applied to UTC
    (37 seconds as of January 2017).
    """

    if len(leapsec_lis_inp) == 0:
        leapsec_lis_use = extract_leapseconds_from_system()
    else:
        leapsec_lis_use = leapsec_lis_inp

    timstp_lis = [e[0] for e in leapsec_lis_use]
    leapsec_lis = [e[1] for e in leapsec_lis_use]

    timstp_lis.append(dtin)
    i_dtin = sorted(timstp_lis).index(dtin)

    leapsec_out = leapsec_lis[i_dtin - 1]
    # -1 because the leap second is the one of the timestamp right before the inserted dtin

    if apply_initial_delta:
        leapsec_out = leapsec_out + 10
    else:
        leapsec_out = leapsec_out

    return leapsec_out

