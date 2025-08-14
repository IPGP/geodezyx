#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:45:36 2024

@author: psakic
"""

import datetime as dt
import logging
import os

import pandas as pd
import xarray as xr

from geodezyx import operational, conv
from geodezyx.marine import obp

log = logging.getLogger('geodezyx')


def dac_get_aviso(date_srt_inp, date_end_inp, out_dir, user, passwd):
    """
    Downloads DAC (Doppler-derived Atmospheric Correction) files from the AVISO FTP server for a given date range.

    Parameters
    ----------
    date_srt_inp : datetime
        The start date of the range for which DAC files are to be downloaded.
    date_end_inp : datetime
        The end date of the range for which DAC files are to be downloaded.
    out_dir : str
        The local directory where the downloaded DAC files should be saved.
    user : str
        The username for the FTP server.
    passwd : str
        The password for the FTP server.

    Returns
    -------
    list
        A list of paths to the downloaded DAC files.

    """

    date_rng = pd.date_range(date_srt_inp, date_end_inp, freq="6H")

    ftpbase = "ftp://ftp-access.aviso.altimetry.fr/auxiliary/dac/dac_delayed_global/"
    p_stk = []
    for d in date_rng:
        jjul_cnes = conv.dt2jjul_cnes(d)
        dac_name = "dac_dif_" + str(jjul_cnes) + "_" + str(d.hour).zfill(2) + ".nc"

        p = ftpbase + str(d.year) + "/" + dac_name

        p_stk.append(p)

    operational.ftp_downld_front(p_stk, out_dir,
                                 user=user,
                                 passwd=passwd,
                                 secure_ftp=False)

    return


def dac_extract(dac_files_lis,
                mlon=45.46,
                mlat=-12.78,
                epoch_min=dt.datetime(1980, 1, 1),
                epoch_max=dt.datetime(2099, 1, 1)):
    """
    Extracts DAC (Doppler-derived Atmospheric Correction) values from a list
    of DAC files for a given location and date range.

    Parameters
    ----------
    dac_files_lis : list
        The list of DAC files from which the DAC values are to be extracted.
    mlon : float, optional
        The longitude of the location for which the DAC values are to be extracted. Default is 45.46.
    mlat : float, optional
        The latitude of the location for which the DAC values are to be extracted. Default is -12.78.
    epoch_min : datetime, optional
        The start date of the range for which the DAC values are to be extracted. Default is January 1, 1980.
    epoch_max : datetime, optional
        The end date of the range for which the DAC values are to be extracted. Default is January 1, 2099.

    Returns
    -------
    DataFrame
        A DataFrame with the extracted DAC values. The index of the DataFrame is the date and the column contains
        the DAC values.
    """

    if not dac_files_lis:
        log.error("No DAC files provided.")
        raise FileNotFoundError("No DAC files provided.")

    ds_stk = []
    dac_stk = []

    for f in dac_files_lis:
        log.debug("read DAC file: %s", f)
        ds = xr.open_dataset(f)
        bn = os.path.basename(f)
        # tname = conv.jjul_cnes2dt(bn[8:13]) + dt.timedelta(hours = int(bn[14:16]))
        dac_out = obp.interp_xy(ds, x=mlon, y=mlat).dac

        dac_val = float(dac_out.values)
        t = dac_out
        t = conv.str_date2dt(dac_out.date)

        ds_stk.append(ds)
        dac_stk.append((t, dac_val))

    df_dac_raw = pd.DataFrame(dac_stk)
    df_dac_raw.set_index(0, inplace=True)

    if epoch_max and epoch_min:
        df_dac = df_dac_raw.loc[epoch_min:epoch_max, 1]
    else:
        df_dac = df_dac_raw

    return df_dac
