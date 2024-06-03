#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 17:41:32 2022

@author: psakic

This module regroups the functions for the resolution estimation of the 
REVOSIMA's OBSCOM-embeded pressure sensor
"""

import itertools
#### Import the logger
import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

import geodezyx.conv as conv
from geodezyx import utils

log = logging.getLogger(__name__)


#  ______               _                 _    __                  _   _
# |  ____|             | |               | |  / _|                | | (_)
# | |__ _ __ ___  _ __ | |_ ___ _ __   __| | | |_ _   _ _ __   ___| |_ _  ___  _ __  ___
# |  __| '__/ _ \| '_ \| __/ _ \ '_ \ / _` | |  _| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# | |  | | | (_) | | | | ||  __/ | | | (_| | | | | |_| | | | | (__| |_| | (_) | | | \__ \
# |_|  |_|  \___/|_| |_|\__\___|_| |_|\__,_| |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/


def obscom_import_pres_df(
        csv_pres_inp,
        dic_coeffs,
        freq_clk=4.096e6,
        count_sensor_temp=168e3,
        count_sensor_pres=30e3,
        integ_timespan_temp=10,
        integ_timespan_pres=100,
        epoc_idx=True,
        out_pres="pa",
):
    """
    Imports pressure data from a CSV file and computes the frequency, temperature, pressure, and depth.

    Parameters
    ----------
    csv_pres_inp : str
        The path to the CSV file containing the pressure data.
    dic_coeffs : dict
        The dictionary of coefficients used for the computations.
    freq_clk : float, optional
        The frequency of the clock in Hertz. Default is 4.096e6.
    count_sensor_temp : float, optional
        The count of the temperature sensor. Default is 168e3.
    count_sensor_pres : float, optional
        The count of the pressure sensor. Default is 30e3.
    integ_timespan_temp : int, optional
        The integration timespan for the temperature sensor. Default is 10.
    integ_timespan_pres : int, optional
        The integration timespan for the pressure sensor. Default is 100.
    epoc_idx : bool, optional
        If True, sets the epoch as the index of the DataFrame. Default is True.
    out_pres : str, optional
        The unit of the output pressure. Default is "pa".

    Returns
    -------
    DataFrame
        A DataFrame with the epoch, frequency, temperature, pressure, and depth.

    Notes
    -----
    This function first reads the CSV file and creates a DataFrame.
    It then computes the frequency from the counter for both the temperature and pressure sensors.
    The temperature is computed from the frequency using the coefficients from the dictionary.
    The pressure is computed from the frequency of both the pressure and temperature sensors.
    The depth is also computed from the frequency of both the pressure and temperature sensors.
    If epoc_idx is True, the epoch is set as the index of the DataFrame.
    """

    dic_coeffs_temp = _get_coeffs_temp(dic_coeffs)

    pres_df_raw = pd.read_csv(csv_pres_inp)

    pres_df = pd.DataFrame()

    pres_df["epoc"] = conv.posix2dt(pres_df_raw["Timestamp(unixtime)"])

    ### rename columns
    pres_df["cnt_temp"] = pres_df_raw["Temperature"]
    pres_df["cnt_pres"] = pres_df_raw["Pressure"]

    ### remove empty values
    pres_df.replace(0, np.nan, inplace=True)
    pres_df.dropna(inplace=True)

    ### compute frequency from counter
    pres_df["frq_temp"] = counter2freq(
        pres_df["cnt_temp"],
        count_sensor=count_sensor_temp,
        freq_clk=freq_clk,
        integ_timespan=integ_timespan_temp,
    )

    pres_df["frq_pres"] = counter2freq(
        pres_df["cnt_pres"],
        count_sensor=count_sensor_pres,
        freq_clk=freq_clk,
        integ_timespan=integ_timespan_pres,
    )
    ### compute temperature from frequency
    pres_df["temp"] = freq2temp(pres_df["frq_temp"], **dic_coeffs_temp)

    ### compute pressure from frequency, using temperature frequency
    pres_df["pres"] = freq2pres(
        pres_df["frq_pres"],
        pres_df["frq_temp"],
        out=out_pres,
        inp_temp="freq",
        **dic_coeffs
    )

    ### compute pressure from frequency, using temperature itself
    ### (for debug only!)
    debug = False
    if debug:
        pres_df["pres2"] = freq2pres(pres_df["frq_pres"],
                                     pres_df["temp"],
                                     out="pa",
                                     inp_temp="temp",
                                     **dic_coeffs)

    ### compute depth
    pres_df["dept"] = freq2pres(
        pres_df["frq_pres"], pres_df["temp"], out="meter", **dic_coeffs
    )

    if epoc_idx:
        pres_df.set_index("epoc", inplace=True)

    return pres_df


def obscom_plot_pres_df(pres_df_in, plot_count=False):
    """
    Plots the pressure and temperature data from a DataFrame.

    Parameters
    ----------
    pres_df_in : DataFrame
        The DataFrame containing the pressure and temperature data. The index of the DataFrame should be the time.
    plot_count : bool, optional
        If True, plots the count of the temperature and pressure sensors. If False, plots the pressure in pa and temperature in °C. Default is False.

    Returns
    -------
    Figure
        The Figure object with the plot of the pressure and temperature data.

    Notes
    -----
    This function creates a plot with two y-axes. The left y-axis corresponds to the pressure data and the right y-axis corresponds to the temperature data.
    The x-axis corresponds to the time.
    The function uses different colors for the pressure and temperature data for better distinction.
    """

    fig, ax1 = plt.subplots()

    if not plot_count:
        col1 = 'pres'
        col2 = 'temp'
        lgd1 = 'pressure (pa)'
        lgd2 = 'temperature (°C)'

    else:
        col1 = 'cnt_temp'
        col2 = 'cnt_pres'
        lgd1 = 'count (#)'
        lgd2 = 'count (#)'

    ax1.set_xlabel('time')
    ax1.set_ylabel(lgd1, color='C1')
    ax1.plot(pres_df_in.index, pres_df_in[col1], color='C1')
    ax1.tick_params(axis='y', labelcolor='C1')

    # Adding Twin Axes
    ax2 = ax1.twinx()
    ax2.plot(pres_df_in.index, pres_df_in[col2], color='C2')
    ax2.set_ylabel(lgd2, color='C2')
    ax2.tick_params(axis='y', labelcolor='C2')

    return fig


#                  __  __ _      _            _                                      _             _
#                 / _|/ _(_)    (_)          | |         ___                        | |           | |
#   ___ ___   ___| |_| |_ _  ___ _  ___ _ __ | |_ ___   ( _ )     ___ ___  _ __  ___| |_ ___ _ __ | |_ ___
#  / __/ _ \ / _ \  _|  _| |/ __| |/ _ \ '_ \| __/ __|  / _ \/\  / __/ _ \| '_ \/ __| __/ _ \ '_ \| __/ __|
# | (_| (_) |  __/ | | | | | (__| |  __/ | | | |_\__ \ | (_>  < | (_| (_) | | | \__ \ ||  __/ | | | |_\__ \
#  \___\___/ \___|_| |_| |_|\___|_|\___|_| |_|\__|___/  \___/\/  \___\___/|_| |_|___/\__\___|_| |_|\__|___/
#

## Wikipedia values
#  https://en.wikipedia.org/wiki/Pound_per_square_inch

PSI_PER_METER = 1.4503773773022  ## 1 METER = 1 DBAR
PSI_TO_BAR = 0.06894757293168
PSI_TO_MBAR = PSI_TO_BAR * 1000
PSI_TO_PA = 6894.757293168


def get_coeffs(sens_type="all", sens_id=158073):
    if sens_id == 0:
        ### from exemple sheets
        dic_coeffs = {
            "U0": 5.766353,
            "Y1": -4025.183,
            "Y2": -11970.76,
            "Y3": 0.0,
            "T1": 29.89307,
            "T2": 0.337614,
            "T3": 56.99511,
            "T4": 157.6942,
            "T5": 0.0,
            "C1": -22682.65,
            "C2": -1143.743,
            "C3": 70903.62,
            "D1": 0.040903,
            "D2": 0.0,
        }

    elif sens_id == 158073 or sens_id == "T3":
        ### for the 1st OBSCOM
        dic_coeffs = {
            "U0": 5.798999,
            "Y1": -3874.954,
            "Y2": -10166.55,
            "Y3": 0,
            "C1": -25657.25,
            "C2": -645.8024,
            "C3": 73515.96,
            "D1": 0.039737,
            "D2": 0,
            "T1": 30.00177,
            "T2": 0.723913,
            "T3": 53.84611,
            "T4": 147.1236,
            "T5": 0,
        }

    elif sens_id == 158073999 or sens_id == "T3gross":  ### 999 suffix = gross
        ### for the 1st OBSCOM
        dic_coeffs = {
            "U0": 5.799,
            "Y1": -3874.95,
            "Y2": -10166.5,
            "Y3": 0,
            "C1": -25657.2,
            "C2": -645.802,
            "C3": 73516,
            "D1": 0.0397368,
            "D2": 0,
            "T1": 30.0018,
            "T2": 0.723913,
            "T3": 53.8461,
            "T4": 147.124,
            "T5": 0,
        }

    elif sens_id == 158076 or sens_id == "T2":
        ### for the 2nd OBSCOM (T2)
        dic_coeffs = {
            "U0": 5.765270,
            "Y1": -3994.585,
            "Y2": -10705.43,
            "Y3": 0,
            "C1": -25561.86,
            "C2": -230.0508,
            "C3": 79516.88,
            "D1": 0.040172,
            "D2": 0,
            "T1": 30.13945,
            "T2": 0.977439,
            "T3": 56.64988,
            "T4": 158.4542,
            "T5": 0,
        }

    elif sens_id == 158076999 or sens_id == "T2gross":  ### 999 suffix = gross
        ### for the 2nd OBSCOM (T2)
        dic_coeffs = {
            "U0": 5.76527,
            "Y1": -3994.58,
            "Y2": -10705.4,
            "Y3": 0,
            "C1": -25561.9,
            "C2": -230.051,
            "C3": 79516.9,
            "D1": 0.0401721,
            "D2": 0,
            "T1": 30.1394,
            "T2": 0.977439,
            "T3": 56.6499,
            "T4": 158.454,
            "T5": 0,
        }

    elif sens_id == 158081 or sens_id == "T4":
        ### for the 2nd OBSCOM (T2)
        dic_coeffs = {
            "U0": 5.797503,
            "Y1": -4022.350,
            "Y2": -12679.85,
            "Y3": 0,
            "C1": -22175.21,
            "C2": -1252.800,
            "C3": 83574.76,
            "D1": 0.037913,
            "D2": 0,
            "T1": 30.23989,
            "T2": 0.124474,
            "T3": 64.95771,
            "T4": 229.3125,
            "T5": 0,
        }

    if sens_type == "all":
        return dic_coeffs
    elif sens_type == "f0":
        return _get_coeffs_f0(dic_coeffs)
    elif sens_type == "temp":
        return _get_coeffs_temp(dic_coeffs)
    else:
        log.error("check sens_type ('pres' or 'temp')")


def _get_coeffs_temp(dic_coeffs_in):
    """
    return only the useful coefficients to compute temp
    """

    return {c: dic_coeffs_in[c] for c in ["U0", "Y1", "Y2", "Y3"]}


def _get_coeffs_f0(dic_coeffs_in):
    """
    return only the useful coefficients to compute f0
    """

    return {c: dic_coeffs_in[c] for c in ["U0", "T1", "T2", "T3", "T4", "T5"]}


def unitary_tests(sens_id=158073):
    """
    return the correct values for pressure/temperature and
    intermediates parameters @ ftemp=172600.0 Hz and fpres=36300.0 Hz
    with the cofficients of the sensor #158073
    (1st OBSCOM)

    Returns
    -------
    ftemp (Hz), fpres (Hz), temp (°C), f0 (Hz), p (psi)
    """

    if sens_id == 158073:
        ftemp = 172600.0  # Hz
        fpres = 36300.0  # Hz
        temp = 20.090562800024895  # °C
        f0 = 33333.93215844317  # Hz # @ 20.090562800024895°C
        p = 4803.3285794411595  # psi # @ 20.090562800024895 °C
    return ftemp, fpres, temp, f0, p


#  _                                      _
# | |                                    | |
# | |_ ___ _ __ ___  _ __   ___ _ __ __ _| |_ _   _ _ __ ___
# | __/ _ \ '_ ` _ \| '_ \ / _ \ '__/ _` | __| | | | '__/ _ \
# | ||  __/ | | | | | |_) |  __/ | | (_| | |_| |_| | | |  __/
#  \__\___|_| |_| |_| .__/ \___|_|  \__,_|\__|\__,_|_|  \___|
#                   | |
#                   |_|


def freq2U(ft_in, U0):
    """
    Convert temperature sensor frequency to coefficient U

    Parameters
    ----------
    ft_in : float
        frequency of the temperature sensor in Hertz.
    U0 : float
        input coefficient U0.

    Returns
    -------
    float
        coefficient U.

    """
    return 1e6 / ft_in - U0


def freq2temp(ft_in, U0, Y1, Y2, Y3):
    """
    Convert temperature sensor frequency to temperature

    Parameters
    ----------
    ft_in :
        frequency of the temperature sensor in Hertz.
    U0, Y1-3 : float
        Calibration coefficients provided by PAROS.

    Returns
    -------
    float
        temperature in Celsius.

    """

    if utils.is_iterable(ft_in):
        return np.array([freq2temp(e, U0, Y1, Y2, Y3) for e in ft_in])
    else:
        if ft_in == 0.0:
            ft_in = np.nan

        # X = (1/ft_in) * 10**6 # (1/30000) * 10**6
        # U = X - U0

        U = freq2U(ft_in, U0)

        Temp = Y1 * U + Y2 * U ** 2 + Y3 * U ** 3

        return Temp


def temp2freq(t_in, U0, Y1, Y2, Y3, out="freq", force_root=False):
    """
    Convert temperature to temperature sensor frequency

    Parameters
    ----------
    t_in :
        temperature in Celsius.
    U0, Y1-3 : float, optional
        Calibration coefficients provided by PAROS.
    out : float, optional
        freqency ('freq') or period ('tau') or coefficient U ('U')
        of the sensor. The default is 'freq'.
    force_root : bool or int
        force the finded root for U
        if integer (1 or 2), use the first or second root
        default is False

    Returns
    -------
    float
        temperature sensor frequency or period.

    Note
    ----
    This function is unstable because of the root finding

    You can force it with `force_root`

    """

    #### Handle T as an interable (array, list....)
    if utils.is_iterable(t_in):
        return np.array([temp2freq(t, U0, Y1, Y2, Y3, out, force_root) for t in t_in])
    #### Handle T as a scalar
    else:
        Uraw = np.roots([Y3, Y2, Y1, -t_in]) + U0
        U = Uraw
        F = 10 ** 6 * U ** -1  ### converted to frequency

        # Fbool = np.logical_and((F < 176000) , (F > 168000))

        Fok = np.min(F)
        Uok = U[Fok == F][0]

        if False and not ((10 ** 6 / U < 176000).all() and (10 ** 6 / U > 168000).all()):
            log.warn("temp U/freq out of range: %f/%f", U, 10 ** 6 / U)

        if out == "freq":
            #            return [ 10**6/u  for u in U  if ( (10**6/u< 176000) and (10**6/u > 168000 ))][0]
            # return [ 10**6/u  for u in U ][ilist]
            # outvals = 10**6 / U
            outvals = Fok
        elif out == "tau":
            # return [ 10**6/u  for u in U  if ( (10**6/u< 176000) and (10**6/u > 168000 ))][ilist]**-1
            # return [ 10**6/u  for u in U ][ilist]**-1
            # outvals = U / 10**6
            outvals = Fok ** -1
        elif out == "U":
            #            return [ u-U0  for u in U  if ( (u> 10**6/176000) and (u < 10**6/168000 ))][0]
            # return [ u-U0  for u in U ][ilist]
            # outvals = Uok - U0
            outvals = freq2U(Fok, U0)
        else:
            log.error("bad out value")
            raise Exception("bad output format")

    return outvals


def _temp2pres_f0(t_or_ft_in, U0, Y1, Y2, Y3, T1, T2, T3, T4, T5, inp_temp="temp"):
    """
    compute the f0 value, i.e. the temperature component for
    the pressure computation

    f0 is homogeneous to a frequency

    Internal function for ``freq2pres``

    Parameters
    ----------
    t_or_ft_in : float
        temperature in Celsius or temperature sensor frequency in Hertz.
        depends on ``inp_temp``
    U0, T1-5 : float, optional
        Calibration coefficients provided by PAROS.
    inp_temp : str
        'freq' temperature sensor frequency in Hertz
        'temp' temperature in Celsius
        default is 'temp'

    Returns
    -------
    float
        f0, the temperature component for the pressure computation.

    Note
    ----
    Using temperature as input is unstable because of the root finding for U

    frequency as input is recommended
    """
    if inp_temp == "temp":
        F = temp2freq(t_or_ft_in, U0, Y1, Y2, Y3, out="freq")
    else:
        F = t_or_ft_in

    # U= (1/F)*10**6 - U0
    U = freq2U(F, U0)

    T0 = T1 + T2 * U + T3 * U ** 2 + T4 * U ** 3 + T5 * U ** 4

    F0 = 10 ** 6 / T0

    return F0


#  _ __  _ __ ___  ___ ___ _   _ _ __ ___
# | '_ \| '__/ _ \/ __/ __| | | | '__/ _ \
# | |_) | | |  __/\__ \__ \ |_| | | |  __/
# | .__/|_|  \___||___/___/\__,_|_|  \___|
# | |
# |_|


def freq2pres(
        fp_in,
        t_or_ft_in,
        C1,
        C2,
        C3,
        D1,
        D2,
        U0,
        Y1,
        Y2,
        Y3,
        T1,
        T2,
        T3,
        T4,
        T5,
        out="psi",
        inp_temp="temp",
):
    """
    Convert pressure sensor frequency to pressure or equivalent distance

    Parameters
    ----------
    fp_in :
        frequency of the pressure sensor in Hertz.
    t_or_ft_in : float
        temperature in Celsius or temperature sensor frequency in Hertz.
        depends on ``inp_temp``
    C1 C2 C3 D1 D2 : float
        Calibration coefficients provided by PAROS.
    U0,Y1,Y2,Y3,T1,T2,T3,T4,T5 : float
        Calibration coefficients provided by PAROS for the temperature sensor.
    out : str, optional
        output as pressure ('psi', 'bar', 'mbar', 'pa') or distance ('meter').
        The default is "psi".
    inp_temp : str
        'freq' temperature sensor frequency in Hertz,
        'temp' temperature in Celsius,
        default is 'temp'

    Returns
    -------
    float
        pressure or equivalent distance.

    Note
    ----
    Using temperature as input is unstable because of the root finding of U
    frequency as input is recommended

    """

    if inp_temp == "temp":
        U = temp2freq(t_or_ft_in, U0, Y1, Y2, Y3, out="U")
    elif inp_temp == "freq":
        U = freq2U(t_or_ft_in, U0)

    C = C1 + C2 * U + C3 * U ** 2
    D = D1 + D2 * U

    T = (1 / fp_in) * 10 ** 6

    ### takes inp_temp = 'freq' or 'temp'
    F0 = _temp2pres_f0(
        t_or_ft_in, U0, Y1, Y2, Y3, T1, T2, T3, T4, T5, inp_temp=inp_temp
    )
    T0 = (1 / F0) * 10 ** 6

    P = C * (1 - (T0 ** 2) / (T ** 2)) * (1 - D * (1 - (T0 ** 2) / (T ** 2)))

    if out == "psi":
        return P
    elif out == "bar":
        return P * PSI_TO_BAR
    elif out == "mbar":
        return P * PSI_TO_MBAR
    elif out == "pa":
        return P * PSI_TO_PA
    elif out == "meter":
        return P / PSI_PER_METER
    else:
        log.error("bad out value")
        raise Exception("bad output format")


def pres2freq(
        p_in,
        t_in,
        C1,
        C2,
        C3,
        D1,
        D2,
        U0,
        Y1,
        Y2,
        Y3,
        T1,
        T2,
        T3,
        T4,
        T5,
        inp="psi",
        return_optimize_object=False,
):
    """
    Convert pressure to  pressure sensor frequency

    There is no analytic solution, thus an root find is required

    Parameters
    ----------
    p_in :
        pressure in psi.
    t_in :
        temperature in Celsius.
    inp : optional
        input for p_in as pressure ('psi', 'bar', 'mbar', 'pa')
        or distance ('meter').
        The default is 'psi'.
    return_optimize_object : bool, optional
        if True, return the object of scipy's optimize.root function.
        if False, return the root value
        The default is False.

    Returns
    -------
    float
        pressure sensor frequency (Hertz).

    Note
    ----

    """

    #### Handle p_in as an interable (array, list....)
    if utils.is_iterable(t_in):
        log.err("t_in has to be a scalar")
        raise Exception

    if utils.is_iterable(p_in):
        return np.array(
            [
                pres2freq(
                    p_in_i,
                    t_in,
                    C1,
                    C2,
                    C3,
                    D1,
                    D2,
                    U0,
                    Y1,
                    Y2,
                    Y3,
                    T1,
                    T2,
                    T3,
                    T4,
                    T5,
                    inp,
                    return_optimize_object,
                )
                for p_in_i in p_in
            ]
        )

    #### Handle p_in as a scalar
    else:
        if inp == "psi":
            puse = p_in
        elif inp == "bar":
            puse = p_in / PSI_TO_BAR
        elif inp == "mbar":
            puse = p_in / PSI_TO_MBAR
        elif inp == "pa":
            puse = p_in / PSI_TO_PA
        elif inp == "meter":
            puse = p_in * PSI_PER_METER

        def wrap_zero(fp):
            return (
                    freq2pres(
                        fp, t_in, C1, C2, C3, D1, D2, U0, Y1, Y2, Y3, T1, T2, T3, T4, T5
                    )
                    - puse
            )

        Opti = scipy.optimize.root(wrap_zero, 30000)

        if return_optimize_object:
            return Opti
        else:
            return Opti.x[0]


#                        _
#                       | |
#   ___ ___  _   _ _ __ | |_ ___ _ __
#  / __/ _ \| | | | '_ \| __/ _ \ '__|
# | (_| (_) | |_| | | | | ||  __/ |
#  \___\___/ \__,_|_| |_|\__\___|_|


def freq2counter(
        freq_or_tau_sensor,
        count_sensor,
        freq_or_tau_clk,
        integ_timespan=1,
        inp="freq",
        round_fct=np.floor,
):
    """
    Convert sensor frequency to count number

    Parameters
    ----------
    freq_or_tau_sensor :
        freqency ('freq') or period ('tau') of the sensor.
    count_sensor : int or float
        number of clock periods counted for one sensor count
    freq_or_tau_clk :
        primary frequency/period of the clock that count period..
    integ_timespan : int or float, optional
        integration timespan in seconds. The default is 1.
    inp : str, optional
        input 'freq' for freqency or 'tau' for period (=1/frequency).
        The default is 'freq'.
    round_fct : function, optional
        the function that round the count values. The default is np.floor.

    Returns
    -------
    n_count :
        number of counts counted by the sensor..

    Note
    ----
    usually, ``count_sensor=30e3``  for pressure, ``count_sensor=168e3``
    for temperature and ``freq_clk=4.096e6``

    """

    if inp == "freq":
        tau_sensor = 1 / freq_or_tau_sensor
        tau_clk = 1 / freq_or_tau_clk
    elif inp == "tau":
        tau_sensor = freq_or_tau_sensor
        tau_clk = freq_or_tau_clk
    else:
        log.error("bad inp value")
        raise Exception("bad inp format")

    n_count = (integ_timespan * tau_sensor * count_sensor) / tau_clk

    if round_fct:
        return round_fct(n_count)
    else:
        return n_count


def counter2freq(
        n_counted_by_clk, count_sensor, freq_clk, integ_timespan=1, out="freq"
):
    """
    Convert count number to sensor frequency

    Parameters
    ----------
    n_counted_by_clk : int or float
        number of counts counted by the sensor.
    count_sensor : int or float
        number of clock periods counted for one sensor count
    freq_clk : int or float
        primary frequency of the clock that count period.
    integ_timespan : int or float, optional
        integration timespan in seconds. The default is 1.
    out : str, optional
        output 'freq' for freqency or 'tau' for period (=1/frequency).
        The default is 'freq'.

    Returns
    -------
    float
        freqency ('freq') or period ('tau') of the sensor

    Note
    ----
    usually, ``count_sensor=30e3``  for pressure, ``count_sensor=168e3``
    for temperature and ``freq_clk=4.096e6``

    """

    tau_sensor = n_counted_by_clk / (freq_clk * integ_timespan * count_sensor)

    if out == "freq":
        return 1 / tau_sensor
    elif out == "tau":
        return tau_sensor
    else:
        log.error("bad out value")
        raise Exception("bad output format")


####### HIGHER LEVEL FCTS #####################################################


def pres_resolution(
        val_presin,
        val_temp_in,
        val_clkin,
        count_presin,
        count_temp_in,
        input_val="freq",
        output_val="meter",
        relative_delta=True,
):
    if input_val == "freq":
        if np.any(val_temp_in < 1) or np.any(val_presin < 1) or np.any(val_clkin < 1):
            log.warn("some vals are < 1 while input is freq. Please double check")
        tau_temp = 1 / val_temp_in
        tau_pres = 1 / val_presin
        tau_clk = 1 / val_clkin
    elif input_val == "tau":
        if np.any(val_temp_in > 1) or np.any(val_presin > 1) or np.any(val_clkin > 1):
            log.warn("some vals are > 1 while input is tau. Please double check")
        tau_temp = val_temp_in
        tau_pres = val_presin
        tau_clk = val_clkin

    # internal values are tau, because at least it simplifies the computation
    # of pres_bias  and avoids the 1/(a+b) = 1/a + 1/b mistake

    freq_clk = 1 / tau_clk

    temp_central = freq2temp(1 / tau_temp)
    pres_central = freq2pres(1 / tau_pres, temp_central)

    tau_pres_bias = 2 / (freq_clk * count_presin)
    tau_temp_bias = 2 / (freq_clk * count_temp_in)

    PresResArr = np.zeros((3, 3))

    if relative_delta:
        divider = pres_central
    else:
        divider = 1.0

    for i, j in [(1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1), (1, 1), (0, 1)]:
        pres_biased = freq2pres(
            1 / (tau_pres + i * tau_pres_bias),
            freq2temp(1 / (tau_temp + j * tau_temp_bias)),
        )

        pres_res = np.abs(pres_biased - pres_central) / divider

        PresResArr[i + 1, j + 1] = pres_res

    if output_val == "psi":
        pass
    elif output_val == "meter":
        pres_central = pres_central / PSI_PER_METER
        if not relative_delta:
            PresResArr = (
                    PresResArr / PSI_PER_METER
            )  ### correction only for absolute values !!!
    else:
        raise Exception("bad output format in pres_resolution")

    return pres_central, temp_central, PresResArr.max(), PresResArr


def resolution_grid_compute(
        Tau_pres,
        Tau_temp,
        fe,
        count_pres,
        count_temp,
        output_val="meter",
        relative_delta=False,
):
    PresVals = np.zeros((len(Tau_pres), len(Tau_temp)))
    TempVals = np.zeros((len(Tau_pres), len(Tau_temp)))
    ResVals = np.zeros((len(Tau_pres), len(Tau_temp)))
    FullResVals = np.zeros((len(Tau_pres), len(Tau_temp), 3, 3))

    for itaupres, itautemp in itertools.product(
            range(len(Tau_pres)), range(len(Tau_temp))
    ):
        taupres, tautemp = Tau_pres[itaupres], Tau_temp[itautemp]

        OutRes = pres_resolution(
            taupres,
            tautemp,
            1 / fe,
            count_pres,
            count_temp,
            input_val="tau",
            output_val=output_val,
            relative_delta=relative_delta,
        )

        PresVals[itaupres, itautemp] = OutRes[0]
        TempVals[itaupres, itautemp] = OutRes[1]
        ResVals[itaupres, itautemp] = OutRes[2]
        FullResVals[itaupres, itautemp, :, :] = OutRes[3]

    return PresVals, TempVals, ResVals, FullResVals


############ PLOT OF THE RESOLUTION VALUES
def resolution_plot_as_gradient_grid(
        PresVals,
        TempVals,
        ResVals,
        temp_ref,
        freq_input=False,
        coef_res=1000,
        relative_delta=False,
):
    #### coef_res = 10**6 if PPM
    #### coef_res = 10**3 if error in meters

    #################### INITIALISATION ##########################

    secx_pres2freq = lambda x: pres2freq(x, temp_ref, "meter")
    secx_freq2pres = lambda x: freq2pres(x, temp_ref, "meter")
    secy_temp2freq = lambda x: temp2freq(x, "freq")
    secy_freq2temp = lambda x: freq2temp(x)

    if not freq_input:
        secx_fcts = (secx_pres2freq, secx_freq2pres)
        secy_fcts = (secy_temp2freq, secy_freq2temp)
        xlab = "Depth (m)"
        ylab = "Temp (°C)"
        secxlab = "Pres Sensor Freq (Hz)"
        secylab = "Temp Sensor Freq (Hz)"
    else:
        secx_fcts = (secx_freq2pres, secx_pres2freq)
        secy_fcts = (secy_freq2temp, secy_temp2freq)
        secxlab = "Depth (m)"
        secylab = "Temp (°C)"
        xlab = "Pres Sensor Freq (Hz)"
        ylab = "Temp Sensor Freq (Hz)"

    if not relative_delta:
        ResVals = ResVals * 10 ** 3
        title = "Error in millimeter max={:.4f}mm".format(ResVals.max())
        format_colorbar = "%.4f"
    else:
        ResVals = ResVals * 10 ** 6
        title = "PPM max={:.4f}".format(ResVals.max())
        format_colorbar = "%.4f"

    ### change values of ResVal

    #################### INITIALISATION ##########################

    fig, ax = plt.subplots()
    ax.ticklabel_format(useOffset=False)
    cm = plt.cm.get_cmap("viridis")
    fig.subplots_adjust(right=0.75, top=0.85)
    im1 = ax.contourf(PresVals, TempVals, ResVals, cmap=cm, levels=200)
    cbar_ax1 = fig.add_axes([0.88, 0.15, 0.04, 0.7])
    cbar1 = fig.colorbar(im1, cax=cbar_ax1, format=format_colorbar)

    secx = ax.secondary_xaxis("top", functions=secx_fcts)
    secy = ax.secondary_yaxis("right", functions=secy_fcts)
    secx.ticklabel_format(useOffset=False)
    secy.ticklabel_format(useOffset=False)

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    secx.set_xlabel(secxlab)
    secy.set_ylabel(secylab)
    ax.set_title(title)

    return fig, ax


#   __                  _   _                                                             _
#  / _|                | | (_)                                                           | |
# | |_ _   _ _ __   ___| |_ _  ___  _ __     __ _ _ __ __ ___   _____ _   _  __ _ _ __ __| |
# |  _| | | | '_ \ / __| __| |/ _ \| '_ \   / _` | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
# | | | |_| | | | | (__| |_| | (_) | | | | | (_| | | | (_| |\ V /  __/ |_| | (_| | | | (_| |
# |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|  \__, |_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
#                                            __/ |                     __/ |
#                                           |___/                     |___/


def freq2pres_old(
        f, U, F0, C1=-22682.65, C2=-1143.743, C3=70903.62, D1=0.040903, D2=0.0
):
    """
    freqence capteur pression > pression
    OLD version
    """

    C = C1 + C2 * U + C3 * U ** 2
    D = D1 + D2 * U

    T = (1 / f) * 10 ** 6
    T0 = (1 / F0) * 10 ** 6

    P = C * (1 - (T0 ** 2) / (T ** 2)) * (1 - D * (1 - (T0 ** 2) / (T ** 2)))
    return P


def temp2freq_old(T, out="freq"):
    """
    temperature > frequence du capteur de temperature (output == "freq")
    OU  coef U (output == "U") OU  periode du capteur U (output == "tau")
    """

    #### Handle T as an interable (array, list....)
    if utils.is_iterable(T):
        return np.array([temp2freq(t, out) for t in T])
    #### Handle T as a scalar
    else:
        U0 = 5.766353
        Y1 = -4025.183
        Y2 = -11970.76
        Y3 = 0.0
        U = np.roots([Y3, Y2, Y1, -T]) + U0
        # X = (1/f) * 10**6 # (1/30000) * 10**6
        # U = X - U0
        # Temp = Y1 * U + Y2 * U**2 + Y3 * U**3
        if out == "freq":
            return [
                10 ** 6 / u for u in U if ((10 ** 6 / u < 176000) and (10 ** 6 / u > 168000))
            ][0]
        elif out == "tau":
            return [
                10 ** 6 / u for u in U if ((10 ** 6 / u < 176000) and (10 ** 6 / u > 168000))
            ][0] ** -1
        elif out == "U":
            return [
                u - U0 for u in U if ((u > 10 ** 6 / 176000) and (u < 10 ** 6 / 168000))
            ][0]
        else:
            return None


def temp2freq_legacy(T, U0, Y, output="freq"):
    # U0 =  5.766353
    # Y1 = -4025.183
    # Y2 = -11970.76
    # Y3 =  0.
    U = np.roots([Y[2], Y[1], Y[0], -T]) + U0
    # X = (1/f) * 10**6 # (1/30000) * 10**6
    # U = X - U0
    # Temp = Y1 * U + Y2 * U**2 + Y3 * U**3
    if output == "freq":
        return [
            10 ** 6 / u for u in U if ((10 ** 6 / u < 176000) and (10 ** 6 / u > 168000))
        ][0]

    else:
        return [u - U0 for u in U if ((u > 10 ** 6 / 176000) and (u < 10 ** 6 / 168000))][0]
