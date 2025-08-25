#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.reffram contains functions for low-level
geometry operations.

it can be imported directly with:
from geodezyx import reffram

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

########## BEGIN IMPORT ##########
#### External modules
import io

#### Import the logger
import logging

import numpy as np
import pandas as pd
import scipy

#### geodeZYX modules
from geodezyx import conv
from geodezyx import stats
from geodezyx import utils

log = logging.getLogger("geodezyx")

##########  END IMPORT  ##########


#  _    _ _       _       _                    _    _____                _      _   _        ______   _
# | |  | (_)     | |     | |                  | |  / ____|              | |    | | (_)      |  ____| | |
# | |__| |_  __ _| |__   | |     _____   _____| | | |  __  ___  ___   __| | ___| |_ _  ___  | |__ ___| |_ ___
# |  __  | |/ _` | '_ \  | |    / _ \ \ / / _ \ | | | |_ |/ _ \/ _ \ / _` |/ _ \ __| |/ __| |  __/ __| __/ __|
# | |  | | | (_| | | | | | |___|  __/\ v /  __/ | | |__| |  __/ (_) | (_| |  __/ |_| | (__  | | | (__| |_\__ \
# |_|  |_|_|\__, |_| |_| |______\___| \_/ \___|_|  \_____|\___|\___/ \__,_|\___|\__|_|\___| |_|  \___|\__|___/
#            __/ |
#           |___/

### High level geodetic function


def itrf_speed_calc(x0, y0, z0, t0, vx, vy, vz, t):
    """
    Translate a position to a given epoch using point velocity

    Parameters
    ----------
    x0,y0,z0 : float
        coordinates at the reference epoch (m).
    t0 : float
        reference epoch (decimal year).
    vx,vy,vz : float
        speed of the point (m/yr).
    t : float
        output epoch.

    Returns
    -------
    xout,yout,zout : floats
        coordinates of the point @ the ref. epoch (m).
    """

    xout = x0 + vx * (t - t0)
    yout = y0 + vy * (t - t0)
    zout = z0 + vz * (t - t0)

    return xout, yout, zout


def itrf_psd_fundamuntal_formula(t, A_l, t_l, tau_l, A_e, t_e, tau_e):
    """
    Get the Post Seismic Deformation of a station

    Parameters
    ----------

    t : float
        epoch in decimal years.
    A_l : float
        Amplitude of the logarithmic term.
    t_l : float
        Earthquake time(date) corresponding to logarithmic term
    tau_l : float
        Relaxation time of the logarithmic term.
    A_e : float
        Amplitude of the exponential term.
    t_e : float
        Earthquake time(date) corresponding to the exponential term.
    tau_e : float
        Relaxation time of the exponential term..

    Returns
    -------
    dL : TYPE
        the total sum of PSD corrections.

    References
    ----------
        http://itrf.ensg.ign.fr/ITRF_solutions/2014/doc/ITRF2014-PSD-model-eqs-IGN.pdf
    """

    dL = A_l * np.log(1 + (t - t_l) / tau_l) + A_e * (1 + (t - t_e) / tau_e)

    return dL


def calc_pos_speed_itrf(x0, y0, z0, t0, vx, vy, vz, t):
    """
    just a wrapper of itrf_speed_calc
    for legacy reasons
    """
    return itrf_speed_calc(x0, y0, z0, t0, vx, vy, vz, t)


def itrf_helmert_get_parameters(
    TRF_Input_name, TRF_Ext_name, verbose=False, convert=True
):
    """
    Get the Helmert parameters for a transformation I/ETRFxx <=> I/ETRFxx

    Parameters
    ----------
    TRF_Input_name : str
        Name of the initial Reference Frame.
        Can be ITRFyyyy, ITRFyy, ETRFyyyy, ETRFyy
        where yyyy or yy is the TRF realization
    TRF_Ext_name : str
        Name of the wished Reference Frame.
        Can be ITRFyyyy, ITRFyy, ETRFyyyy, ETRFyy
        where yyyy or yy is the TRF realization
    verbose : bool, optional
        Print the parameters. The default is True.
    convert : bool, optional
        Gives directly useful values in the good units
        for a proper transformation

    Returns
    -------
    T : 3-Array
        Translation.
    Tdot : 3-Array
        Translation rate.
    D : float
        Scale factor.
    Ddot : float
        Scale factor rate.
    R : 3-Array
        Rotation.
    Rdot : 3-Array
        Rotation rte.
    epoch_ref_Xe : float
        Reference epoch of TRF_Ext_name.

    Notes
    -----
    Based on the values of
    http://etrs89.ensg.ign.fr/pub/EUREF-TN-1.pdf

    Warning: a validation unsing the EUREF converter is highly recomended
    https://www.epncb.oma.be/_productsservices/coord_trans/


    """

    ### Get short name if necessary
    if len(TRF_Input_name) == 8:
        TRF_Input_name = TRF_Input_name[:4] + TRF_Input_name[6:]
    if len(TRF_Ext_name) == 8:
        TRF_Ext_name = TRF_Ext_name[:4] + TRF_Ext_name[6:]

    COEFFS = """TRF_I,TRF_E,TRF_E_epoch,T1,T2,T3,D,R1,R2,R3,Tdot1,Tdot2,Tdot3,Ddot,Rdot1,Rdot2,Rdot3
ITRF14,ETRF14,1989.0,0.0,0.0,0.0,0.00,0.000,0.000,0.000,0.0,0.0,0.0,0.00,0.085,0.531,-0.770
ITRF05,ETRF05,1989.0,56.0,48.0,-37.0,0.00,0.000,0.000,0.000,0.0,0.0,0.0,0.00,0.054,0.518,-0.781
ITRF00,ETRF00,1989.0,54.0,51.0,-48.0,0.00,0.000,0.000,0.000,0.0,0.0,0.0,0.00,0.081,0.490,-0.792
ITRF97,ETRF97,1989.0,41.0,41.0,-49.0,0.00,0.000,0.000,0.000,0.0,0.0,0.0,0.00,0.200,0.500,-0.650
ITRF96,ETRF96,1989.0,41.0,41.0,-49.0,0.00,0.000,0.000,0.000,0.0,0.0,0.0,0.00,0.200,0.500,-0.650
ITRF94,ETRF94,1989.0,41.0,41.0,-49.0,0.00,0.000,0.000,0.000,0.0,0.0,0.0,0.00,0.200,0.500,-0.650
ITRF93,ETRF93,1989.0,19.0,53.0,-21.0,0.00,0.000,0.000,0.000,0.0,0.0,0.0,0.00,0.320,0.780,-0.670
ITRF92,ETRF92,1989.0,38.0,40.0,-37.0,0.00,0.000,0.000,0.000,0.0,0.0,0.0,0.00,0.210,0.520,-0.680
ITRF91,ETRF91,1989.0,21.0,25.0,-37.0,0.00,0.000,0.000,0.000,0.0,0.0,0.0,0.00,0.210,0.520,-0.680
ITRF90,ETRF90,1989.0,19.0,28.0,-23.0,0.00,0.000,0.000,0.000,0.0,0.0,0.0,0.00,0.110,0.570,-0.710
ITRF89,ETRF89,1989.0,0.0,0.0,0.0,0.00,0.000,0.000,0.000,0.0,0.0,0.0,0.00,0.110,0.570,-0.710
ITRF14,ETRF14,2010.0,0.0,0.0,0.0,0.00,1.785,11.151,-16.170,0.0,0.0,0.0,0.00,0.085,0.531,-0.770
ITRF08,ETRF14,2010.0,-1.6,-1.9,-2.4,0.02,1.785,11.151,-16.170,0.0,0.0,0.1,-0.03,0.085,0.531,-0.770
ITRF05,ETRF14,2010.0,-2.6,-1.0,2.3,-0.92,1.785,11.151,-16.170,-0.3,0.0,0.1,-0.03,0.085,0.531,-0.770
ITRF00,ETRF14,2010.0,-0.7,-1.2,26.1,-2.12,1.785,11.151,-16.170,-0.1,-0.1,1.9,-0.11,0.085,0.531,-0.770
ITRF97,ETRF14,2010.0,-7.4,0.5,62.8,-3.80,1.785,11.151,-16.430,-0.1,0.5,3.3,-0.12,0.085,0.531,-0.790
ITRF96,ETRF14,2010.0,-7.4,0.5,62.8,-3.80,1.785,11.151,-16.430,-0.1,0.5,3.3,-0.12,0.085,0.531,-0.790
ITRF94,ETRF14,2010.0,-7.4,0.5,62.8,-3.80,1.785,11.151,-16.430,-0.1,0.5,3.3,-0.12,0.085,0.531,-0.790
ITRF93,ETRF14,2010.0,50.4,-3.3,60.2,-4.29,4.595,14.531,-16.570,2.8,0.1,2.5,-0.12,0.195,0.721,-0.840
ITRF92,ETRF14,2010.0,-15.4,-1.5,70.8,-3.09,1.785,11.151,-16.430,-0.1,0.5,3.3,-0.12,0.085,0.531,-0.790
ITRF91,ETRF14,2010.0,-27.4,-15.5,76.8,-4.49,1.785,11.151,-16.430,-0.1,0.5,3.3,-0.12,0.085,0.531,-0.790
ITRF90,ETRF14,2010.0,-25.4,-11.5,92.8,-4.79,1.785,11.151,-16.430,-0.1,0.5,3.3,-0.12,0.085,0.531,-0.790
ITRF89,ETRF14,2010.0,-30.4,-35.5,130.8,-8.19,1.785,11.151,-16.430,-0.1,0.5,3.3,-0.12,0.085,0.531,-0.790
ITRF14,ETRF00,2010.0,54.7,52.2,-74.1,2.12,1.701,10.290,-16.632,0.1,0.1,-1.9,0.11,0.081,0.490,-0.792
ITRF08,ETRF00,2010.0,53.1,50.3,-76.5,2.14,1.701,10.290,-16.632,0.1,0.1,-1.8,0.08,0.081,0.490,-0.792
ITRF05,ETRF00,2010.0,52.1,51.2,-71.8,1.20,1.701,10.290,-16.632,-0.2,0.1,-1.8,0.08,0.081,0.490,-0.792
ITRF00,ETRF00,2010.0,54.0,51.0,-48.0,0.00,1.701,10.290,-16.632,0.0,0.0,0.0,0.00,0.081,0.490,-0.792
ITRF97,ETRF00,2010.0,47.3,52.7,-11.3,-1.68,1.701,10.290,-16.892,0.0,0.6,1.4,-0.01,0.081,0.490,-0.812
ITRF96,ETRF00,2010.0,47.3,52.7,-11.3,-1.68,1.701,10.290,-16.892,0.0,0.6,1.4,-0.01,0.081,0.490,-0.812
ITRF94,ETRF00,2010.0,47.3,52.7,-11.3,-1.68,1.701,10.290,-16.892,0.0,0.6,1.4,-0.01,0.081,0.490,-0.812
ITRF93,ETRF00,2010.0,105.1,48.9,-13.9,-2.17,4.511,13.670,-17.032,2.9,0.2,0.6,-0.01,0.191,0.680,-0.862
ITRF92,ETRF00,2010.0,39.3,50.7,-3.3,-0.97,1.701,10.290,-16.892,0.0,0.6,1.4,-0.01,0.081,0.490,-0.812
ITRF91,ETRF00,2010.0,27.3,36.7,2.7,-2.37,1.701,10.290,-16.892,0.0,0.6,1.4,-0.01,0.081,0.490,-0.812
ITRF90,ETRF00,2010.0,29.3,40.7,18.7,-2.67,1.701,10.290,-16.892,0.0,0.6,1.4,-0.01,0.081,0.490,-0.812
ITRF89,ETRF00,2010.0,24.3,16.7,56.7,-6.07,1.701,10.290,-16.892,0.0,0.6,1.4,-0.01,0.081,0.490,-0.812
ITRF14,ITRF08,2010.0,1.6,1.9,2.4,-0.02,0.00,0.00,0.00,0.0,0.0,-0.1,0.03,0.00,0.00,0.00
ITRF14,ITRF05,2010.0,2.6,1.0,-2.3,0.92,0.00,0.00,0.00,0.3,0.0,-0.1,0.03,0.00,0.00,0.00
ITRF14,ITRF00,2010.0,0.7,1.2,-26.1,2.12,0.00,0.00,0.00,0.1,0.1,-1.9,0.11,0.00,0.00,0.00
ITRF14,ITRF97,2010.0,7.4,-0.5,-62.8,3.80,0.00,0.00,0.26,0.1,-0.5,-3.3,0.12,0.00,0.00,0.02
ITRF14,ITRF96,2010.0,7.4,-0.5,-62.8,3.80,0.00,0.00,0.26,0.1,-0.5,-3.3,0.12,0.00,0.00,0.02
ITRF14,ITRF94,2010.0,7.4,-0.5,-62.8,3.80,0.00,0.00,0.26,0.1,-0.5,-3.3,0.12,0.00,0.00,0.02
ITRF14,ITRF93,2010.0,-50.4,3.3,-60.2,4.29,-2.81,-3.38,0.40,-2.8,-0.1,-2.5,0.12,-0.11,-0.19,0.07
ITRF14,ITRF92,2010.0,15.4,1.5,-70.8,3.09,0.00,0.00,0.26,0.1,-0.5,-3.3,0.12,0.00,0.00,0.02
ITRF14,ITRF91,2010.0,27.4,15.5,-76.8,4.49,0.00,0.00,0.26,0.1,-0.5,-3.3,0.12,0.00,0.00,0.02
ITRF14,ITRF90,2010.0,25.4,11.5,-92.8,4.79,0.00,0.00,0.26,0.1,-0.5,-3.3,0.12,0.00,0.00,0.02
ITRF14,ITRF89,2010.0,30.4,35.5,-130.8,8.19,0.00,0.00,0.26,0.1,-0.5,-3.3,0.12,0.00,0.00,0.02
ITRF14,ITRF88,2010.0,25.4,-0.5,-154.8,11.29,0.10,0.00,0.26,0.1,-0.5,-3.3,0.12,0.00,0.00,0.02
"""

    COEFFS = io.StringIO(COEFFS)
    # p="/home/psakicki/Downloads/coeff_helmert2.csv"
    DF = pd.read_csv(COEFFS)

    DF2 = DF[(DF["TRF_I"] == TRF_Input_name) & (DF["TRF_E"] == TRF_Ext_name)]
    inver = 1.0
    inver_bool = False

    if len(DF2) == 0:
        DF2 = DF[(DF["TRF_E"] == TRF_Input_name) & (DF["TRF_I"] == TRF_Ext_name)]
        inver = -1.0
        inver_bool = True

    DF3 = DF2.iloc[0]

    if convert:
        convT = 0.001
        convD = 10**-9
        convR = 4.84813681109536e-09
    else:
        convT = 0.0
        convD = 0.0
        convR = 0.0

    T = inver * convT * DF3[["T1", "T2", "T3"]].values
    Tdot = inver * convT * DF3[["Tdot1", "Tdot2", "Tdot3"]].values

    D = inver * convD * DF3["D"]
    Ddot = inver * convD * DF3["Ddot"]

    R = inver * convR * DF3[["R1", "R2", "R3"]].values
    Rdot = inver * convR * DF3[["Rdot1", "Rdot2", "Rdot3"]].values

    epoch_ref_Xe = DF3["TRF_E_epoch"]

    if verbose:
        log.info("%s => %s", TRF_Input_name, TRF_Ext_name)
        log.info("   T: %s", T, Tdot)
        log.info("   D: %s", D, Ddot)
        log.info("   R: %s", R, Rdot)
        log.info("t Xe: %s", epoch_ref_Xe)
        log.info(" Inv: %s", inver_bool)
        log.info("DataFrame:")
        log.info(DF2)

    return T, Tdot, D, Ddot, R, Rdot, epoch_ref_Xe


def itrf_helmert_trans(Xi, epoch_Xi, T, Tdot, D, Ddot, R, Rdot, epoch_ref_Xe):
    """
    Do the Helmert transformation I/ETRFxx <=> I/ETRFxx

    The seven Helmert Parameters are generated with ``itrf_helmert_get_parameters``

    Parameters
    ----------
    Xi : 3-Array, Nx3-Array OR 6-Array, Nx6-Array
        Points in the initial Reference Frame.
        3 columns if positions only
        6 columns if positions+velocities
    epoch_Xi : float
        epoch of the initial Reference Frame points (decimal year).
    T : 3-Array
        Translation parameter.
    Tdot : 3-Array
        Translation rate parameter.
    D : float
        Scale factor parameter.
    Ddot : float
        Scale factor rate parameter.
    R : 3-Array
        Rotation parameter.
    Rdot : 3-Array
        Rotation rate parameter.
    epoch_ref_Xe : float
        Reference epoch of the destination TRF.

    Returns
    -------
    Xe : 3-Array or Nx3-Array
        Points in the destination Reference Frame.



    Notes
    -----

    ** Warning **
    This function does not the velocity shift from one epoch to another \n
    The output coordinates will be provided at the same epoch as the input ones \n
    Use ``itrf_speed_calc`` function to perform this potential step
    before (with the initial velocities) or after (with the transformed velocities)
    calling the present function.


    Based on the theory and values of *EUREF Technical Note 1:
    Relationship and Transformation between the International
    and the European Terrestrial Reference Systems, Z. Altamimi, 2018*
    http://etrs89.ensg.ign.fr/pub/EUREF-TN-1.pdf

    We recommend to confirm the values with the official EUREF converter
    http://www.epncb.oma.be/_productsservices/coord_trans/index.php


    ** Notes for the French geodetic system **
    By definition, since the 2021/01/04
    RGF93(v2b) = ETRF2000@2019.0

    References
    ----------
    https://geodesie.ign.fr/index.php?page=rgf93 \n
    https://geodesie.ign.fr/contenu/fichiers/RGF93v2b-RAF18b.pdf \n
    https://geodesie.ign.fr/contenu/fichiers/rgf93v2b_information_cnig.pdf

    """
    ## prelimiary warning
    if utils.is_iterable(epoch_Xi):
        log.warning("epoch_Xi is an iterable !!!")
        log.info("The function works correctly with a single epoch only !!!")
        log.info("if several initial Reference Frame epoch")
        log.info("use a loop outside the function")

    if not type(Xi) is np.array:
        Xi = np.array(Xi)

    # manage the case where only one point is given
    if len(np.shape(Xi)) == 1:
        Xi = Xi[..., np.newaxis].T

    ##### Here we manage if there is velocities or not
    if Xi.shape[1] == 6:
        velocity_mode = True
    else:
        velocity_mode = False

    #### Trick if we are in velocity mode, we recall the fct
    #### with velocity_mode=False to have the positions
    if velocity_mode:
        Xe_pos = itrf_helmert_trans(
            Xi[:, :3], epoch_Xi, T, Tdot, D, Ddot, R, Rdot, epoch_ref_Xe
        )

    #### Translation
    if velocity_mode:
        Topera = Tdot
    else:
        Topera = T + Tdot * (epoch_Xi - epoch_ref_Xe)

    #### Scale Factor
    D_mat = np.eye(3).dot(D)
    Ddot_mat = np.eye(3).dot(Ddot)
    if velocity_mode:
        D_mat_opera = Ddot_mat
    else:
        D_mat_opera = Ddot_mat * (epoch_Xi - epoch_ref_Xe) + D_mat

    ### Rotation
    r_1, r_2, r_3 = R
    r_dot_1, r_dot_2, r_dot_3 = Rdot

    R_mat = np.array([[0.0, -r_3, r_2], [r_3, 0.0, -r_1], [-r_2, r_1, 0.0]])

    Rdot_mat = np.array(
        [[0.0, -r_dot_3, r_dot_2], [r_dot_3, 0.0, -r_dot_1], [-r_dot_2, r_dot_1, 0.0]]
    )

    if velocity_mode:
        R_mat_opera = Rdot_mat
    else:
        R_mat_opera = Rdot_mat * (epoch_Xi - epoch_ref_Xe) + R_mat

    #### Computation
    Xe_stk = []
    for xi in Xi:
        if velocity_mode:
            xe = xi[3:] + Topera + D_mat_opera.dot(xi[:3]) + R_mat_opera.dot(xi[:3])
        else:
            xe = xi[:3] + Topera + D_mat_opera.dot(xi[:3]) + R_mat_opera.dot(xi[:3])

        Xe_stk.append(xe)

    #### Final Aestetics
    Xe = np.array(Xe_stk)
    Xe = np.squeeze(Xe)
    Xe = Xe.astype(np.float64)

    ### final concat of the position and the velocities
    if velocity_mode:
        if len(np.shape(Xe)) == 1:
            Xe = np.hstack((Xe_pos, Xe))
        else:
            Xe = np.column_stack((Xe_pos, Xe))

    return Xe


def seasonal_geocenter_motion(year_frac_inp, reverse=False):
    """
    Calculate the ITRF2020 seasonal geocenter motion along X, Y, and Z axes.

    The sign and phase conventions used are such that the seasonal
    motion of the Earth's center of figure (CF) with respect to
    the Earth's center of mass (CM).

    Parameters
    ----------
    year_frac_inp : float
        Time in fractional years.
    reverse : bool, optional
        If True, reverse the phase of 180deg.
        and then describes the motion of the CM w.r.t. CF.
        Default is False.

    Returns
    -------
    cf_x : float
        Seasonal geocenter motion in meters along the X axis.
    cf_y : float
        Seasonal geocenter motion in meters along the Y axis.
    cf_z : float
        Seasonal geocenter motion in meters along the Z axis.

    Source
    ------
    https://itrf.ign.fr/ftp/pub/itrf/itrf2020/ITRF2020-geocenter-motion.dat

    Altamimi, Z., Rebischung, p., Collilieux, X. et al. ITRF2020:
    an augmented reference frame refining the modeling of nonlinear station motions.
    J Geod 97, 47 (2023). https://doi.org/10.1007/s00190-023-01738-w

    https://link.springer.com/article/10.1007/s00190-023-01738-w
    """

    year_frac = (year_frac_inp - np.floor(year_frac_inp))

    # Annual and semi-annual amplitudes and phases for X, Y, Z axes
    a1_x, phi1_x = 1.23, -123.2
    a1_y, phi1_y = 3.48, 152.9
    a1_z, phi1_z = 2.76, -139.5

    a2_x, phi2_x = 0.49, 107.2
    a2_y, phi2_y = 0.22, 1.6
    a2_z, phi2_z = 1.19, 30.5

    if not reverse:
        phi_delta = 0
    else:
        # Altamimi et al. (2023), end of sec. 5.3
        phi_delta = np.pi

    # Convert phases from degrees to radians
    phi1_x_rad = np.deg2rad(phi1_x) + phi_delta
    phi1_y_rad = np.deg2rad(phi1_y) + phi_delta
    phi1_z_rad = np.deg2rad(phi1_z) + phi_delta

    phi2_x_rad = np.deg2rad(phi2_x) + phi_delta
    phi2_y_rad = np.deg2rad(phi2_y) + phi_delta
    phi2_z_rad = np.deg2rad(phi2_z) + phi_delta

    # Calculate the CF motion along each axis
    cf_x = a1_x * np.cos(2 * np.pi * year_frac - phi1_x_rad) + a2_x * np.cos(
        4 * np.pi * year_frac - phi2_x_rad
    )
    cf_y = a1_y * np.cos(2 * np.pi * year_frac - phi1_y_rad) + a2_y * np.cos(
        4 * np.pi * year_frac - phi2_y_rad
    )
    cf_z = a1_z * np.cos(2 * np.pi * year_frac - phi1_z_rad) + a2_z * np.cos(
        4 * np.pi * year_frac - phi2_z_rad
    )

    cf_x, cf_y, cf_z = cf_x * 10**-3, cf_y * 10**-3, cf_z * 10**-3

    return cf_x, cf_y, cf_z


#  _    _      _                     _     _                        __                           _   _
# | |  | |    | |                   | |   | |                      / _|                         | | (_)
# | |__| | ___| |_ __ ___   ___ _ __| |_  | |_ _ __ __ _ _ __  ___| |_ ___  _ __ _ __ ___   __ _| |_ _  ___  _ __
# |  __  |/ _ \ | '_ ` _ \ / _ \ '__| __| | __| '__/ _` | '_ \/ __|  _/ _ \| '__| '_ ` _ \ / _` | __| |/ _ \| '_ \
# | |  | |  __/ | | | | | |  __/ |  | |_  | |_| | | (_| | | | \__ \ || (_) | |  | | | | | | (_| | |_| | (_) | | | |
# |_|  |_|\___|_|_| |_| |_|\___|_|   \__|  \__|_|  \__,_|_| |_|___/_| \___/|_|  |_| |_| |_|\__,_|\__|_|\___/|_| |_|
#


def _helmert_trans_estim_matrixs_maker(X1, X2):
    """
    internal function for helmert_trans_estim
    """
    x1, y1, z1 = X1
    x2, y2, z2 = X2

    block_1 = np.eye(3)

    block_2 = np.array([[0.0, -z1, y1, x1], [z1, 0.0, -x1, y1], [-y1, x1, 0.0, z1]])

    l = X2 - X1
    A = np.hstack((block_1, block_2))

    return l, A


# def helmert_trans_frontend(X1list , X2list, Weights=[]):
#     HParam , A , l = reffram.helmert_trans_estim(X1list,X2list)
#     X2_trans_out = reffram.helmert_trans_apply(X1list,HParam)


def helmert_trans_estim(X1list, X2list, Weights=[]):
    """
    estimates 7 parameters of a 3D Helmert transformation between a set of points
    X1 and a set of points X2 (compute transformation X1 => X2)

    Parameters
    ----------

    X1list & X2list : list or np.array
        Input point set
        list of N (x,y,z) points, or an (N,3)-shaped numpy array

    Weights : list of N Weights,
        or an numpy array of shape (N,)

    Returns
    -------

    HParam :
        7 Helmert params. : x,y,z translations, x,y,z rotations, scale
        translations are given in meters, rotations in arcsec,
        scale in unitless ratio
    A :
        Design matrix
    l :
        Differences X2 - X1 (before transformation !!!)

    References
    ----------
    https://elib.uni-stuttgart.de/bitstream/11682/9661/1/BscThesis_GaoYueqing.pdf
    """

    l_stk = []
    A_stk = []
    Bool_stk = []

    log.info("we have %s points", len(X1list))

    for X1, X2 in zip(X1list, X2list):

        if np.sum(np.isnan(X1)) or np.sum(np.isnan(X2)):
            log.warning("one point component is nan, skipping")
            Bool_stk.append(False)
            continue

        lmono, Amono = _helmert_trans_estim_matrixs_maker(X1, X2)

        l_stk.append(lmono)
        A_stk.append(Amono)
        Bool_stk.append(True)

    A = np.vstack(A_stk)
    l = np.hstack(l_stk)
    Bool = np.array(Bool_stk)

    if len(Weights) == 0:
        W = np.eye(len(l))
    else:
        Weights_corr = np.array(Weights)
        Weights_corr = Weights_corr[Bool]
        W = np.repeat(Weights_corr, 3) * np.eye(len(l))

    N = (A.T).dot(W).dot(A)
    AtWB = (A.T).dot(W).dot(l)

    HParam = np.linalg.inv(N).dot(AtWB)

    return HParam, A, l


def helmert_trans_apply(Xin, SevenParam_in, legacy_mode=False):
    """
    Apply an Helmert transformation (7-parameters)
    to a set of points

    Parameters
    ----------
    Xin : list or np.array
        input point set
        list of N (x,y,z) points, or an (N,3)-shaped numpy array.

    SevenParam_in : 7 element list or array
        7 Helmert params. : x,y,z translations, x,y,z rotations, scale.

    legacy_mode : bool, optional
        Use a non-optimized and slow computation approach (but same result).
        This option should be removed in the Future.
        The default is False.

    Returns
    -------
    Xout : list/array of N (x,y,z) points
        output transformed points. Same type as the input.

    """

    tx, ty, tz, rx, ry, rz, scal = SevenParam_in

    R = np.array([[1.0, rz, -ry], [-rz, 1.0, rx], [ry, -rx, 1.0]])

    T = np.array([tx, ty, tz])
    S = 1.0 + scal

    typ = utils.get_type_smart(Xin)

    #### Apply the transformation here
    if legacy_mode:  #### SLOW !!!
        Xout = []
        for X1 in Xin:
            X2 = S * np.dot(R, X1) + T
            Xout.append(X2)
    else:  ##### 100x faster with the Einstein sum
        Xout = S * np.einsum("ij,kj->ki", R, Xin) + np.tile(T, (len(Xin), 1))

    Xout = typ(Xout)

    return Xout


def helmert_trans_estim_minimisation(
    X1in,
    X2in,
    HParam_apri=np.zeros(7),
    L1norm=True,
    tol=10**-9,
    full_output=False,
    method="Powell",
):
    """
    estimates 7 parameters of a 3D Helmert transformation between a set of points
    X1 and a set of points X2 (compute transformation X1 => X2)
    using a Minimization approach (and not a Least Square inversion)

    Parameters
    ----------

    X1in & X2in :  list or np.array
        Input point set
        list of N (x,y,z) points, or an (N,3)-shaped numpy array

    HParam_apri : list
        list of 7 values,
        The Apriori for the Helmert parameter

    L1norm : bool
        Use the L1-norm as a criteria, use the quadratic sum instead if False

    tol : float
        tolerence for the convergence

    full_output : bool
        return only the result if True, return the scipy optimize result if False

    method : str, optional
        minimization method.
        see scipy.optimize.minimize for details
        The default is "Powell".

    Returns
    -------
    Res :
        7 Helmert params. : x,y,z translations, x,y,z rotations, scale
        translations are given in meters, rotations in arcsec,
        scale in unitless ratio
    """

    def minimiz_helmert_fct(HParam_mini_in, X1in, X2in, L1norm_mini=L1norm):
        """
        This fct is the input for the scipy optimization fct
        """
        HParam_mini_wk = HParam_mini_in.copy()
        X12wrk = helmert_trans_apply(X1in, HParam_mini_wk)

        if not L1norm_mini:  # return a L2 norm
            SUM = np.sum((np.sum(np.power(X2in - X12wrk, 2), axis=1)))
        else:  # Return L1 norm
            SUM = np.sum(np.sum(np.abs(X2in - X12wrk), axis=1))
        return SUM

    RES = scipy.optimize.minimize(
        minimiz_helmert_fct,
        HParam_apri,
        (X1in, X2in, L1norm),
        method=method,
        tol=tol,
        options={"maxiter": 1000, "xtol": tol, "ftol": tol},
    )

    if RES.status != 0:
        log.warning("something went wrong (status != 0)")
        log.warning("here is the scipy.optimize.minimize message")
        log.warning(" > " + RES.message)

    if not full_output:
        return RES.x
    else:
        return RES


def helmert_trans_estim_minimisation_scalar(
    X1, X2, HParam_opti_apriori, L1norm=True, itera=2
):
    """
    Estimates the Helmert parameters but based on a minimisation approach
    between X1 and X2 (as suggested in the IGS combination software)

    NOT STABLE AVOID THE USE
    (it does the optimization of the 1st parameter,
    and then the optimization of the 2nd one etc, etc...)

    Parameters
    ----------
    X1 & X2 : list of N (x,y,z) points, or an (N,3)-shaped numpy array
        Input point sets
    HParam_opti_apriori : list of 7 values
        The Apriori for the Helmert parameter
    L1norm : bool
        Use the L1-norm as a criteria, use the quadratic sum instead if False
    itera : int, optional
        number of iterations. The default is 2.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    log.warning("unstable approach, avoid !!!!!")

    def minimiz_helmert_fct_scalar(
        hparam_mono_in, hparam_mono_id, HParam_mini_in, X1in, X2in
    ):
        """
        This fct is the input for the scipy optimization fct
        """
        HParam_mini_wk = HParam_mini_in.copy()
        HParam_mini_wk[hparam_mono_id] = hparam_mono_in
        X12wrk = helmert_trans_apply(X1in, HParam_mini_wk)

        if not L1norm:  # return a L2 norm
            return np.sum(np.sqrt(np.sum(np.power(X2in - X12wrk, 2), axis=1)))
        else:  # Return L1 norm
            return np.sum(np.sum(np.abs(X2in - X12wrk), axis=1))

    HParam_opti_wrk = HParam_opti_apriori.copy()

    for j in range(itera):  # iter iterations
        for i in range(7):  # 7 Helmert Parameters:
            RES = scipy.optimize.minimize_scalar(
                minimiz_helmert_fct_scalar,
                args=(i, HParam_opti_wrk, X1, X2),
                tol=10**-20,
            )
            HParam_opti_wrk[i] = RES.x

        if not j:
            HParam_opti_prev = HParam_opti_wrk.copy()
        else:
            log.info(
                "Helmert Param minimisation iter %s %s",
                j + 1,
                HParam_opti_wrk - HParam_opti_prev,
            )
            HParam_opti_prev = HParam_opti_wrk.copy()

    return HParam_opti_wrk


#### ASTRONOMY FUNCTION
def semi_major_axis_from_mean_motion(n):
    """
    source : https://space.stackexchange.com/questions/18289/how-to-get-semi-major-axis-from-tle
    """
    mu = 3.9860044189 * 10**14
    a = (mu ** (1.0 / 3.0)) / ((2 * n * np.pi / 86400) ** (2.0 / 3.0))
    return a


#  _                     _                    _    _____                           _        _        ______                _   _
# | |                   | |                  | |  / ____|                         | |      (_)      |  ____|              | | (_)
# | |     _____      __ | |     _____   _____| | | |  __  ___  ___  _ __ ___   ___| |_ _ __ _  ___  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
# | |    / _ \ \ /\ / / | |    / _ \ \ / / _ \ | | | |_ |/ _ \/ _ \| '_ ` _ \ / _ \ __| '__| |/ __| |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# | |___| (_) \ v  v /  | |___|  __/\ v /  __/ | | |__| |  __/ (_) | | | | | |  __/ |_| |  | | (__  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
# |______\___/ \_/\_/   |______\___| \_/ \___|_|  \_____|\___|\___/|_| |_| |_|\___|\__|_|  |_|\___| |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/


def BL_from_points(listpointin):
    """
    From a list of 2-D or 3-dD points, returns the a matrix with distance
    between each points

    Parameters
    ----------
    listpointin : list or numpy.array
        List of N 2D or 3D points [[x1,y1,z1] ... [xn , yn , zn]]

    Returns
    -------
    BL : numpy.array
        matrix with distances between each points

    """

    N = len(listpointin)
    BL = np.empty((N, N))

    for i, pt1 in enumerate(listpointin):
        for j, pt2 in enumerate(listpointin):

            if i == j:
                BL[i, j] = 0
            else:
                BL[i, j] = np.linalg.norm(pt1 - pt2)

    return BL


def rotmat2(theta, angtype="deg"):

    if angtype == "deg":
        theta = np.deg2rad(theta)

    rotmat = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])

    return rotmat


def rotmat3(
    alpha, beta, gamma, xyzreftuple=([1, 0, 0], [0, 1, 0], [0, 0, 1]), angtype="deg"
):

    xaxis, yaxis, zaxis = xyzreftuple

    if angtype == "deg":
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)

    Rx = trans.rotation_matrix(alpha, xaxis)
    Ry = trans.rotation_matrix(beta, yaxis)
    Rz = trans.rotation_matrix(gamma, zaxis)
    R = trans.concatenate_matrices(Rz, Ry, Rx)[:3, :3]

    return R


def rotate_points(
    alphal,
    betal,
    gammal,
    pointlin,
    Rtype="R1",
    xyzreftuple=([1, 0, 0], [0, 1, 0], [0, 0, 1]),
    angtype="deg",
    fullout=False,
):
    """
    R1  = Rz(g) * Ry(b) * Rx(a)
         si les RPY sont donnés dans le NED
         alors les positions résultantes sont dans le NED

    R2  =  matrice RPY2ENU
        si les RPY sont donnés dans le NED
        alors les  résultantes sont DANS LE ENU
        pas besoin de rotation NED2ENU

        Grewal et al. 2007

    Entrée :
        Angles n = A
        liste de listes de p * [ points ]

    Sortie :
        liste de listes [ [ xA ] [ xA ] ... xP [ xA ] ]
    """

    xaxis, yaxis, zaxis = xyzreftuple

    if not utils.is_iterable(alphal):
        alphal = np.array([alphal])
        betal = np.array([betal])
        gammal = np.array([gammal])
        boolnotiterable = True
    else:
        boolnotiterable = False

    pointlout = []
    R_out = []

    for pt in pointlin:

        if not utils.is_iterable(pt) or len(pt) != 3:
            log.error("pts != 3 coords")
            return 0

        pointltmp = []

        for a, b, g in zip(alphal, betal, gammal):

            R1 = rotmat3(a, b, g, angtype=angtype, xyzreftuple=xyzreftuple)
            R2 = conv.c_rpy2enu(a, b, g, angtype=angtype)

            if Rtype == "R1":
                R = R1
            elif Rtype == "R2":
                R = R2
            R_out.append(R)

            pointltmp.append(np.dot(R, pt))

        pointlout.append(pointltmp)

        if boolnotiterable:
            pointlout = pointltmp

        pointlout = np.array(pointlout)

    if fullout:
        return pointlout, R_out
    else:
        return pointlout


def guess_seq_len(seq):
    # source
    # http://stackoverflow.com/questions/11385718/python-finding-repeating-sequence-in-list-of-integers
    guess = 1

    if len(set(seq)) == 1:
        return 1

    max_len = len(seq) / 2
    for x in range(2, max_len):
        if seq[0:x] == seq[x : 2 * x]:
            return x

    return guess


def wrapTo2Pi(lon):
    """
     wrapTo2Pi Wrap angle in radians to [0 2*pi]

    lambdaWrapped = wrapTo2Pi(LAMBDA) wraps angles in LAMBDA, in radians,
    to the interval [0 2*pi] such that zero maps to zero and 2*pi maps
    to 2*pi. (In general, positive multiples of 2*pi map to 2*pi and
    negative multiples of 2*pi map to zero.)

    See also wrapToPi, wrapTo180, wrapTo360.

    """
    lon = np.array(lon)
    positiv = lon > 0
    outlon = np.mod(lon, 2 * np.pi)
    outlon[np.logical_and(outlon == 0, positiv)] = 2 * np.pi
    return outlon


def wrapToPi(lon):
    """
    wrapToPi Wrap angle in radians to [-pi pi]

       lambdaWrapped = wrapToPi(LAMBDA) wraps angles in LAMBDA, in radians,
       to the interval [-pi pi] such that pi maps to pi and -pi maps to
       -pi.  (In general, odd, positive multiples of pi map to pi and odd,
       negative multiples of pi map to -pi.)

       See also wrapTo2Pi, wrapTo180, wrapTo360.

    """

    outlon = np.array(lon)
    q = np.logical_and((outlon < -np.pi), (np.pi < outlon))
    outlon[q] = wrapTo2Pi(outlon[q] + np.pi) - np.pi
    return outlon


def wrapTo180(lonin):
    """
    wrapTo180 Wrap angle in degrees to [-180 180]

    lonWrapped = wrapTo180(LON) wraps angles in LON, in degrees, to the
    interval [-180 180] such that 180 maps to 180 and -180 maps to -180.
    (In general, odd, positive multiples of 180 map to 180 and odd,
    negative multiples of 180 map to -180.)

    See also wrapTo360, wrapTo2Pi, wrapToPi.
    """
    lon = np.array(lonin)
    q = (lon < -180) and (180 < lon)
    lon[q] = wrapTo360(lon[q] + 180) - 180

    return lon


def wrapTo360(lonin):
    """
    wrapTo360 Wrap angle in degrees to [0 360]

    lonWrapped = wrapTo360(LON) wraps angles in LON, in degrees, to the
    interval [0 360] such that zero maps to zero and 360 maps to 360.
    (In general, positive multiples of 360 map to 360 and negative
    multiples of 360 map to zero.)

    See also wrapTo180, wrapToPi, wrapTo2Pi.
    """

    lon = np.array(lonin)

    positiveInput = lon > 0
    lon = np.mod(lon, 360)
    lon[(lon == 0) & positiveInput] = 360
    return lon


# Pas convaincu de son utilité
def unwrap180(anglist, angtype="deg"):

    if angtype == "deg":
        seuil = 360

    angout = []

    for a in anglist:
        if a > seuil / 2:
            a = a - seuil
        angout.append(a)

    return angout


def wrap360(anglist, angtype="deg"):

    angout = []

    if angtype == "deg":
        seuil = 360
    elif angtype == "rad":
        seuil = 2 * np.pi

    for a in anglist:
        if a < 0:
            a = a + seuil

        angout.append(a)

    return angout


class interp1d_ang:

    def __init__(self, T, A, angtype="deg", kind="linear", bounds_error=False):

        if angtype == "deg":
            A = np.deg2rad(A)

        self.A = A
        self.T = T
        self.C = np.cos(A)
        self.S = np.sin(A)

        self.CfT = scipy.interpolate.interp1d(
            T, self.C, kind=kind, bounds_error=bounds_error
        )
        self.SfT = scipy.interpolate.interp1d(
            T, self.S, kind=kind, bounds_error=bounds_error
        )

    def __call__(self, T, angtype="deg"):

        I = np.arctan2(self.SfT(T), self.CfT(T))
        I = wrap360(I, angtype="rad")

        if angtype == "deg":
            return np.rad2deg(I)
        else:
            return I


def group_consecutives(vals, step=1):
    """
    Return list of consecutive lists of numbers from vals (number list).
    """
    run = []
    result = [run]
    expect = None
    for v in vals:
        if (v == expect) or (expect is None):
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step

        result2 = []
        for r in result:
            if len(r) > 1:
                result2.append([r[0], r[-1]])
            else:
                result2.append(r)

    return result2


def randomwalk_normal(N=100, d=2, moy=0, sigma=1):
    """
    d = dimension
    """
    return np.cumsum(moy + np.random.randn(N, d) * sigma)


def randomwalk_uniform(N=100, d=2, bound=0.5):
    """
    d = dimension
    bound = contraint of the random walk
    """
    return np.cumsum(np.random.uniform(-bound, bound, (N, d)))


def circle_draw(xc, yc, R, N):
    theta = np.linspace(0, 2 * np.pi, N)
    X = np.cos(theta) * R + xc
    Y = np.sin(theta) * R + yc
    return X, Y


def random_walk_in_a_circle(
    x0,
    y0,
    xc,
    yc,
    R,
    N,
    step_size,
    param=1,
    polar=True,
    uniform_or_normal="n",
    rand_seed=-1,
):
    """
    generate random walk in a circle

    Parameters
    ----------
    x0 : float
        x start of the random walk.
    y0 : float
        y start of the random walk.
    xc : float
        x center of the circle.
    yc : float
        y center of the circle.
    R : float
        radius of the circle.
    N : int
        number of points.
    step_size : float
        size of the step.
    param : int, optional
        control the random. The default is 1.
    polar : TYPE, optional
        DESCRIPTION. The default is True.
    uniform_or_normal : 'u' or 'n', optional
        uniform of normal walk. The default is 'n'.
    rand_seed : int, optional
        control the random. The default is -1.

    Returns
    -------
    X : np.array
        X-coordinates of the random walk.
    Y : np.array
        Y-coordinates of the random walk.
    Xcircle : np.array
        X-coordinates of the circle.
    Ycircle : np.array
        Y-coordinates of the circle.

    Exemple
    -------
    for un in ('u','n'):
    for pol in range(2):
        X,Y , Xcircle , Ycircle = random_walk_in_a_circle(10,10,0,0,50,10000,polar = pol,uniform_or_normal=un)

        plt.figure()
        plt.plot(Xcircle,Ycircle)
        plt.plot(X,Y)
        plt.axis('equal')
        plt.suptitle(un + str(pol))

    """

    X = [x0]
    Y = [y0]

    if rand_seed > -1:
        RAND = np.random.RandomState(rand_seed)
    else:
        RAND = np.random.RandomState(np.random.randint(10**6))

    Xcircle, Ycircle = circle_draw(xc, yc, R, 500)

    for i in range(N - 1):
        D = R + 1
        iwhil = 0
        while D > R:
            iwhil += 1
            if iwhil > 500:
                log.warning("infinite loop in random_walk_in_a_circle ... %s")
            if polar:
                if uniform_or_normal == "u":
                    dalpha = RAND.uniform(-param, param) * 2 * np.pi
                else:
                    dalpha = RAND.normal(0, param) * 2 * np.pi
                drho = step_size
                dx = drho * np.cos(dalpha)
                dy = drho * np.sin(dalpha)
                # print dx,dy
            else:
                if uniform_or_normal == "u":
                    dx = np.random.uniform(-param, param)
                    dy = np.random.uniform(-param, param)
                else:
                    dx = np.random.normal(0, param)
                    dy = np.random.normal(0, param)

            xtemp = X[-1] + dx
            ytemp = Y[-1] + dy
            D = np.sqrt((xtemp - xc) ** 2 + (ytemp - yc) ** 2)

        X.append(xtemp)
        Y.append(ytemp)

    X = np.array(X)
    Y = np.array(Y)

    return X, Y, Xcircle, Ycircle


def randn_bool(N, true_ratio=0.5, RandGene=None):
    if RandGene is None:
        RandGene = np.random.RandomState()
    if type(RandGene) is int:
        RandGene = np.random.RandomState(RandGene)
    try:
        randlis = RandGene.uniform(size=N)
    except AttributeError:
        "ERR : AttributeError : RandGene  may be an int32/int64, but it's not an authentic int as required ..."
    boolout_lis = []
    for r in randlis:
        if r < true_ratio:
            boolout_lis.append(True)
        else:
            boolout_lis.append(False)

    return boolout_lis


def points_circle_border(
    Npts,
    r,
    r_sigma,
    az_type_normal=True,
    main_dir=3.14159,
    dir_range=3.14159,
    seed=None,
):
    if not seed:
        seed = np.random.randint(10000)

    S = np.random.RandomState(seed)

    if not az_type_normal:
        Az = S.rand(Npts) * 2 * np.pi
    else:
        Az = S.randn(Npts) * dir_range + main_dir

    R = np.array(Npts * [r]) - np.abs(S.randn(Npts) * r_sigma)

    X, Y = conv.polar2cartesian(R, Az, "rad")

    return X, Y


def estimated_autocorrelation(x):
    """
    http://stackoverflow.com/q/14297012/190597
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    """
    n = len(x)
    variance = x.var()
    x = x - x.mean()
    r = np.correlate(x, x, mode="full")[-n:]
    assert np.allclose(
        r, np.array([(x[: n - k] * x[-(n - k) :]).sum() for k in range(n)])
    )
    result = r / (variance * (np.arange(n, 0, -1)))
    return result


def savage_buford_formula(Vs, X, d):
    """
    X : distance à la faille , un iterable pour toutes le profil,
    un nombre pour la longeur max

    d : profondeur de la faille
    retourne X , et Vdeform(X)

    X et d doivent être dans la même unité, Vs pas forcément
    """

    if not utils.is_iterable(X):
        X = np.arange(-X, X, 1)
    return X, (Vs / np.pi) * np.arctan2(X, d)


def R2_calc(y_obs, y_fit, with_r2_bis=False):
    # https://en.wikipedia.org/wiki/Coefficient_of_determination
    ybar = np.mean(y_obs)
    SStot = np.sum((y_obs - ybar) ** 2)
    SSreg = np.sum((y_fit - ybar) ** 2)
    SSres = np.sum((y_obs - y_fit) ** 2)

    r2 = 1.0 - (SSres / SStot)
    r2bis = SSreg / SStot

    if not with_r2_bis:
        return r2
    else:
        return r2, r2bis


def R2_from_a_line_regress(Xobs, Yobs, a, b):
    # https://en.wikipedia.org/wiki/Coefficient_of_determination
    Xfit, Yfit = stats.linear_reg_getvalue(Xobs, a, b)
    r2 = R2_calc(Yobs, Yfit)
    return r2


def project_point_on_plan(N, M, A):
    """
    Project a point on a Plan

    Parameters
    ----------
    N : 3-Array (Vector)
        Vector describing the plan.
    M : 3-Array (Vector)
        Point of the plan.
    A : 3-Array (Vector)
        Point we want to project.

    Returns
    -------
    p : 3-Array (Vector)
        Projected point.

    Note
    ----
    It can be an ambiguity on the sign

    References
    ----------
    https://fr.wikipedia.org/wiki/Distance_d%27un_point_%C3%A0_un_plan
    https://www.mathematex.fr/viewtopic.php?t=923
    """

    ##dist
    d = np.abs(N.dot(A - M)) / np.linalg.norm(N)

    ## Normalized N
    Nnorm = N / np.linalg.norm(N)

    P = Nnorm * d

    return P


#  ______                _   _                _____                                         _
# |  ____|              | | (_)              / ____|                                       | |
# | |__ _   _ _ __   ___| |_ _  ___  _ __   | |  __ _ __ __ ___   _____ _   _  __ _ _ __ __| |
# |  __| | | | '_ \ / __| __| |/ _ \| '_ \  | | |_ | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
# | |  | |_| | | | | (__| |_| | (_) | | | | | |__| | | | (_| |\ v /  __/ |_| | (_| | | | (_| |
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \_____|_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
#                                                                        __/ |
#                                                                       |___/


def helmert_trans_legacy(
    Xa, params="itrf2008_2_etrf2000", invert=True, workepoc=2009.0
):
    """
    This function is highly unstable and should be used only for legacy reason


    Parameters
    ----------
    Xa : TYPE
        DESCRIPTION.
    params : TYPE, optional
        DESCRIPTION. The default is 'itrf2008_2_etrf2000'.
    invert : TYPE, optional
        DESCRIPTION. The default is True.
    workepoc : TYPE, optional
        DESCRIPTION. The default is 2009..

    Returns
    -------
    Xb : TYPE
        DESCRIPTION.

    Notes
    -----
    NB 1 : http://etrs89.ensg.ign.fr/memo-V8.pdf
    NB 2 :
    Transformation inverse : * -1 pour les paramètres
    cf https://en.wikipedia.org/wiki/Helmert_transformation

    optimisé pour RGF93 => ITRF2008 (d'ou le invert = True & workepoc = 2009.)

    NB3 : Attention losque l'on compare avec
    le convertisseur EUREF avec une vitesse
    parce que elle aussi est modifiée dans la conversion ETRS => ITRS.
    Conclusion : le faire en 2 étapes:
    ETRS => ITRS dans la même epoc
    ITRS epoc 1 => ITRS epoc 2

    NE MARCHE PAS PARCE BESOIN DU RATE(TAUX) DES PARAMS D'HELMERT !!!!!!
    (160923)

    if manual
    then params is a tuple
    params = (t1,t2,t3,dab,r1,r2,r3)

    """

    if invert:
        inver = -1.0
    else:
        inver = 1.0

    mas2rad = 0.0000000048481368
    mas2rad = 4.8481368111e-6 * 1e-3

    if params == "itrf2008_2_etrf2000":

        t1rate = 0.1 * 10**-3
        t2rate = 0.1 * 10**-3
        t3rate = -1.8 * 10**-3
        dabrate = 0.08 * 10**-9
        r1rate = 0.081 * mas2rad
        r2rate = 0.490 * mas2rad
        r3rate = -0.792 * mas2rad

        t1 = (52.1 * 10**-3 + t1rate * (workepoc - 2000.0)) * inver
        t2 = (49.3 * 10**-3 + t2rate * (workepoc - 2000.0)) * inver
        t3 = (-58.5 * 10**-3 + t3rate * (workepoc - 2000.0)) * inver
        dab = (1.34 * 10**-9 + dabrate * (workepoc - 2000.0)) * inver
        r1 = (0.891 * mas2rad + r1rate * (workepoc - 2000.0)) * inver
        r2 = (5.39 * mas2rad + r2rate * (workepoc - 2000.0)) * inver
        r3 = (-8.712 * mas2rad + r3rate * (workepoc - 2000.0)) * inver

    elif params == "itrf2000_2_etrf2000":
        t1 = 54.0 * 10**-3 * inver
        t2 = 51.0 * 10**-3 * inver
        t3 = -48.0 * 10**-3 * inver
        dab = 0.0 * 10**-9 * inver
        r1 = 0.891 * mas2rad * inver
        r2 = 5.390 * mas2rad * inver
        r3 = -8.712 * mas2rad * inver

    R = np.matrix([[dab, -r3, r2], [r3, dab, -r1], [-r2, r1, dab]])

    Xb = Xa + np.matrix([t1, t2, t3]) + np.dot(R, Xa)
    Xb = np.squeeze(np.array(Xb))

    return Xb


def mat_poids(Sinp, Ninp, fuvinp=1):
    """
    discontinued
    """
    # Sinp : liste des Sigmas sig = sqrt(var)
    # Ninp : liste de la taille de chaque blocs (obs)
    # fuvinp = 1 : facteur unitaire de variance

    if len(Sinp) != len(Ninp):
        raise Exception("S et N de taille differente")

    Ktemp = []

    for i in range(len(Sinp)):
        Ktemp.append(np.eye(Ninp[i]) * Sinp[i] ** 2)

    K = np.linalg.block_diag(*Ktemp)
    Q = (1 / fuvinp) * K
    P = np.linalg.inv(Q)

    return K, Q, P
