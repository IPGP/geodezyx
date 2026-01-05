# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.geodyn contains functions to determine the
Euler pole of a tectonic plate.

it can be imported directly with:
from geodezyx import geodyn

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox

This module is based on the work of :
Goudarzi, M. A., Cocard, M., & Santerre, R. (2014).
EPC: Matlab software to estimate Euler pole parameters. GPS Solutions, 18(1),
153–162. https://doi.org/10.1007/s10291-013-0354-4
"""

import numpy as np

#  ______      _             _____      _
# |  ____|    | |           |  __ \    | |
# | |__  _   _| | ___ _ __  | |__) |__ | | ___
# |  __|| | | | |/ _ \ '__| |  ___/ _ \| |/ _ \
# | |___| |_| | |  __/ |    | |  | (_) | |  __/
# |______\__,_|_|\___|_|    |_|   \___/|_|\___|
#

### Euler Pole determination


def euler_pole_calc(
    lat_ref, long_ref, vn_ref, ve_ref, incvn_ref=None, incve_ref=None, rt=6.3710088e6
):
    """
    Compute the Euler pole of a set of reference points

    Parameters
    ----------
    lat_ref,long_ref : list or numpy.array
        latitude and longitude of the reference points (deg)

    vn_ref,ve_ref : list or numpy.array
        north and east velocities of the reference points (m/yr)

    incvn_ref,incve_ref : list or numpy.array
        uncertainties on north and east velocities of the reference points (m/yr)

    rt : float
        Earth Radius
        IUGG value = 6.3710088e6
        IAU value  = 6.378e6

    Returns
    -------
    w : numpy.array
        Euler vector (rad/yr)

    wratedeg : float
        Rate of rotation (deg/Myr)

    wlat,wlong : float
        latitude and longitude of the Euler pole (deg)

    wwmat : numpy.array
        weight matrix (for debug)

    desmat : numpy.array
        design matrix (for debug)

    nrmatinv : 2-tuple
        output of scipy's lstsq fct (for debug)

    Source
    ------
    based on
    Goudarzi, M. a., Cocard, M., & Santerre, R. (2014).
    EPC: Matlab software to estimate Euler pole parameters. GPS Solutions, 18(1), 153–162.
    https://doi.org/10.1007/s10291-013-0354-4

    Notes
    -----
    w is the common element for all Euler pole functions
    Should remain in rad/yr
    Written by C. Geisert (ENSTA/LIENSs) - 2017
    """
    if (incve_ref is None) or (incvn_ref is None):
        incve_ref = np.ones(len(ve_ref))
        incvn_ref = np.ones(len(vn_ref))

    all_pos_ref = np.column_stack((lat_ref, long_ref))
    coords = (
        all_pos_ref * np.pi / 180
    )  # [rad] position des stations de la plaque de ref lat,long [°] > lat,long [rad]
    dm = topo2dm(coords)  # [rad]
    desmat = rt * dm  # design matrix [m] [rad]
    tdesmat = np.transpose(desmat)  # [m] [rad]
    n = len(all_pos_ref)
    vel_vector = list(
        np.zeros([2 * n, 1])
    )  # vn,ve ; vn, ve alternes  vecteur (vn,ve)0 puis (vn,ve)1, puis stations suivantes..
    for index in range(0, 2 * n - 1, 2):
        vel_vector[index] = vn_ref[int(index / 2)]
        vel_vector[index + 1] = ve_ref[int(index / 2)]  # vn_ref, ve_ref > vel_vector
    s = [0] * len(incvn_ref)
    covars = np.transpose(np.array([incvn_ref, incve_ref, s]))  # [m/yr]
    wwmat = covarvec2wtmat(covars)  # weight matrix of observations  #[m/yr]
    wwmat = np.array(wwmat[0])
    nrmat = tdesmat.dot(wwmat).dot(desmat)  # [m] [rad]
    lennrmat = np.shape(nrmat)  # []
    nrmatinv = np.linalg.lstsq(nrmat, np.eye(lennrmat[0]), rcond=None)  # [m] [│rad]
    a = (tdesmat.dot(wwmat)).dot(vel_vector)  # [m, rad] * [m/yr]
    w = nrmatinv[0].dot(a)  # taille 3x3  # [m, rad] * [m/yr]  #   [rad/yr]
    #    wx=w[0]  #[rad/ yr]
    #    wy=w[1]  #[rad/ yr]
    #    wz=w[2]  #[rad/ yr]
    #    wrate=np.sqrt(wx**2+wy**2+wz**2)*1e6    # [rad/yr]> [rad/Myr]
    #    wratedeg=wrate*180/np.pi  # [deg/Myr]
    #    wlat=(np.arctan(wz/(np.sqrt(wx**2+wy**2))))*180/np.pi # [°]
    #    wlong=(np.arctan(wy/wx))*180/np.pi # [°]

    wlat, wlong, wrate, wratedeg = euler_pole_vector_to_latlongrate(w)

    return w, wratedeg, wlat, wlong, wwmat, desmat, nrmatinv


def euler_pole_vector_to_latlongrate(w):
    """
    Convert Euler pole vector to latitude longitude and rate

    Parameters
    ----------
    w : numpy.array
        Pole of rotation vector computed by euler_pole_calc (rad/yr)

    Returns
    -------
    wlat,wlong : float
        latitude and longitude of the Euler pole (deg)

    wrate : float
        Rate of rotation (rad/Myr)

    wratedeg : float
        Rate of rotation (deg/Myr)

    Notes
    -----
    w is the common element for all Euler pole functions
    Should remain in rad/yr
    Written by P. Sakic based on C. Geisert's work
    """
    wx = w[0]  # [rad/ yr]
    wy = w[1]  # [rad/ yr]
    wz = w[2]  # [rad/ yr]
    wrate = np.sqrt(wx**2 + wy**2 + wz**2) * 1e6  # [rad/yr]> [rad/Myr]
    wratedeg = wrate * 180 / np.pi  # [deg/Myr]
    wlat = (np.arctan(wz / (np.sqrt(wx**2 + wy**2)))) * 180 / np.pi  # [°]
    wlong = (np.arctan(wy / wx)) * 180 / np.pi  # [°]
    return wlat, wlong, wrate, wratedeg


def euler_vels_relative_to_ref(
    w,
    lat_itrf,
    long_itrf,
    vn_itrf,
    ve_itrf,
    incvn_itrf=None,
    incve_itrf=None,
    rt=6.3710088e6,
):
    """
    Compute relative velocities of points with respect to a reference plate/Euler pole

    Parameters
    ----------
    w : numpy.array
        Pole of rotation vector computed by euler_pole_calc (rad/yr)

    lat_itrf,long_itrf : list or numpy.array
        latitude and longitude of the points (deg)

    vn_itrf,ve_itrf : list or numpy.array
        north and east velocities of the points (m/yr)

    incvn_itrf,incve_itrf : list or numpy.array
        uncertainties on north and east velocities of the points (m/yr)

    rt : float
        Earth Radius
        IUGG value = 6.3710088e6
        IAU value  = 6.378e6


    Returns
    -------
    vel_reltoref : numpy.array
        relative velocities of points with respect to the Euler pole w

    Notes
    -----
    w is the common element for all Euler pole functions
    Should remain in rad/yr
    Written by C. Geisert (ENSTA/LIENSs) - 2017
    """
    if (incvn_itrf is None) or (incve_itrf is None):
        incvn_itrf = np.ones(len(vn_itrf))
        incve_itrf = np.ones(len(ve_itrf))

    all_pos = np.column_stack((lat_itrf, long_itrf))
    nstation = len(all_pos)
    vel_rot = [0] * nstation
    for i in range(nstation):  #
        mat = rt * np.array(
            [
                np.array(
                    [
                        np.sin(all_pos[i][1] * np.pi / 180),
                        -np.cos(all_pos[i][1] * np.pi / 180),
                        0,
                    ]
                ),
                np.array(
                    [
                        -np.sin(all_pos[i][0] * np.pi / 180)
                        * np.cos(all_pos[i][1] * np.pi / 180),
                        -np.sin(all_pos[i][0] * np.pi / 180)
                        * np.sin(all_pos[i][1] * np.pi / 180),
                        np.cos(all_pos[i][0] * np.pi / 180),
                    ]
                ),
            ]
        )  # [°] > [rad] pour le calcul

        vel_rot[i] = mat.dot(w)  #   [m] cross [rad/yr]
    # vel_rot=np.array(vel_rot)*1000  # [m/an]  > [mm/an]
    v_itrf = np.transpose(
        np.array([vn_itrf, ve_itrf])
    )  # les vitesses de toutes les stations [m/yr]
    vel_reltoref = v_itrf - vel_rot  # [m/an]

    return vel_reltoref


def euler_pole_vector_from_latlonrate(
    wlat, wlong, wrate, return_w_in_deg_per_megayr=False
):
    """
    Compute the Euler vector from the pole latitude, longitude and rate

    Parameters
    ----------
    wlat,wlong : float
        latitude and longitude of the Euler pole (deg)

    wrate : float
        rate of the Euler pole (rad/Myr)

    return_w_in_deg_per_megayr : bool
        if True, return w in deg/Myr instead of rad/yr

    Returns
    -------
    w : numpy.array
        Pole of rotation vector (rad/yr or deg/Myr)

    Notes
    -----
    w is the common element for all Euler pole functions
    Should remain in rad/yr per default
    Written by C. Geisert (ENSTA/LIENSs) - 2017
    """
    # OUTPUT : w [rad/year]

    ### Conversion deg => rad coordinates
    wlong = np.deg2rad(wlong)
    wlat = np.deg2rad(wlat)

    wx = (wrate / 1e6) * np.cos(wlat) * np.cos(wlong)  # *np.pi/180  # [rad/yr]
    wy = (wrate / 1e6) * np.cos(wlat) * np.sin(wlong)  # *np.pi/180  # [rad/yr]
    wz = (wrate / 1e6) * np.sin(wlat)  # *np.pi/180                # [rad/yr]
    w = [0] * 3

    if return_w_in_deg_per_megayr:
        k = (180 / np.pi) * 10**6
    else:
        k = 1.0

    w[0] = wx * k
    w[1] = wy * k
    w[2] = wz * k

    return w


# %% ____Estimation de l'incertitude lors du calcul des EPP
def euler_pole_quality(w, vn_ref, ve_ref, nrmatinv, desmat, wwmat, pretty_output=True):
    """
    Compute the uncertainties of the Euler pole determination

    Parameters
    ----------
    w : numpy.array
        Pole of rotation vector computed by euler_pole_calc

    vn_ref,ve_ref : list or numpy.array
        north and east velocities of the reference points (m/yr)

    nrmatinv : 2-tuple
        output of scipy's lstsq fct from euler_pole_calc

    wwmat : numpy.array
        weight matrix from euler_pole_calc

    desmat : numpy.array
        design matrix from euler_pole_calc

    pretty_output : bool
        if True, convert sigma_ww_latlon to pertinent units directly,
        returns raw units instead

    Returns
    -------
    sigma_ww : numpy.array
        Uncertainty on the Euler vector

    sigma_ww_latlon : numpy.array
        Uncertainty on the Euler pole : [rate_sigma, lat_sigma, long_sigma]
        if pretty_output == True :  [deg/Myr,deg,deg]
        if pretty_output == False : [rad/yr,rad,rad]

    dv_topo3 : numpy.array
        Residual velocities for the references points

    wrmse : float
        weigthed RMS on Residual velocities (m)

    wrmse_norm : float
        nomalized weigthed RMS on Residual velocities (m)

    rmse : float
        unweigthed RMS on Residual velocities (m)

    apost_sigma : float
        a-posteriori sigma (m)

    Notes
    -----
    w is the common element for all Euler pole functions
    Should remain in rad/yr
    """
    v_ref = np.transpose(
        np.array([vn_ref, ve_ref])
    )  # les vitesses de toutes les stations [mm/yr]
    estm_vel = desmat.dot(w)
    estm_vel = np.transpose(np.array([estm_vel[0::2], estm_vel[1::2]]))
    estm_vel_diff = v_ref - estm_vel
    dv_topo = np.transpose(estm_vel_diff)
    dv_topo = np.transpose(dv_topo)
    nd = len(dv_topo)
    dv_topo = np.append(dv_topo[:, 0], dv_topo[:, 1])
    dv_topo2 = np.zeros([2 * nd])
    for index in range(0, 2 * nd - 1, 2):
        dv_topo2[index] = dv_topo[int(index / 2)]
        dv_topo2[index + 1] = dv_topo[int(index / 2 + nd)]
    df_value = len(dv_topo2) - 3
    s0_2 = (np.transpose(dv_topo2).dot(wwmat).dot(dv_topo2)) / df_value
    c_ww = s0_2 * (nrmatinv[0])
    omega_length = np.sqrt(np.transpose(w).dot(w))
    # H
    h00 = w[0] / omega_length
    h01 = w[1] / omega_length
    h02 = w[2] / omega_length
    h10 = (-w[0] * w[2]) / (omega_length**2 * np.sqrt(w[0] ** 2 + w[1] ** 2))
    h11 = (-w[1] * w[2]) / (omega_length**2 * np.sqrt(w[0] ** 2 + w[1] ** 2))
    h12 = -np.sqrt(w[0] ** 2 + w[1] ** 2) / omega_length**2
    h20 = -w[1] / (w[0] ** 2 + w[1] ** 2)
    h21 = w[0] / (w[0] ** 2 + w[1] ** 2)
    h22 = 0
    H = np.array([[h00, h01, h02], [h10, h11, h12], [h20, h21, h22]])
    c_ww_latlong = H.dot(c_ww).dot(np.transpose(H))

    ### OUTPUT GENERATION
    sigma_ww = np.sqrt(np.diag(c_ww))
    sigma_ww_latlong = np.sqrt(np.diag(c_ww_latlong))
    if pretty_output:
        rate_sigma = (sigma_ww_latlong[0]) * 1e6 * 180 / np.pi  # [rad/yr]> [°/Myr]
        lat_sigma = (sigma_ww_latlong[1]) * 180 / np.pi  # [rad] >[°]
        long_sigma = (sigma_ww_latlong[2]) * 180 / np.pi  # [rad]>[°]
        sigma_ww_latlong = np.array([rate_sigma, lat_sigma, long_sigma])
    dv_topo3 = np.reshape(dv_topo2, (int(len(dv_topo2) / 2), 2))

    wrmse = np.sqrt((np.transpose(dv_topo2).dot(wwmat)).dot(dv_topo2) / len(dv_topo2))
    wrmse_norm = np.sqrt(
        (np.transpose(dv_topo2).dot(wwmat / (np.sum(np.diag(wwmat))))).dot(dv_topo2)
        / len(dv_topo2)
    )
    rmse = np.sqrt(np.transpose(dv_topo2).dot(dv_topo2) / len(dv_topo2))
    apost_sigma = np.sqrt(s0_2)

    return sigma_ww, sigma_ww_latlong, dv_topo3, wrmse, wrmse_norm, rmse, apost_sigma


def topo2dm(coords):
    """
    Internal fuction for euler_pole_calc
    """
    n = len(coords)
    dm = np.zeros([2 * n, 3])
    for index in range(0, 2 * n - 1, 2):
        dm[index, 0] = np.sin(coords[int((index + 1) / 2)][1])
        dm[index, 1] = -np.cos(coords[int((index + 1) / 2)][1])
        dm[index + 1, 0] = -np.sin(coords[int((index + 1) / 2)][0]) * np.cos(
            coords[int((index + 1) / 2)][1]
        )
        dm[index + 1, 1] = -np.sin(coords[int((index + 1) / 2)][0]) * np.sin(
            coords[int((index + 1) / 2)][1]
        )
        dm[index + 1, 2] = np.cos(coords[int((index + 1) / 2)][0])
    return dm


def covarvec2wtmat(covars):
    """
    Internal fuction for euler_pole_calc
    """
    wtmat = []
    aa = np.shape(covars)
    cc = np.zeros([aa[0] * 2, aa[0] * 2])
    a = np.shape(cc)
    for idx in range(0, a[0], 2):
        cc[idx, idx] = covars[int(idx / 2)][0] ** 2
        cc[idx + 1, idx + 1] = covars[int(idx / 2)][1] ** 2

        cc[idx, idx + 1] = covars[int(idx / 2)][2]
        cc[idx + 1, idx] = covars[int(idx / 2)][2]

    cc[np.isnan(cc)] = 0
    wtmat = np.linalg.lstsq(cc, np.eye(a[0]), rcond=None)
    return wtmat
