#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 10:30:57 2024

Filière ING3 - PPMD - Traitement de la mesure de phase

@author: Samuel Nahmani (1,2)
https://www.ipgp.fr/annuaire/nahmani/)
contact : nahmani@ipgp.fr ou samuel.nahmani@ign.fr
(1) Université Paris Cité, Institut de physique du globe de Paris, CNRS, IGN, F-75005 Paris, France.
(2) Univ Gustave Eiffel, ENSG, IGN, F-77455 Marne-la-Vallée, France. 

Version: 1.0
Dépendances: numpy, geodezyx, datetime

"""
#%%

"""
Conversion du code Matlab de Mahooti (2024) en Python
Meysam Mahooti (2024). Klobuchar Ionospheric Delay Model (https://www.mathworks.com/matlabcentral/fileexchange/59530-klobuchar-ionospheric-delay-model), MATLAB Central File Exchange. Retrieved February 6, 2024. 

"""


#%%
import numpy as np
from geodezyx import conv
import datetime as dt


def klobuchar(phi, lambda_, elev, azimuth, tow, alpha, beta):
    """
    Compute ionospheric range correction for GPS L1 frequency using Klobuchar model.
    
    Translation to Python of Meysam Mahooti (2024). Klobuchar Ionospheric Delay Model 
    https://www.mathworks.com/matlabcentral/fileexchange/59530-klobuchar-ionospheric-delay-model
    MATLAB Central File Exchange. Retrieved February 7, 2024.
    
    Parameters
    ----------
    phi : float
        Geodetic latitude of receiver (degrees).
    lambda_ : float
        Geodetic longitude of receiver (degrees).
    elev : float
        Elevation angle of satellite (degrees).
    azimuth : float
        Geodetic azimuth of satellite (degrees).
    tow : float
        Time of Week (seconds).
    alpha : array_like
        The coefficients of a cubic equation representing the amplitude 
        of the vertical delay (4 coefficients - 8 bits each).
    beta : array_like
        The coefficients of a cubic equation representing the period 
        of the model (4 coefficients - 8 bits each).
    
    Returns
    -------
    d_ion1 : float
        Ionospheric slant range correction for the L1 frequency (metres).
    
    References
    ----------
    Klobuchar, J.A., (1996) "Ionospheric Effects on GPS", in
    Parkinson, Spilker (ed), "Global Positioning System Theory and
    Applications, pp.513-514.
    
    ICD-GPS-200, Rev. C, (1997), pp. 125-128
    
    NATO, (1991), "Technical Characteristics of the NAVSTAR GPS",
    pp. A-6-31   -   A-6-33
    """

    # Constants
    c = 2.99792458e8  # Speed of light
    deg2semi = 1 / 180  # Degrees to semicircles
    semi2rad = np.pi  # Semicircles to radians
    deg2rad = np.pi / 180  # Degrees to radians

    # Convert angles to radians or semicircles
    a = azimuth * deg2rad  # Azimuth in radians
    e = elev * deg2semi    # Elevation angle in semicircles

    # Earth Centered angle
    psi = 0.0137 / (e + 0.11) - 0.022
    lat_i = phi * deg2semi + psi * np.cos(a)  # Subionospheric latitude
    lat_i = np.clip(lat_i, -0.416, 0.416)    # Clipping lat_i

    # Subionospheric longitude
    long_i = lambda_ * deg2semi + (psi * np.sin(a) / np.cos(lat_i * semi2rad))

    # Geomagnetic latitude
    lat_m = lat_i + 0.064 * np.cos((long_i - 1.617) * semi2rad)

    # Time
    t = 4.32e4 * long_i + tow
    t = t % 86400  # Seconds of day

    # Slant factor
    s_f = 1 + 16 * (0.53 - e) ** 3

    # Period of model
    per = beta[0] + beta[1] * lat_m + beta[2] * lat_m ** 2 + beta[3] * lat_m ** 3
    per = max(per, 72000)

    # Phase of the model
    x = 2 * np.pi * (t - 50400) / per

    # Amplitude of the model
    amp = alpha[0] + alpha[1] * lat_m + alpha[2] * lat_m ** 2 + alpha[3] * lat_m ** 3
    amp = max(amp, 0)

    # Ionospheric correction
    if np.abs(x) > 1.57:
        d_ion1 = s_f * 5e-9
    else:
        d_ion1 = s_f * (5e-9 + amp * (1 - x ** 2 / 2 + x ** 4 / 24))

    # Ionospheric slant range correction
    d_ion1 = c * d_ion1
    return d_ion1

if __name__ == "__main__":
    
   # Coefficients alpha et beta
    alpha = np.array([0.382e-7, 0.149e-7, -0.179e-6, 0.0])
    beta = np.array([0.143e6, 0.0, -0.328e6, 0.113e6])
    
    # Temps UTC
    year = 2000
    month = 1
    day = 1
    hour = 20
    minute = 45
    second = 0
    
    # Autres paramètres
    elev = 20
    azimuth = 210
    phi = 40
    lambda_ = 260
    
    # Calcul du TOW
    # Exemple d'utilisation
    t = dt.datetime(year, month, day, hour, minute, second)
    week, tow = conv.dt2gpstime(t,secinweek=True)
    print("GPS Week:", week)
    print("GPS Seconds:", tow)
 
    
    # Correction ionosphérique L1
    dIon1 = klobuchar(phi, lambda_, elev, azimuth, tow, alpha, beta)
    print("Ionospheric L1 Correction: ", dIon1)
