#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 13:39:48 2019

@author: psakicki
"""

from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

def ECI_2_kepler_elts(P,V,rad2deg=True,
                      mu=3.9860044188e14):
    Pnorm = np.linalg.norm(P)
    Vnorm = np.linalg.norm(V)
    
    ### orbital momentum vector H
    H = np.cross(P,V)
    Hnorm = np.linalg.norm(H)
    
    ### eccentricity vector Ecc
    Ecc = (np.cross(V,H)/mu) - (P/Pnorm)
    Ecc_norm = np.linalg.norm(Ecc)
    
    ### true anomaly ν (nu)
    N = np.cross(np.array([0,0,1]),H)
    Nnorm = np.linalg.norm(N)
    if np.dot(P,V) >= 0:
        nu = np.arccos(np.dot(Ecc,P)/(Ecc_norm*Pnorm))
    else:
        nu = 2*np.pi - np.arccos(np.dot(Ecc,P)/(Ecc_norm*Pnorm))        
    
    ### orbit inclination i
    i = np.arccos(H[2]/Hnorm)
    
    ### orbit eccentricity ecc
    ecc = np.linalg.norm(Ecc)
    
    ### eccentric anomaly E (called e here)
    #e1 = 2*np.arctan(np.tan(nu/2)/np.sqrt((1+ecc)/(1-ecc)))
    e2 = 2*np.arctan2(np.tan(nu/2),np.sqrt((1+ecc)/(1-ecc)))
    e=e2
    
    ### longitude of the ascending node Ω (omega_LAN)
    if N[1] >= 0:
        omega_LAN = np.arccos(N[0]/Nnorm)
    else:
        omega_LAN = 2*np.pi - np.arccos(N[0]/Nnorm)
    
    ### argument of periapsis ω (omega_periarg)
    if Ecc[2] >= 0:
        omega_periarg = np.arccos((np.dot(N,Ecc))/(Nnorm * Ecc_norm))
    else:
        omega_periarg = 2*np.pi - np.arccos((np.dot(N,Ecc))/(Nnorm * Ecc_norm))
    
    ### mean anomaly m, (Kepler’s Equation)
    m = e - ecc*np.sin(e)
    
    ### semi-major axis a
    a = ((2/Pnorm) - (Vnorm**2 / mu))**-1

    if rad2deg:
        i=np.rad2deg(i)
        omega_periarg=np.rad2deg(omega_periarg)
        omega_LAN=np.rad2deg(omega_LAN)
        m=np.rad2deg(m)

    return a,ecc,i,omega_periarg,omega_LAN,m