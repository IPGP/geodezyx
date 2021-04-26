# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 13:37:40 2019

@author: psakicki

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics

Copyright (C) 2019 Pierre Sakic (GFZ, pierre.sakic@gfz-postdam.de)
GitHub repository :
https://github.com/PierreS1/GeodeZYX-Toolbox-Lite

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as np
import scipy
import re
    
import conversion_general as cnv_gen

import sys
sys.dont_write_bytecode = True

#### Coordinates conversion
    
#   _____                    _ _             _               _____                              _             
#  / ____|                  | (_)           | |             / ____|                            (_)            
# | |     ___   ___  _ __ __| |_ _ __   __ _| |_ ___  ___  | |     ___  _ ____   _____ _ __ ___ _  ___  _ __  
# | |    / _ \ / _ \| '__/ _` | | '_ \ / _` | __/ _ \/ __| | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \ 
# | |___| (_) | (_) | | | (_| | | | | | (_| | ||  __/\__ \ | |___| (_) | | | \ V /  __/ |  \__ \ | (_) | | | |
#  \_____\___/ \___/|_|  \__,_|_|_| |_|\__,_|\__\___||___/  \_____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
#                                                                                                             
#                                                                                                             

def wnorm(phi,a=6378137.,e2=0.00669438003):
    """
    Compute the Ellipsoid "Grande Normale"
    
    Source
    ------
    ALG0021 in NTG_71 (IGN Lambert Projection Tech document)
    Based on PYACS of J.-M. Nocquet
    """
    from numpy import sqrt,sin
    #e=sqrt(e2)
    wd=sqrt(1 - e2*sin(phi)**2)
    result=a/wd
    return result

def normal_vector(phi,llambda,angle='deg',normalized=True):
    """
    Compute the Ellipsoid "Normale"
    
    Source
    ------
    P. Bosser (2012), Géométrie de l'Ellipsoïde, p27
    """

    if angle == 'deg':
        llambda = np.deg2rad(llambda)
        phi = np.deg2rad(phi)

    a = np.cos(llambda) * np.cos(phi)
    b = np.sin(llambda) * np.cos(phi)
    c = np.sin(phi)

    n = np.array([a,b,c])

    if normalized:
        n = n / np.linalg.norm(n)

    return n


def GEO2XYZ(phi,llambda,he,angle='deg',
            a=6378137.,e2=0.00669438003):
    """
    Coordinates conversion

    FLH Geographic => XYZ ECEF Cartesian

    Parameters
    ----------
    
    phi,llambda,he : numpy.array of floats
        latitude (deg/rad), longitude (deg/rad), height (m)   
    
    angle : string
        describe the angle input type : 'deg' or 'rad'
    
    a, e2 : floats
        ellipsoid parameters (WGS84 per default)
    

    Returns
    -------
    X,Y,Z : numpy.array of floats
        cartesian coordinates X,Y,Z in meters
        
    Source
    ------
    Based on PYACS of J.-M. Nocquet
    """

    if angle == 'deg':
        llambda = np.deg2rad(llambda)
        phi = np.deg2rad(phi)

    #f=1.0 - sqrt(1-e2)

    wn=wnorm(phi)

    X=(wn+he)*np.cos(phi)*np.cos(llambda)
    Y=(wn+he)*np.cos(phi)*np.sin(llambda)
    Z=(wn*(1.0-e2)+he)*np.sin(phi)

    return X,Y,Z


def XYZ2GEO(x,y,z,outdeg=True,
            A=6378137.,E2=0.00669438003):
    """
    Coordinates conversion

    XYZ ECEF Cartesian => FLH Geographic

    Parameters
    ----------
    x,y,z : numpy.array of floats
        cartesian coordinates X,Y,Z in meters

    outdeg : bool
        if True, degrees output for the angle, else radian
    
    A, E2 : floats
        ellipsoid parameters (WGS84 per default)
    

    Returns
    -------
    latitude, longitude, height : numpy.array of floats
        latitude (deg/rad), longitude (deg/rad), height (m)   
        
    Source
    ------
    Based on PYACS of J.-M. Nocquet
    """

    F=1.0 - np.sqrt(1-E2)

    TP = np.sqrt(x**2 + y**2)
    R = np.sqrt(TP**2 + z**2)

    TMU = np.arctan2( (z/TP)*((1. - F) + (E2*A)/R) ,1)
    RLBDA = np.arctan2(y,x)

    S3 = (np.sin(TMU))**3
    C3 = (np.cos(TMU))**3
    T1 = z*(1 - F) + E2*A*S3
    T2 = (1 - F)*( TP - E2*A*C3 )

    RPHI = np.arctan2(T1,T2)

    RHE = TP*(np.cos(RPHI)) + z*(np.sin(RPHI))
    RHE = RHE - A*( np.sqrt(1-E2*(np.sin(RPHI)**2)) )

    if outdeg:
        RPHI = np.rad2deg(RPHI)
        RLBDA = np.rad2deg(RLBDA)

    return RPHI,RLBDA,RHE


def XYZ2ENU(dX,dY,dZ,lat0,lon0):
    """
    Coordinates conversion

    XYZ ECEF Geocentric => ENU Topocentic
    
    core function for XYZ2ENU_2, use the later in priority

    dXYZ = XYZrover - XYZref
    
    Parameters
    ----------
    dX,dY,dZ : floats or numpy.array of floats
        cartesian coordinates difference between the considered point(s)
        and the reference point 
        
    lat0,lon0 : floats or numpy.array of floats
        if they are iterable arrays (of the same size as dX,dY,dZ) 
        a different x0,y0,z0 will be applied for each dX,dY,dZ element

    Returns
    -------
    E,N,U : numpy.array of floats
        East North Up Component (m) w.r.t. x0,y0,z0
        
    Source
    ------
    https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
    
    Note
    ----
    This recursive fuction should be improved
    """
    
    ## Case one ref point per dXYZ
    if cnv_gen.is_iterable(lat0):
        E,N,U = [] , [] , []
        for dX_m,dY_m,dZ_m,lat0_m,lon0_m in zip(dX,dY,dZ,lat0,lon0):
            E_m , N_m , U_m = XYZ2ENU(dX_m,dY_m,dZ_m,lat0_m,lon0_m)
            E.append(E_m)
            N.append(N_m)
            U.append(U_m)
            
        return np.squeeze(np.array(E)) , \
               np.squeeze(np.array(N)) , \
               np.squeeze(np.array(U)) 
    
    # case onle ref point for all dXYZ
    else:
        f0 = np.deg2rad(lat0)
        l0 = np.deg2rad(lon0)
    
        R = np.array([[ -np.sin(l0) , np.cos(l0) , 0 ] ,
        [ -np.sin(f0)*np.cos(l0) , -np.sin(f0)*np.sin(l0) , np.cos(f0) ] ,
        [  np.cos(f0)*np.cos(l0) , np.cos(f0)*np.sin(l0)  , np.sin(f0)]])
    
        enu = np.dot(R,np.vstack((dX,dY,dZ)))
    
        E = enu[0,:]
        N = enu[1,:]
        U = enu[2,:]
    
        return E,N,U

def XYZ2ENU_2(X,Y,Z,x0,y0,z0):
    """
    Coordinates conversion

    XYZ ECEF Geocentric => ENU Topocentic
    
    Parameters
    ----------
    X,Y,Z : numpy.array of floats
        cartesian coordinates X,Y,Z in meters
        
    x0,y0,z0 : floats or numpy.array of floats
        coordinate of the wished topocentric origin point in the geocentric frame
        if they are iterable arrays (of the same size as X,Y,Z) 
        a different x0,y0,z0 will be applied for each X,Y,Z element

    Returns
    -------
    E,N,U : numpy.array of floats
        East North Up Component (m) w.r.t. x0,y0,z0
        
    Source
    ------
    https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
    """
    
    f0,l0,h0 = XYZ2GEO(x0,y0,z0,outdeg=1)
    dX = np.array(X) - x0
    dY = np.array(Y) - y0
    dZ = np.array(Z) - z0
    E,N,U = XYZ2ENU(dX,dY,dZ,f0,l0)
    return E,N,U


def XYZ2ENU_around_fix_pos(X,Y,Z):
    """
    Computes the mean of the X,Y,Z and
    return the topocentric ENU position around this mean position
    
    Parameters
    ----------
    X,Y,Z : numpy.array of floats
        cartesian coordinates X,Y,Z in meters
        
    Returns
    -------
    E,N,U : numpy.array of floats
        East North Up Component (m) w.r.t. the mean position
        
    x0 , y0 , z0 : floats
        coordinates of the mean position in the geocentic frame
    """
    x0 = np.nanmean(X)
    y0 = np.nanmean(Y)
    z0 = np.nanmean(Z)
    return XYZ2ENU_2(X,Y,Z,x0,y0,z0) , np.array([x0,y0,z0])


def ENU2XYZ(E,N,U,Xref,Yref,Zref):
    """
    Coordinates conversion

    ENU Topocentic => XYZ ECEF Geocentric 
    
    Parameters
    ----------
    X,Y,Z : numpy.array of floats
        cartesian coordinates X,Y,Z in meters
        
    x0,y0,z0 : floats
        coordinate of the topocentric origin point in the geocentric frame    

    Returns
    -------
    E,N,U : numpy.array of floats
        East North Up Component (m) w.r.t. x0,y0,z0
        
    Source
    ------
    https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
    """
    
    if cnv_gen.is_iterable(E):
        Xlist , Ylist , Zlist = [] , [] , []
        for e,n,u in zip(E,N,U):
            x,y,z = ENU2XYZ(e,n,u,Xref,Yref,Zref)
            Xlist.append(x)
            Ylist.append(y)
            Zlist.append(z)
        return np.array(Xlist) , np.array(Ylist) , np.array(Zlist)
    
    else:
        fr,lr,hr = XYZ2GEO(Xref,Yref,Zref)
        f0 = np.deg2rad(fr)
        l0 = np.deg2rad(lr)
    
        R = np.array([[ -np.sin(l0)            , np.cos(l0)             , 0          ],
                      [ -np.sin(f0)*np.cos(l0) , -np.sin(f0)*np.sin(l0) , np.cos(f0) ],
                      [  np.cos(f0)*np.cos(l0) , np.cos(f0)*np.sin(l0)  , np.sin(f0)]])
    
        R3 = R.T
        R3 = scipy.linalg.inv(R)
    
        ENU = np.vstack((E,N,U))
    
        xyz = np.dot(R3,ENU) #+ np.vstack((Xref,Yref,Zref))
    
        X = float(xyz[0]) + Xref
        Y = float(xyz[1]) + Yref
        Z = float(xyz[2]) + Zref
    
        return X,Y,Z


def ENU2XYZ_legacy(E,N,U,Xref,Yref,Zref):
    """
    KEPT FOR LEGACY REASONS, use ENU2XYZ

    diffère de ENU2XYZ pour d'obscure raisons, à investiguer !!!
    est laissé pour des scripts de conversion de GINS (170119)
    """

    fr,lr,hr = XYZ2GEO(Xref,Yref,Zref)
    f0 = np.deg2rad(fr)
    l0 = np.deg2rad(lr)

    R = np.array([[ -np.sin(l0)            ,  np.cos(l0)            ,          0 ] ,
                  [ -np.sin(f0)*np.cos(l0) , -np.sin(f0)*np.sin(l0) , np.cos(f0) ] ,
                  [  np.cos(f0)*np.cos(l0) ,  np.cos(f0)*np.sin(l0) , np.sin(f0) ]])

    R3 = R.T

    ENU = np.vstack((E,N,U))

    xyz = np.dot(R3,ENU) #+ np.vstack((Xref,Yref,Zref))

    X = float(xyz[0])
    Y = float(xyz[1])
    Z = float(xyz[2])

    return X,Y,Z


def sFLH2sXYZ(F,L,H,sF,sL,sH,ang='deg'):
    """
    Convert standard deviation
    Geographic FLH => Cartesian ECEF XYZ
    
    WARNING
    -------
    Inputs values are assumed as uncorrelated, which is not accurate
    Must be improved
    
    Source
    ------
    Linear Algebra, Geodesy, and GPS p332
    """
    
    
    if ang == 'deg':
        F  = np.deg2rad(F)
        L  = np.deg2rad(L)
        sF = np.deg2rad(sF)
        sL = np.deg2rad(sL)
    X,Y,Z  = GEO2XYZ(F,L,H)
    X2,Y2,Z2 = GEO2XYZ(F+sF,L+sL,H+sH)

    return np.abs(X-X2) , np.abs(Y-Y2) , np.abs(Z-Z2)

#def sXYZ2sFLH(X,Y,Z,sX,sY,sZ):
#    F,L,H    = XYZ2GEO(X,Y,Z)
#    F2,L2,H2 = XYZ2GEO(X+sX,Y+sY,Z+sZ)
#
#    return np.abs(F-F2) , np.abs(L-L2) , np.abs(H-H2)
#

def sFLH2sENU(F,L,H,sF,sL,sH,ang='deg'):
    """
    Convert standard deviation
    Geographic FLH => Cartesian Topocentric ENU
    
    WARNING
    -------
    Inputs values are assumed as uncorrelated, which is not accurate
    Must be improved
    
    Source
    ------
    Linear Algebra, Geodesy, and GPS p332
    """
    # conversion batarde du sigma FLH => sigma ENU
    # Par conversion des angles en distance
    if ang == 'deg':
        F = np.deg2rad(F)
        L = np.deg2rad(L)
        sF = np.deg2rad(sF)
        sL = np.deg2rad(sL)

    E2=0.00669438003
    #f=1.0 - np.sqrt(1-E2)
    A=6378137.

    # Je présume que Rpsi est le rayon du cercle tangent à l'ellipsoide à une
    # lattitude donnée (comm du 150710)
    Psi = np.arctan((1 - E2) * np.tan(F))
    Rpsi = A * np.sqrt(1 - E2) / np.sqrt(1 - E2 * np.cos(Psi)**2)

    # r est le rayon d'un petit cercle (un parallèle)
    r = np.cos(Psi) * (Rpsi + H)

    # On pourait simplifier par 2pi mais autant avoir toute la démarche
    sE = (np.pi * 2 * r    * sL) / (2 * np.pi)
    sN = (np.pi * 2 * Rpsi * sF) / (2 * np.pi)
    sU = sH
    return sE,sN,sU

def sENU2sFLH(F,L,H,sE,sN,sU,ang='deg',
              A=6378137.,E2=0.00669438003):
    """
    Convert standard deviation
    Cartesian Topocentric ENU => Geographic FLH 
    
    WARNING
    -------
    Inputs values are assumed as uncorrelated, which is not accurate
    Must be improved
    
    Source
    ------
    Linear Algebra, Geodesy, and GPS p332
    """
    # conversion batarde du sigma ENU => sigma FLH
    # Par conversion des angles en distance
    
    #f=1.0 - np.sqrt(1-E2)

    # Je présume que Rpsi est le rayon du cercle tangent à l'ellipsoide à une
    # lattitude donnée (comm du 150710)
    Psi = np.arctan((1 - E2) * np.tan(F))
    Rpsi = A * np.sqrt(1 - E2) / np.sqrt(1 - E2 * np.cos(Psi)**2)

    # r est le rayon d'un petit cercle (un parallèle)
    r = np.cos(Psi) * (Rpsi + H)

    # On pourait simplifier par 2pi mais autant avoir toute la démarche
    sL = sE *  (2 * np.pi) / (np.pi * 2 * r    )
    sF = sN *  (2 * np.pi) / (np.pi * 2 * Rpsi )
    sH = sU

    if ang=='deg':
        sL = np.rad2deg(sL)
        sF = np.rad2deg(sF)
    return sF,sL,sH


def sXYZ2sENU(X,Y,Z,sX,sY,sZ,sXY=0,sYZ=0,sXZ=0):
    """
    Convert standard deviation
    Cartesian ECEF XYZ => Cartesian Topocentric ENU 
    
    WARNING
    -------
    Inputs values are assumed as uncorrelated, which is not accurate
    Must be improved
    
    Source
    ------
    Linear Algebra, Geodesy, and GPS p332
    """

    #conversion "batarde" en assumant par défaut
    #que XYZ ne sont pas corrélé (pas bien ...)
    #Linear Algebra, Geodesy, and GPS p332

    SIGMAxyz = np.array([[sX**2,sXY,sXZ],
                         [sXY,sY**2,sYZ],
                         [sXZ,sYZ,sZ**2]])

    F,L,H = XYZ2GEO(X,Y,Z)

    C = C_ecef2enu(F,L,angtype='deg')

    SIGMAenu = np.dot(np.dot(C,SIGMAxyz),C.T)

    sE = np.sqrt(SIGMAenu[0,0])
    sN = np.sqrt(SIGMAenu[1,1])
    sU = np.sqrt(SIGMAenu[2,2])

    return sE,sN,sU


def ECI2RTN_or_RPY(P,V,C,out_rpy=False,rpy_theo_mode=False):
    """
    convert ECI coordinates in RTN (RIC) or RPY (Roll Pitch Yaw)

    Parameters
    ----------
    P : numpy.array
        3D vector, position of the ref object in ECI frame
    V : numpy.array
        3D vector, velocity of the ref object in ECI frame
    C : numpy.array
        3D Vector, coordinates in ECF frame that will be transformed
    out_rpy : bool
        if True output in RPY frame, RTN instead
    rpy_theo_mode : bool
        use the theoretical matrix composition, but wrong ans.
        empirically, only for debug !!

    Returns
    -------
    Cout : conversion of P in RTN or RPY ref. frame

    Source
    ------
    "Coordinate Systems", ASEN 3200 1/24/06 George H. Born
    """

    C_eci2rtn_mat  = C_eci2rtn(P,V)

    if not out_rpy:
        TransMat = C_eci2rtn_mat
    else:

        if not rpy_theo_mode:
            # Pour de très obscures raisons la composition est inversée
            # par rapport à l'ordre standard ... (241017)
            TransMat = np.dot(C_rtn2rpy().T,C_eci2rtn_mat.T)
        else:
            print("WARN : using the theoretical mode for RPY conversion, UNSTABLE & WRONG !")
            TransMat = np.dot(C_rtn2rpy(),C_eci2rtn_mat)

        # Mais reste compatible avec Wikipedia
        # https://en.wikipedia.org/wiki/Permutation_matrix#Permutation_of_rows_and_columns

    Cout = TransMat.dot( C )

    # EXEMPLE POUR DEBUG
    #     AAAAAAAAAAAAA yawerr,aux,nomi
    #       -1.606542300481316
    #       -1.606542300481316
    # CCCCCCCCCCCC sunors orientation du soleil
    #   2.379858191647978E-002   6.654859888132475E-001  -7.460308480163662E-001
    # DDDDDDDDDDDDD Inputs : xsat , vsat , xsun
    #   -15779.437215560150435   -23508.964967410494864     8644.453520237560952
    #        1.296220384065493       -1.917314252639830       -2.847293074880089
    #  -1.342404641197934E+008  -5.934080206107079E+007  -2.572342068927875E+007
    # EEEEEEEEEEEE Rcrs2ors
    #   3.532662370257567E-001  -7.688096793507927E-001   5.330195492501879E-001
    #  -5.225364448455557E-001   3.104619621186746E-001   7.941181766553120E-001
    #  -7.759888076421138E-001  -5.590572803191509E-001  -2.920042493231377E-001
    #    Rcrs2ors = np.array([[ 0.35330134, -0.76880968,  0.53301955],
    #                         [-0.52248415,  0.31046196,  0.79411818],
    #                         [-0.77600804, -0.55905728, -0.29200425]])

    #    Psat = np.array([   -15779.437215560150435 ,  -23508.964967410494864  ,   8644.453520237560952])
    #    Vsat = np.array([        1.296220384065493  ,     -1.917314252639830  ,     -2.847293074880089])
    #    Psun = np.array([  -1.342404641197934E+008 , -5.934080206107079E+007 , -2.572342068927875E+007])

    return Cout


def ECI2RTN(P,V,C):
    """
    legacy wrapper of ECI2RTN_or_RPY
    """
    return ECI2RTN_or_RPY(P,V,C,out_rpy=False)


def ECEF2ECI(xyz,utc_times):
    """ 
    Convert ECEF (Earth Centered Earth Fixed) positions to ECI (Earth Centered Inertial)
    positions

    Parameters
    ----------
    xyz : numpy.array of floats
        XYZ are cartesian positions in ECEF. Should have shape (N,3)
    
    utc_times : numpy.array of floats
        UTC_times are UTC timestamps, as datetime objects. Sould have shape (N)
    
    Returns
    -------
    eci : numpy.array of floats
        Earth Centered Inertial UTC_times are UTC_times, as datetime objects. Sould have shape (N)    

    Note
    ----
    Requires pyorbital module

    Theory
    ------   
    
     [X]    [C -S 0][X]
     [Y]  = [S  C 0][Y]
     [Z]eci [0  0 1][Z]ecf

     C and S are cos() and sin() of gmst (Greenwich Meridian Sideral Time)

    Source
    ------
    http://ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
    Inspired from satellite-js (https://github.com/shashwatak/satellite-js)
    """
    from pyorbital import astronomy
    # XYZ and utc_time must have the same shape
    #if not xyz.shape[:-1] == utc_times.shape:
    #    raise ValueError("shape mismatch for XYZ and utc_times (got {} and {})".format(xyz.shape[:-1],utc_times.shape))

    #    gmst = -1 * astronomy.gmst(utc_times) # EDIT 180430 : Why this -1 ??? removed because wrong ! ...
    gmst = 1 * astronomy.gmst(utc_times)

    eci = xyz.copy()
    eci[:,0] = xyz[:,0]*np.cos(gmst) - xyz[:,1]*np.sin(gmst)
    eci[:,1] = xyz[:,0]*np.sin(gmst) + xyz[:,1]*np.cos(gmst)
    return eci


#  _____       _        _   _               __  __       _        _               
# |  __ \     | |      | | (_)             |  \/  |     | |      (_)              
# | |__) |___ | |_ __ _| |_ _  ___  _ __   | \  / | __ _| |_ _ __ _  ___ ___  ___ 
# |  _  // _ \| __/ _` | __| |/ _ \| '_ \  | |\/| |/ _` | __| '__| |/ __/ _ \/ __|
# | | \ \ (_) | || (_| | |_| | (_) | | | | | |  | | (_| | |_| |  | | (_|  __/\__ \
# |_|  \_\___/ \__\__,_|\__|_|\___/|_| |_| |_|  |_|\__,_|\__|_|  |_|\___\___||___/
#                                                                                 
                                                                                 
### Rotation matrix

### The rotation matrices conventions are described in details in 
### Sakic (2016, PhD manuscript) p126
    
### It is based on the book of Grewal, M. S., Weill, L. R., & Andrews, A. P. (2007). 
### Global Positioning Systems, Inertial Navigation, and Integration. 
### Hoboken, NJ, USA : John Wiley and Sons, Inc.

def C_2D(theta,angtype='deg'):

    if angtype == 'deg':
        theta = np.deg2rad(theta)

    return np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta), np.cos(theta)]])

def C_ned2enu():
    return np.array([[0,1,0],[1,0,0],[0,0,-1]])

def C_enu2ned():
    return np.array([[0,1,0],[1,0,0],[0,0,-1]])

def C_ned2ecef(phi,llambda,angtype='deg'):

    if angtype == 'deg':
        phi = np.deg2rad(phi)
        llambda = np.deg2rad(llambda)

    In_ned = np.array([1,0,0])
    Ie_ned = np.array([0,1,0])
    Id_ned = np.array([0,0,1])

    Ix_ned = np.array([-np.cos(llambda) * np.sin(phi),-np.sin(llambda),-np.cos(llambda)*np.cos(phi)])
    Iy_ned = np.array([-np.sin(llambda) * np.sin(phi),np.cos(llambda),-np.sin(llambda)*np.cos(phi)])
    Iz_ned = np.array([np.cos(phi),0,-np.sin(phi)])

    C_ned2ecef = np.array([[np.dot(Ix_ned,In_ned),np.dot(Ix_ned,Ie_ned),np.dot(Ix_ned,Id_ned)],
                    [np.dot(Iy_ned,In_ned),np.dot(Iy_ned,Ie_ned),np.dot(Iy_ned,Id_ned)],
                    [np.dot(Iz_ned,In_ned),np.dot(Iz_ned,Ie_ned),np.dot(Iz_ned,Id_ned)]])

    return C_ned2ecef

def C_ecef2ned(phi,llambda,angtype='deg'):
    return np.linalg.inv(C_ned2ecef(phi,llambda,angtype))


def C_enu2ecef(phi,llambda,angtype='deg'):
    if angtype == 'deg':
        phi = np.deg2rad(phi)
        llambda = np.deg2rad(llambda)


    Ie_enu = np.array([1,0,0])
    In_enu = np.array([0,1,0])
    Iu_enu = np.array([0,0,1])

    Ix_enu = np.array([-np.sin(llambda),-np.cos(llambda) * np.sin(phi),np.cos(llambda)*np.cos(phi)])
    Iy_enu = np.array([np.cos(llambda),-np.sin(llambda)*np.sin(phi),np.sin(llambda)*np.cos(phi)])
    Iz_enu = np.array([0,np.cos(phi),np.sin(phi)])


    C_enu2ecef = np.array([[np.dot(Ix_enu,Ie_enu),np.dot(Ix_enu,In_enu),np.dot(Ix_enu,Iu_enu)],
                [np.dot(Iy_enu,Ie_enu),np.dot(Iy_enu,In_enu),np.dot(Iy_enu,Iu_enu)],
                [np.dot(Iz_enu,Ie_enu),np.dot(Iz_enu,In_enu),np.dot(Iz_enu,Iu_enu)]])

    return C_enu2ecef

def C_ecef2enu(phi,llambda,angtype='deg'):
    return C_enu2ecef(phi,llambda,angtype).T


def C_rpy2enu(roll,pitch,yaw,angtype='deg'):

    if angtype == 'deg':
        roll = np.deg2rad(roll)
        pitch = np.deg2rad(pitch)
        yaw = np.deg2rad(yaw)

    R = roll
    P = pitch
    Y = yaw

    Ir = np.array([np.sin(Y)*np.cos(P),np.cos(Y)*np.cos(P),np.sin(P)])
    Ip = np.array([np.cos(R) * np.cos(Y) + np.sin(R)*np.sin(Y)*np.sin(P),
                   -np.cos(R)*np.sin(Y)+np.sin(R)*np.cos(Y)*np.sin(P),
                   -np.sin(R)*np.cos(P)])
    Iy = np.array([-np.sin(R) *  np.cos(Y) + np.cos(R) * np.sin(Y) * np.sin(P) ,
                   np.sin(R) * np.sin(Y) + np.cos(R) * np.cos(Y) * np.sin(P) ,
                   - np.cos(R) * np.cos(P)    ])

#    Ie = np.array([ np.sin(Y) * np.cos(P) ,
#    np.cos(R) * np.cos(Y) + np.sin(R) * np.sin(Y) * np.sin(P) ,
#    -np.sin(R) * np.cos(Y) + np.cos(R) * np.sin(Y) * np.sin(P) ])
#
#    In = np.array([ np.cos(Y) * np.cos(P),
#     -np.cos(R) * np.sin(Y) + np.sin(R) * np.cos(Y) * np.sin(P),
#     np.sin(R) * np.sin(Y) + np.cos(R) * np.cos(Y) * np.sin(P) ])
#
#    Iu = np.array([ np.sin(P),
#     -np.sin(R) * np.cos(P),
#     -np.cos(R) * np.cos(P) ])

    C_rpy2enu = np.hstack((Ir[np.newaxis].T,Ip[np.newaxis].T,Iy[np.newaxis].T))

    return C_rpy2enu


def C_rpy2enu2(roll,pitch,yaw,angtype='deg'):

    if angtype == 'deg':
        roll = np.deg2rad(roll)
        pitch = np.deg2rad(pitch)
        yaw = np.deg2rad(yaw)

    R = roll
    P = pitch
    Y = yaw

    Sr = np.sin(R)
    Cr = np.cos(R)
    Sp = np.sin(P)
    Cp = np.cos(P)
    Sy = np.sin(Y)
    Cy = np.cos(Y)

    C_rpy2enu = np.array([[Sy*Cp,  Cr*Cy + Sr*Sy*Sp , - Sr*Cy + Cr*Sy*Sp],
                         [Cy*Cp , -Cr*Sy + Sr*Cy*Sp ,   Sr*Sy + Cr*Cy*Sp],
                         [Sp    , -Sr*Cp , -Cr * Cp]])

    return C_rpy2enu


def C_eci2rtn(P,V):
    """
    gives the transformation matrix between ECI ref. frame and RTN (aka RIC or RSW)

    Parameters
    ----------
    P : numpy.array
        3D vector, position of the object in ECI frame

    V : numpy.array
        3D vector, velocity of the object in ECI frame
        
    Returns
    -------
    C_eci2rtn : numpy.array
        transformation matrix
            
    Source
    ------
    "Coordinate Systems", ASEN 3200 1/24/06 George H. Born
    https://www.colorado.edu/ASEN/asen3200/handouts/Coordinate%20System.pdf
    """


    H = np.cross(P,V)

    R = P / np.linalg.norm(P)
    C = H / np.linalg.norm(H)
    I = np.cross(C,R)

    C_eci2rtn = np.column_stack((R,I,C))

    return C_eci2rtn


def C_rtn2rpy():
    C_rtn2rpy = np.array([[ 0.,  0  ,  -1.],
                          [ 1.,  0. ,   0.],
                          [ 0,  -1. ,   0.]])
    return C_rtn2rpy


def C_cep2itrs(xpole , ypole):
    """
    xpole and ypole are given in mas, miliarcsecond 
    
    Source
    ------
        Hofmann-Wellenhof,et al. 2008 - GNSS - p21
        Xu - 2007 - GPS - p17
        Xu - Orbits - p15
    """
    
    xpole2 = np.deg2rad(xpole * (1./3600.) * (10** -3))
    ypole2 = np.deg2rad(ypole * (1./3600.) * (10** -3))
    
    C_cep2itrs = np.array([[1.      , 0.     ,  xpole2],
                           [0.      , 1.     , -ypole2],
                           [-xpole2 , ypole2 , 1.     ]])
    
    return C_cep2itrs



def rot_quelconq(theta,Vx,Vy,Vz,angtype='deg'):

    N = np.linalg.norm([Vx,Vy,Vz])

    Vx = Vx / N
    Vy = Vy / N
    Vz = Vz / N

    if angtype == 'deg':
        theta = np.deg2rad(theta)

    c = np.cos(theta)
    s = np.sin(theta)

    R = np.array([[Vx ** 2 + (1-Vx**2)*c , Vx * Vy * (1 - c) - Vz * s , Vx * Vz * (1-c) +Vy * s],
                 [Vx * Vy * (1-c) + Vz * s , Vy**2 + (1-Vy**2)*c , Vy*Vz * (1-c) - Vx * s ],
                 [Vx * Vz * (1-c) - Vy * s , Vy * Vz * ( 1-c) +Vx * s , Vz**2 + (1-Vz**2) * c]])

    return R



def vector_RPY(V,out_deg=True,ad_hoc_mode=False):
    """
    Donne le "Roll Pitch Yaw" pour un Point/Vecteur
    (les coordonnées spheriques plus exactement, le roll sera tjrs nul)

    le ad_hoc_mode est la pour corriger un bug que je n'explique pas
    pour l'utilisation de ads_offset
    (cette methode a été trouvé par les scripts ENURPY2 et ENURPY3)
    """

    x = V[0]
    y = V[1]
    z = V[2]

    n = np.linalg.norm(V)

    roll  = 0
    pitch = np.arcsin(z/n)
    yaw   = np.arctan2(y,x)

    if out_deg:
        roll  = np.rad2deg(roll)
        pitch = np.rad2deg(pitch)
        yaw   = np.rad2deg(yaw)

    if ad_hoc_mode:
        roll  = -roll
        pitch = -pitch
        yaw   = -yaw

        roll , pitch , yaw = pitch , roll , yaw

    return roll ,  pitch , yaw

def add_offset(Direction , Delta , Point = None ,
               out_delta_enu = False ):
    """
    Un mode adhoc a du être implémenté dans la fonction élémentaire
    vector_RPY pour que ça puisse marcher correctement, reste à
    comprendre pourquoi ...

    out_delta_enu est prioritaire sur  les autres modes de sortie
    et renvoie le  delta

    par défaut on a le delta + le vecteur Direction
    et si Point est précisé : Point + Direction
    """

    rpy = vector_RPY(Direction,ad_hoc_mode=True)

    Crpy2enu = C_rpy2enu(*rpy)
    #Cenu2rpy = np.linalg.inv(C_rpy2enu(*rpy))

    Delta_enu = Crpy2enu.dot(C_enu2ned().dot(Delta))

    if out_delta_enu:
        print(out_delta_enu)

    if not np.isclose(np.linalg.norm(Delta) , np.linalg.norm(Delta_enu)):
        print("WARN : np.linalg.norm(Delta) != np.linalg.norm(Delta_enu)")
        print(np.linalg.norm(Delta) , np.linalg.norm(Delta_enu) , \
        np.linalg.norm(Delta) - np.linalg.norm(Delta_enu))


    if not type(Point) is None:
        outpoint = Point + Delta_enu
    else:
        outpoint = Direction + Delta_enu

    return outpoint




#                       _         _____                              _             
#     /\               | |       / ____|                            (_)            
#    /  \   _ __   __ _| | ___  | |     ___  _ ____   _____ _ __ ___ _  ___  _ __  
#   / /\ \ | '_ \ / _` | |/ _ \ | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \ 
#  / ____ \| | | | (_| | |  __/ | |___| (_) | | | \ V /  __/ |  \__ \ | (_) | | | |
# /_/    \_\_| |_|\__, |_|\___|  \_____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
#                  __/ |                                                           
#                 |___/                                                            

### Angle conversion
    
def dms2dec_num(deg,minn=0,sec=0):
    """
    Angle conversion
    
    Convert :
    Degree Minute Second `float` Angle => decimal Degree `float` Angle
        
    Parameters
    ----------
    deg & minn & sec : float
        degres, minutes and seconds of the input Angle

    Returns
    -------
    dd_float : float
        Decimal degree Angle
    """
    return deg + minn * (1./60.) +  sec * (1./3600.)


def deg_dec2dms(deg_in):
    """
    Angle conversion
    
    Convert :
    Decimal Degree => Degree Minute Seconds
        
    Parameters
    ----------
    deg_in : float or iterable of floats
        Decimal degrees

    Returns
    -------
    deg , minu , sec : numpy.array of float
        3 arrays for Degrees, Minutes, Seconds
    """
    deg              = np.floor(deg_in)
    decimal_part     = deg_in - deg 
    decimal_part_sec = decimal_part * 3600
    minu             = np.floor_divide(decimal_part_sec,60)
    sec              = decimal_part_sec - minu * 60
    sec              = np.round(sec,8)
    return deg , minu , sec



def dms2dec(dms_str , onlyDM=False):
    """   
    Angle conversion

    Convert :
    Degree Minute Second `string` Angle => decimal Degree `float` Angle
    
    can manage only DM angle in input  
    
    NB : Too Complicated .... must be simplified
    
    Parameters
    ----------
    dms_str : str
        string of DMS angle
        e.g.
        "2°20'35.09"E"

    onlyDM : bool
        True if the string is only a Degree Minute Angle 
        e.g.
        "40°52.0931'N"
        "28°31.4136'E"
                
    Returns
    -------
    dd_float : float
        Decimal degree Angle
    """
    
    if not onlyDM:
        print("WARN : DMS mode not well implemented yet !!! ")
    
    
    dms_str = dms_str.strip()

    if re.match('[swoSWO]', dms_str):
        sign = -1
    else:
        sign = 1

    # former regex
    lis = re.findall(r'\d+|\D+', dms_str)

    ipt = lis.index('.')
    decimal = lis[ipt-1] +  lis[ipt] + lis[ipt+1]
    lis[ipt-1] = decimal
    lis = lis[:ipt]

    for i,e in enumerate(lis):
        try:
            lis[i] = float(e)
        except:
            lis.remove(e)

    lis = [float(e) for e in lis]
    print(lis)

    deg = lis[0]
    minu = lis[1]
    if onlyDM:
        sec = 0.
    else:
        try:
            sec  = lis[2]
        except IndexError:
            print("ERR : did you forgot to activate the DM only mode ? ;) ?")
            return None

    dd_float = sign * (deg + minu / 60. + sec / 3600.)
    return dd_float

def arcsec2deg(arcsec_in):
    """
    Angle conversion
    
    Convert :
    Arcsecond => Degrees
    
    NB : Not really useful, should be avoided
    """
        
    return arcsec_in / 3600.

def deg2arcsec(deg_in):
    """
    Angle conversion
    
    Convert :
    Degrees => Arcsecond
    
    NB : Not really useful, should be avoided
    """
    return deg_in * 3600.


def angle2equivalent_earth_radius(angle_in,angtype='deg',
                                  earth_radius = 6371008.8):
    """
    Quick and simple function which gives the equivalent distance on a 
    Earth great circle of an angle
    
    angle can be : "deg", "rad", "mas"
    """
    earth_circum = earth_radius * 2 * np.pi
    if angtype == "deg":
        equiv_out = (angle_in * earth_circum) / 360.
    elif angtype == "rad":
        equiv_out = (angle_in * earth_circum) / (np.pi *2)
    elif angtype == "mas":
        equiv_out = (angle_in * 10**-3 * earth_circum) / (86400.)  
        
    return equiv_out
        
def anglesfromvects(xa,ya,xb,yb,angtype='deg'):
    """
    Determines the angle between the points A(xa,ya) and B(xb,yb)
    
    angle can be : "deg", "rad"
    """
    A = np.array([xa,ya])
    B = np.array([xb,yb])

    ps = np.inner(A,B)
    a = np.arccos( ps / (np.linalg.norm(A) * np.linalg.norm(B)) )

    if angtype == 'deg':
        return np.rad2deg(a)
    elif angtype == 'rad':
        return a
    else:
        print('ERR : angfromvects : mauvais angtype')

def angle_interpolation_quick(A,B,w):
    """
    Determine the interpolation between angle A & B
    by conversion to the cartesian space
    the parameter w € [0,1] define the interpoled angle C(w)
    where C(w=0) = A  &  C(w=1) = B

    Source
    ------
    https://stackoverflow.com/questions/2708476/rotation-interpolation
    """
    CS = (1-w)*np.cos(A) + w*np.cos(B)
    SN = (1-w)*np.sin(A) + w*np.sin(B)
    C = np.atan2(SN,CS)

    return C

def angle_from_3_pts(p1,p2,p3):
    """
    Gives angle between 3 points
    p3 is the vertex (sommet)
    
    Source
    ------
    http://www.les-mathematiques.net/phorum/read.php?8,596072,596231
    """

    p1,p2,p3 = [np.array(e) for e in (p1,p2,p3)]
    x1 , y1 = p1
    x2 , y2 = p2
    x3 , y3 = p3

    kos = ((x1 - x3)*(x2 - x3) + (y1 - y3)*(y2 - y3)) / (np.linalg.norm(p1 - p3) * np.linalg.norm(p2 - p3))

    return np.arccos(kos)


def cartesian2polar(x,y):
    """
    Coordinates conversion 
    
    cartesian => polar conversion
    
    Parameters
    ----------
    x , y : numpy.array of float
        cartesian coordinates

    Returns
    -------
    r , theta : float 
        polar coordinates (radius / angle in radians)
    """
    if cnv_gen.is_iterable(x):
        x = np.array(x)
    if cnv_gen.is_iterable(y):
        y = np.array(y)

    theta = np.arctan2(y,x)
    r     = np.sqrt(x** + y**2)
    return r , theta

def polar2cartesian(r,theta,ang='deg'):
    """
    Coordinates conversion 
    
    polar => cartesian conversion
    
    Parameters
    ----------
    r , theta : float or iterable of floats
        polar coordinates

    ang : string
        'deg' (degrees) or 'rad' (radian)

    Returns
    -------
    x , y : numpy.array of float
        cartesian coordinates
    """
    if cnv_gen.is_iterable(r):
        r = np.array(r)
    if cnv_gen.is_iterable(theta):
        theta = np.array(theta)

    if ang == 'deg':
        theta = np.deg2rad(theta)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x , y

#  _                     _                _                                                 
# | |                   | |              | |                                                
# | |     _____      __ | | _____   _____| |                                                
# | |    / _ \ \ /\ / / | |/ _ \ \ / / _ \ |                                                
# | |___| (_) \ V  V /  | |  __/\ V /  __/ |                                                
# |______\___/ \_/\_/   |_|\___| \_/ \___|_|_         __                  _   _             
#                                 | |      (_)       / _|                | | (_)            
#   __ _  ___  ___  _ __ ___   ___| |_ _ __ _  ___  | |_ _   _ _ __   ___| |_ _  ___  _ __  
#  / _` |/ _ \/ _ \| '_ ` _ \ / _ \ __| '__| |/ __| |  _| | | | '_ \ / __| __| |/ _ \| '_ \ 
# | (_| |  __/ (_) | | | | | |  __/ |_| |  | | (__  | | | |_| | | | | (__| |_| | (_) | | | |
#  \__, |\___|\___/|_| |_| |_|\___|\__|_|  |_|\___| |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|
#   __/ |                                                                                   
#  |___/       

### Low level geometric function

def dist(A,B):
    """
    Cartesian distance between 2 points A & B

    Parameters
    ----------
    A,B : numpy.array of float
        Description param1
        
    Returns
    -------
    D : float
        distance between A & B
    """
    A = np.array(A)
    B = np.array(B)
    return np.linalg.norm(A - B)

def dist_diff(A,B):
    """
    First derivative of cartesian distance between 2 points A & B

    Parameters
    ----------
    A, B : numpy.array of float
        "points", can be 2D or 3D vectors (list, np.array ...)

    Returns
    -------
    diffA : numpy.array of float 
        [ dD/dxa , dD/dya , dD/dza ]

    diffB : numpy.array of float
        [ dD/dxb , dD/dyb , dD/dzb ] = -diffA
    """

    dAB   = A-B
    dist  = scipy.linalg.norm(dAB)

    diffA =   dAB / dist
    diffB = - dAB / dist

    return diffA, diffB

def relative_orientation(x1,y1,x2,y2,trigo_orient = True):
    """
    Return the angle between a point 1 and a point 2
    (reference vector : to the right)
    """
    if trigo_orient:
        ang1 = np.mod(360 - np.rad2deg(np.arctan2((x2 - x1),(y2 - y1))),360)
    else:
        ang1 = np.mod(np.rad2deg(np.arctan2((x2 - x1),(y2 - y1))),360)
    return ang1

def barycenter(points_list_in):
    """
    Determines the barycenter of a list of points
    """
    points_arr = np.array(points_list_in)
    return np.mean(points_arr[:,-3:], axis=0)

def pythagore(a,b,c=0):
    """
    Computes Pythagore's formula
    """
    return np.sqrt(a**2 + b**2 + c**2)

def equilateral_triangle(side):
    """
    Gives coordinates of an equilateral triangle of a given side
    """
    hauteur_len = np.sqrt(side**2 - (side * .5) **2)
    A = np.array([0,(2./3.) * hauteur_len])
    B = np.array([side / 2. , - (1./3.) * hauteur_len])
    C = np.array([- side / 2. , - (1./3.) * hauteur_len])

    return A,B,C


def vincenty_full(point1, point2, miles=False,full=True,azimuth_in_deg=True):
    """
    
    Vincenty's formula (inverse method) to calculate the distance (in
    kilometers or miles) between two points on the surface of a spheroid. 
    Gives also the Azimuth between 2 points

    Parameters
    ----------
    point1, point2 : iterable of float
        coordinates of the points

    miles : bool
        kilometers if True
        
    full : bool
        Description param3

    azimuth_in_deg : bool
        azimut in Rad if False
                
    Returns
    -------
    s : float
        Distance betwwen the 2 points
    
    fwdAz,revAz : float
        Forward and Reverse Azimuth between the 2 points 
        
    Source
    ------
    https://github.com/maurycyp/vincenty/blob/master/vincenty/
    
    Examples
    --------
    >>> vincenty((0.0, 0.0), (0.0, 0.0))  # coincident points
    0.0
    >>> vincenty((0.0, 0.0), (0.0, 1.0))
    111.319491
    >>> vincenty((0.0, 0.0), (1.0, 0.0))
    110.574389
    >>> vincenty((0.0, 0.0), (0.5, 179.5))  # slow convergence
    19936.288579
    >>> vincenty((0.0, 0.0), (0.5, 179.7))  # failure to converge
    >>> boston = (42.3541165, -71.0693514)
    >>> newyork = (40.7791472, -73.9680804)
    >>> vincenty(boston, newyork)
    298.396057
    >>> vincenty(boston, newyork, miles=True)
    185.414657
    """
    
    # WGS 84
    a = 6378137  # meters
    f = 1 / 298.257223563
    b = 6356752.314245  # meters; b = (1 - f)a
    
    MILES_PER_KILOMETER = 0.621371
    
    MAX_ITERATIONS = 200
    CONVERGENCE_THRESHOLD = 1e-12  # .000,000,000,001

    # short-circuit coincident points
    if point1[0] == point2[0] and point1[1] == point2[1]:
        return 0.0

    import math 
    
    U1 = math.atan((1 - f) * math.tan(math.radians(point1[0])))
    U2 = math.atan((1 - f) * math.tan(math.radians(point2[0])))
    L = math.radians(point2[1] - point1[1])
    Lambda = L

    sinU1 = math.sin(U1)
    cosU1 = math.cos(U1)
    sinU2 = math.sin(U2)
    cosU2 = math.cos(U2)
    sinL  = math.sin(L)
    cosL  = math.cos(L)

    for iteration in range(MAX_ITERATIONS):
        sinLambda = math.sin(Lambda)
        cosLambda = math.cos(Lambda)
        sinSigma = math.sqrt((cosU2 * sinLambda) ** 2 +
                             (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ** 2)
        if sinSigma == 0:
            return 0.0  # coincident points
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
        sigma = math.atan2(sinSigma, cosSigma)
        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
        cosSqAlpha = 1 - sinAlpha ** 2
        try:
            cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha
        except ZeroDivisionError:
            cos2SigmaM = 0
        C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
        LambdaPrev = Lambda
        Lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma *
                                               (cos2SigmaM + C * cosSigma *
                                                (-1 + 2 * cos2SigmaM ** 2)))
        if abs(Lambda - LambdaPrev) < CONVERGENCE_THRESHOLD:
            break  # successful convergence
    else:
        return None  # failure to converge

    uSq = cosSqAlpha * (a ** 2 - b ** 2) / (b ** 2)
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma *
                 (-1 + 2 * cos2SigmaM ** 2) - B / 6 * cos2SigmaM *
                 (-3 + 4 * sinSigma ** 2) * (-3 + 4 * cos2SigmaM ** 2)))
    s = b * A * (sigma - deltaSigma)

    s /= 1000  # meters to kilometers
    if miles:
        s *= MILES_PER_KILOMETER  # kilometers to miles


    if not full:
        return s
    else:        
        fwdAz = np.arctan2(cosU2*sinL,  cosU1*sinU2-sinU1*cosU2*cosL)
        revAz = np.arctan2(cosU1*sinL, -sinU1*cosU2+cosU1*sinU2*cosL)
       
        if azimuth_in_deg:
            fwdAz = np.rad2deg(fwdAz)
            revAz = np.rad2deg(revAz)
       
        return s,fwdAz,revAz
    
def orthogonal_projection(Xa,Xb,Xv):
    """
    Project a point A on a line defined with a vector V and a point B
    
    Parameters
    ----------
    Xa : list/numpy.array of float
        Coordinates of A point, we want to project

    Xb : list/numpy.array of float
        Coordinates of B point, the origin point of the vector
        
    Xv : list/numpy.array of float
        Coordinates of the line direction vector
                
    Returns
    -------
    Xh : numpy.array of float
        Coordinates of H point, projection of A
    
    D : float
        Distance between A and H
        
    Note
    ----
    Misc. Notes

    Source
    ------
    https://fr.wikipedia.org/wiki/Projection_orthogonale
    """
    
    xa , ya = Xa
    xb , yb = Xb    
    xv , yv = Xv
    
    D =np.sqrt(xv**2 + yv **2) 
    BH = ((xa - xb)*xv + (ya - yb)*yv) / D
    
    xh = xb + (BH / D) * xv
    yh = yb + (BH / D) * yv
    
    Xh = np.array([xh,yh])
    return Xh , dist(Xa,Xh)


def line_maker(x1,y1,x2,y2,nbpts=10000):
    """
    Determine points of a line easily
    
    Parameters
    ----------
    x1,y1 : float
        Coordinates of the start point

    x2,y2 : float
        Coordinates of the end point
                
    Returns
    -------
    X,Y : numpy.array of float
        points of the line
    """

    X = np.linspace(x1,x2,nbpts)
    Y = np.linspace(y1,y2,nbpts)
    return X,Y
    

#  _____           _           _   _                _____           _        
# |  __ \         (_)         | | (_)              / ____|         | |       
# | |__) | __ ___  _  ___  ___| |_ _  ___  _ __   | |     __ _ _ __| |_ ___  
# |  ___/ '__/ _ \| |/ _ \/ __| __| |/ _ \| '_ \  | |    / _` | '__| __/ _ \ 
# | |   | | | (_) | |  __/ (__| |_| | (_) | | | | | |___| (_| | |  | || (_) |
# |_|   |_|  \___/| |\___|\___|\__|_|\___/|_| |_|  \_____\__,_|_|   \__\___/ 
#                _/ |                                                        
#               |__/        

def latitude_isometric(phi,e):
    """
    ALG0001 in NTG_71 (IGN Lambert Projection Tech document)
    """
    A = np.tan( (np.pi/4) + (phi/2) )
    B1 = 1 - e * np.sin(phi)
    B2 = 1 + e * np.sin(phi)
    
    L = np.log(A * ((B1/B2)**(.5 * e)))
    
    return L

def lambert_projection(long,lat,n,c,e,lc,Xs,Ys):
    """
    ALG0003 in NTG_71 (IGN Lambert Projection Tech document)
    """
    L = latitude_isometric(lat,e)
    X = Xs + c*np.exp(-n * L) * np.sin(n * (long - lc))
    Y = Ys - c*np.exp(-n * L) * np.cos(n * (long - lc))
    return X,Y

def lambert_secante_parameter(a , e , l0 , phi0 , phi1 , phi2 , X0 , Y0):
    """
    ALG0054 in NTG_71 (IGN Lambert Projection Tech document)
    """
    lc = l0
    A = np.log((wnorm(phi2,a,e**2) * np.cos(phi2)) / (wnorm(phi1,a,e**2) * np.cos(phi1)))
    B = latitude_isometric(phi1,e) - latitude_isometric(phi2,e)
    
    n = A/B
    
    C = ((wnorm(phi1,a,e**2) * np.cos(phi1))/n) * np.exp(n * latitude_isometric(phi1,e))
        
    if np.isclose(phi0,np.pi/2):
        Xs = X0
        Ys = Y0  
    else:
        Xs = X0
        Ys = Y0 + C * np.exp(n * latitude_isometric(phi0,e))
            
    return e , n , C , lc , Xs , Ys

def lambert_projection_CC_frontend(long,lat,NZ = 93):
    """
    Project WGS84/ITRF Latitude/Longitude to Lambert Conique Conforme

    Parameters
    ----------
    long,lat : float
        WGS84/ITRF Longitude/Latitude
        
    NZ : int
        Lambert93, NZ = 0 or = 93
                
    Returns
    -------
    X,Y : float
        Projected coordinates
    
    Source
    ------
    https://fr.wikipedia.org/wiki/Projection_conique_conforme_de_Lambert   
    geodesie.ign.fr/contenu/fichiers/documentation/algorithmes/notice/NTG_71.pdf
    https://geodesie.ign.fr/contenu/fichiers/documentation/rgf93/Lambert-93.pdf
 
    """
    if NZ in (0,93):
        phi0 = np.deg2rad(46.5)
        phi1 = np.deg2rad(44)
        phi2 = np.deg2rad(49)
        l0 = np.deg2rad(3)
        X0 = 700000.
        Y0 = 6600000.
        a = 6378137.
        e = 0.0818191910435
    else:        
        phi0 = np.deg2rad(41 + NZ)
        phi1 = np.deg2rad(phi0 - 0.75)
        phi2 = np.deg2rad(phi0 + 0.75)
        l0 = np.deg2rad(3)
        X0 = 1700000.
        Y0 = (NZ * 1000000.) + 200000.
        a = 6378137.
        e = 0.0818191910435
    
    e , n , c , lc , Xs , Ys = lambert_secante_parameter(a, e, l0, phi0, phi1, phi2, X0, Y0)
    X,Y = lambert_projection(np.deg2rad(long),np.deg2rad(lat),n,c,e,lc,Xs,Ys)
    
    return X,Y


#  _    _ _       _       _                    _    _____                _      _   _        ______   _       
# | |  | (_)     | |     | |                  | |  / ____|              | |    | | (_)      |  ____| | |      
# | |__| |_  __ _| |__   | |     _____   _____| | | |  __  ___  ___   __| | ___| |_ _  ___  | |__ ___| |_ ___ 
# |  __  | |/ _` | '_ \  | |    / _ \ \ / / _ \ | | | |_ |/ _ \/ _ \ / _` |/ _ \ __| |/ __| |  __/ __| __/ __|
# | |  | | | (_| | | | | | |___|  __/\ V /  __/ | | |__| |  __/ (_) | (_| |  __/ |_| | (__  | | | (__| |_\__ \
# |_|  |_|_|\__, |_| |_| |______\___| \_/ \___|_|  \_____|\___|\___/ \__,_|\___|\__|_|\___| |_|  \___|\__|___/
#            __/ |                                                                                            
#           |___/                                                                         
    
### High level geodetic function

def itrf_speed_calc(x0,y0,z0,t0,vx,vy,vz,t):
    """
    Args :
        x0,y0,z0 (floats) : coordinates at the reference epoch (m)

        t0 (float) : reference epoch (decimal year)

        vx,vy,vz (floats) : speed of the point (m/yr)

        t (float) : output epoch
    Returns :
        xout,yout,zout : coordinates of the point @ the ref. epoch (m)
    """

    xout = x0 + vx * ( t - t0 )
    yout = y0 + vy * ( t - t0 )
    zout = z0 + vz * ( t - t0 )

    return xout,yout,zout


def itrf_psd_fundamuntal_formula(t,A_l,t_l,tau_l,A_e,t_e,tau_e):
    """
    http://itrf.ensg.ign.fr/ITRF_solutions/2014/doc/ITRF2014-PSD-model-eqs-IGN.pdf
    """
    
    dL = A_l * np.log(1 + (t - t_l)/ tau_l) + A_e * (1 + (t - t_e)/ tau_e)


def calc_pos_speed_itrf(x0,y0,z0,t0,vx,vy,vz,t):
    """
    just a wrapper of itrf_speed_calc
    for legacy reasons
    """
    return itrf_speed_calc(x0,y0,z0,t0,vx,vy,vz,t)


def helmert_trans(Xa,params='itrf2008_2_etrf2000',invert=True,workepoc=2009.):
    """
    NB 1 : http://etrs89.ensg.ign.fr/memo-V8.pdf
    NB 2 :
    Transformation inverse : * -1 pour les paramètres
    cf https://en.wikipedia.org/wiki/Helmert_transformation

    optimisé pour RGF93 => ITRF2008 (d'ou le invert = True & workepoc = 2009.)

    NB3 : Attention losque l'on compare avec
          le convertisseur EUREF avec une vitesse
          parce que elle aussi est modifiée dans la conversion ETRS => ITRS.
          Conclusion : le faire en 2 étapes ETRS = > ITRS dans la même epoc
                                            ITRS epoc 1 => ITRS epoc 2

    if manual
    then params is a tuple
    params = (t1,t2,t3,dab,r1,r2,r3)

    NE MARCHE PAS PARCE BESOIN DU RATE(TAUX) DES PARAMS D'HELMERT !!!!!!
    (160923)
    """
    if invert:
        inver = -1.
    else:
        inver = 1.

    mas2rad = 0.0000000048481368
    mas2rad = 4.8481368111e-6 * 1e-3

    if params == 'itrf2008_2_etrf2000':

        t1rate = .1 * 10**-3
        t2rate = .1 * 10**-3
        t3rate = -1.8 * 10**-3
        dabrate = .08 * 10**-9
        r1rate =  .081 * mas2rad
        r2rate =  .490 * mas2rad
        r3rate = -.792 * mas2rad

        t1rate = .1 * 10**-3
        t2rate = .1 * 10**-3
        t3rate = -1.8 * 10**-3
        dabrate = .08 * 10**-9
        r1rate =  .081 * mas2rad
        r2rate =  .490 * mas2rad
        r3rate = -.792 * mas2rad


        t1  =(52.1  * 10**-3   + t1rate  * ( workepoc - 2000.)) *inver
        t2  =(49.3  * 10**-3   + t2rate  * ( workepoc - 2000.)) *inver
        t3  =(-58.5 * 10**-3   + t3rate  * ( workepoc - 2000.)) *inver
        dab =( 1.34 * 10**-9   + dabrate * ( workepoc - 2000.)) *inver
        r1  =( 0.891 * mas2rad + r1rate  * ( workepoc - 2000.)) *inver
        r2  =( 5.39  * mas2rad + r2rate  * ( workepoc - 2000.)) *inver
        r3  =( -8.712* mas2rad + r3rate  * ( workepoc - 2000.)) *inver


    elif params == 'itrf2000_2_etrf2000':
        t1  =54.0  * 10**-3   *inver
        t2  =51.0  * 10**-3   *inver
        t3  =-48.0 * 10**-3   *inver
        dab = 0.0  * 10**-9   *inver
        r1  = 0.891 * mas2rad *inver
        r2  = 5.390 * mas2rad *inver
        r3  = -8.712* mas2rad *inver




    R = np.matrix([[dab,-r3,r2],
                   [r3,dab,-r1],
                   [-r2,r1,dab]])

    Xb = Xa + np.matrix([t1,t2,t3]) + np.dot(R,Xa)
    Xb = np.squeeze(np.array(Xb))

    return Xb


def helmert_trans_estim_matrixs_maker(X1 , X2):
    """
    internal function for helmert_trans_estim
    """
    x1 , y1 , z1 = X1
    x2 , y2 , z2 = X2
    
    
    block_1 = np.eye(3)
    
    block_2 = np.array([[ 0. , -z1,  y1, x1],
                       [ z1,   0., -x1, y1],
                       [-y1,  x1,   0., z1]])
    
    l = X2 - X1
    A = np.hstack((block_1 , block_2))

    return l,A


def helmert_trans_estim(X1list , X2list):
    """
    estimates 7 parameters of a 3D Helmert transformation between a set of points
    X1 and a set of points X2 (compute transformation X1 => X2)
    
    Parameters
    ----------
    
    X1list & X2list : list of N (x,y,z) points ,
        or an numpy array of shape (3, N)
    
    Returns
    -------
    Res :
        7 Helmert params. : x,y,z translations, x,y,z rotations, scale
    A :
        Design matrix    
    l :
        X2 - X1
        
    Source
    ------
    https://elib.uni-stuttgart.de/bitstream/11682/9661/1/BscThesis_GaoYueqing.pdf
    """

    l_stk = []
    A_stk = []
    
    
    for X1 , X2 in zip(X1list , X2list):
        lmono , Amono = helmert_trans_estim_matrixs_maker(X1,X2)
            
        l_stk.append(lmono)
        A_stk.append(Amono)
        
    A = np.vstack(A_stk)
    l = np.hstack(l_stk)
    
    Res = scipy.linalg.inv((A.T).dot(A)).dot(A.T).dot(l)
    
    return Res , A , l





#### ASTRONOMY FUNCTION
def semi_major_axis_from_mean_motion(n):
    """
    source : https://space.stackexchange.com/questions/18289/how-to-get-semi-major-axis-from-tle
    """
    mu = 3.9860044189 * 10**14
    a  = (mu**(1./3.)) / ((2*n*np.pi/86400)**(2./3.))
    return a    

