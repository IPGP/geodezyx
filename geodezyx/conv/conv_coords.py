#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:11:12 2020

@author: psakic

This sub-module of geodezyx.conv deals with coordinate conversion.

it can be imported directly with:
from geodezyx import conv

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""


#### Import the logger
import logging

########## BEGIN IMPORT ##########
#### External modules
import numpy as np
import scipy
from pyorbital import astronomy

#### geodeZYX modules
from geodezyx import utils
from geodezyx.conv import conv_rotation_matrices as rotmat

# import re
log = logging.getLogger(__name__)

#### Import star style
# from geodezyx import *                   # Import the GeodeZYX modules
# from geodezyx.externlib import *         # Import the external modules
# from geodezyx.megalib.megalib import *   # Import the legacy modules names

##########  END IMPORT  ##########
#### Coordinates conversion
    
#   _____                    _ _             _               _____                              _             
#  / ____|                  | (_)           | |             / ____|                            (_)            
# | |     ___   ___  _ __ __| |_ _ __   __ _| |_ ___  ___  | |     ___  _ ____   _____ _ __ ___ _  ___  _ __  
# | |    / _ \ / _ \| '__/ _` | | '_ \ / _` | __/ _ \/ __| | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \ 
# | |___| (_) | (_) | | | (_| | | | | | (_| | ||  __/\__ \ | |___| (_) | | | \ V /  __/ |  \__ \ | (_) | | | |
#  \_____\___/ \___/|_|  \__,_|_|_| |_|\__,_|\__\___||___/  \_____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
#                                                                                                             
#                                                                                                             


def vector_separator(ABC):
    """
    Split a Nx3 Array/DataFrame in three separated 1-component vectors
    To simplify the usage of the conversion functions 
    (which take single component vectors as input)

    Parameters
    ----------
    ABC : Array or DataFrame
        Nx3 XYZ, ENU... array.

    Returns
    -------
    A : Array
        1st component.
    B : Array
        2nd component.
    C : Array
        3rd component.

    """
    ABC = np.array(ABC)
    return ABC[:,0],ABC[:,1],ABC[:,2]



def wnorm(phi,a=6378137.,e2=0.00669438003):
    """
    Compute the Ellipsoid "Grande Normale"
    
    References
    ----------
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
    
    References
    ----------
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
        
    References
    ----------
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
        
    References
    ----------
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
    
    **use XYZ2ENU_2 in priority**
    
    (XYZ2ENU is a core function for XYZ2ENU_2)

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
        
    References
    ----------
    https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
    
    Note
    ----
    This recursive fuction should be improved
    """
    
    ## Case one ref point per dXYZ
    if utils.is_iterable(lat0):
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
        
    References
    ----------
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
    return the topocentric ENU position around its mean position
    
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


def ENU2XYZ(E,N,U,x0,y0,z0,velocity_mode=False):
    """
    Coordinates conversion

    ENU Topocentic => XYZ ECEF Geocentric 
    
    Parameters
    ----------
    E,N,U : numpy.array of floats
        cartesian coordinates E,N,U in meters
        
    x0,y0,z0 : floats or numpy.array of floats
        coordinate of the topocentric origin point in the geocentric frame   
        if they are iterable arrays (of the same size as E,N,U) 
        a different x0,y0,z0 will be applied for each E,N,U element
    velocity_mode : bool
        For the Velocity mode, coordinates of the topocentric origin point 
        are NOT added at the end of the conversion
        Default is False

    Returns
    -------
    X,Y,Z : numpy.array of floats
        ECEF X,Y,Z Component (m)
        
    References
    ----------
    https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
    """
    
    if utils.is_iterable(E):
        if not utils.is_iterable(x0): ### The input x0,y0,z0 is a single point 
            n = len(E)
            X0,Y0,Z0 = [x0]*n,[y0]*n,[z0]*n
        else:                         ### The input x0,y0,z0 is a vector = different x0,y0,z0 for each point
            X0,Y0,Z0 = x0,y0,z0
                        
        Xlist , Ylist , Zlist = [] , [] , []
        
        for e,n,u,x0use,y0use,z0use in zip(E,N,U,X0,Y0,Z0):
            x,y,z = ENU2XYZ(e,n,u,x0use,y0use,z0use,
                            velocity_mode=velocity_mode)
            Xlist.append(x)
            Ylist.append(y)
            Zlist.append(z)
        return np.array(Xlist) , np.array(Ylist) , np.array(Zlist)
    
    else:
        fr,lr,hr = XYZ2GEO(x0,y0,z0)
        f0 = np.deg2rad(fr)
        l0 = np.deg2rad(lr)
    
        R = np.array([[ -np.sin(l0)            , np.cos(l0)             , 0          ],
                      [ -np.sin(f0)*np.cos(l0) , -np.sin(f0)*np.sin(l0) , np.cos(f0) ],
                      [  np.cos(f0)*np.cos(l0) , np.cos(f0)*np.sin(l0)  , np.sin(f0)]])
    
        R3 = R.T
        R3 = scipy.linalg.inv(R)
    
        ENU = np.vstack((E,N,U))
    
        xyz = np.dot(R3,ENU)
    
        if velocity_mode:
            X = float(xyz[0])
            Y = float(xyz[1])
            Z = float(xyz[2])
        else:
            X = float(xyz[0]) + x0
            Y = float(xyz[1]) + y0
            Z = float(xyz[2]) + z0
        
        return X,Y,Z


def ENU2XYZ_legacy(E,N,U,Xref,Yref,Zref):
    """
    KEPT FOR LEGACY REASONS, use ENU2XYZ

    diffère de ENU2XYZ pour d'obscure raisons, à investiguer !!!
    est laissé pour des scripts de conversion de GINS (170119)
    
    this fct compute the dXYZ and not the final XYZ

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

    dX = float(xyz[0])
    dY = float(xyz[1])
    dZ = float(xyz[2])

    return dX,dY,dZ

def GEO2XYZ_vector(FLH,angle='deg',
                   a=6378137.,e2=0.00669438003):
    
        
    #### if Nx3 array => 3xN array
    FLH = utils.transpose_vector_array(FLH)
    
    X,Y,Z = GEO2XYZ(FLH[0],
                    FLH[1],
                    FLH[2],
                    angle=angle,
                    a=a,
                    e2=e2)
    XYZ = np.column_stack((X,Y,Z))
    
    return XYZ
    
def XYZ2ENU_vector(XYZ,xyz0):
    
    XYZ = utils.transpose_vector_array(XYZ)
    
    E,N,U = XYZ2ENU_2(XYZ[0], XYZ[1], XYZ[2],
                      xyz0[0], xyz0[1], xyz0[2])
    ENU = np.column_stack((E,N,U))
    
    return ENU
    
def XYZ2GEO_vector(XYZ,outdeg=True,
            A=6378137.,E2=0.00669438003):
    
    XYZ = utils.transpose_vector_array(XYZ)

    F,L,H = XYZ2GEO(XYZ[0], XYZ[1], XYZ[2],
                    outdeg=outdeg,
                    A=A,
                    E2=E2)
            
    FLH = np.column_stack((F,L,H))
    
    return FLH

    
def ENU2XYZ_vector(ENU,xyz_ref):
    
    ENU = utils.transpose_vector_array(ENU)
    
    X,Y,Z = ENU2XYZ(ENU[0],
            ENU[1],
            ENU[2],
            xyz_ref[0],
            xyz_ref[1],
            xyz_ref[2])
    
    XYZ = np.column_stack((X,Y,Z))
    
    return XYZ



def sFLH2sXYZ(F,L,H,sF,sL,sH,ang='deg'):
    """
    Convert standard deviation
    Geographic FLH => Cartesian ECEF XYZ
    
    WARNING
    -------
    Inputs values are assumed as uncorrelated, which is not accurate
    Must be improved
    
    References
    ----------
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
    
    References
    ----------
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
    
    References
    ----------
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
    Inputs values are now assumed as correlated (241105)
    
    References
    ----------
    Linear Algebra, Geodesy, and GPS p332
    """

    #conversion "batarde" en assumant par défaut
    #que XYZ ne sont pas corrélé (pas bien ...)
    #Linear Algebra, Geodesy, and GPS p332

    SIGMAxyz = np.array([[sX**2,sXY,sXZ],
                         [sXY,sY**2,sYZ],
                         [sXZ,sYZ,sZ**2]])

    F,L,H = XYZ2GEO(X,Y,Z)

    # old and bad rotation matrix (bofore 20241104)
    # C = rotmat.C_ecef2enu(F,L,angtype='deg')
    # new and good rotation matrix (after 20241104)
    C = rotmat.C_ecef2enu_sigma(F,L,angtype='deg')

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
        use the theoretical matrix composition, but wrong 
        ans. empirically, only for debug !!

    Returns
    -------
    Cout : conversion of P in RTN or RPY ref. frame

    References
    ----------
    "Coordinate Systems", ASEN 3200 1/24/06 George H. Born
    """

    C_eci2rtn_mat  = rotmat.C_eci2rtn(P,V)

    if not out_rpy:
        TransMat = C_eci2rtn_mat
    else:

        if not rpy_theo_mode:
            # Pour de très obscures raisons la composition est inversée
            # par rapport à l'ordre standard ... (241017)
            TransMat = np.dot(rotmat.C_rtn2rpy().T,C_eci2rtn_mat.T)
        else:
            log.warning("using the theoretical mode for RPY conversion, UNSTABLE & WRONG !")
            TransMat = np.dot(rotmat.C_rtn2rpy(),C_eci2rtn_mat)

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
        Earth Centered Inertial coordinates. will have shape (N,3)    

    Note
    ----
    Requires pyorbital module 
    
     [X]    [C -S 0][X]
     [Y]  = [S  C 0][Y]
     [Z]eci [0  0 1][Z]ecf

     C and S are cos() and sin() of gmst (Greenwich Meridian Sideral Time)

    References
    ----------
    http://ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
    Inspired from satellite-js (https://github.com/shashwatak/satellite-js)
    """

    # XYZ and utc_time must have the same shape
    #if not xyz.shape[:-1] == utc_times.shape:
    #    raise ValueError("shape mismatch for XYZ and utc_times (got {} and {})".format(xyz.shape[:-1],utc_times.shape))

    #    gmst = -1 * astronomy.gmst(utc_times) # EDIT 180430 : Why this -1 ??? removed because wrong ! ...
    gmst = 1 * astronomy.gmst(utc_times)

    eci = xyz.copy()
    eci[:,0] = xyz[:,0]*np.cos(gmst) - xyz[:,1]*np.sin(gmst)
    eci[:,1] = xyz[:,0]*np.sin(gmst) + xyz[:,1]*np.cos(gmst)
    return eci



def ECI2ECEF(xyz,utc_times):
    """ 
    Convert ECI (Earth Centered Inertial) positions to ECEF (Earth Centered Earth Fixed)
    positions

    Parameters
    ----------
    xyz : numpy.array of floats
        XYZ are cartesian positions in Earth Centered Inertial. Should have shape (N,3)
    
    utc_times : numpy.array of floats
        UTC_times are UTC timestamps, as datetime objects. Sould have shape (N)
    
    Returns
    -------
    ecef : numpy.array of floats
        Earth Centered Earth Fixed coordinates. will have shape (N,3)        

    Note
    ----
    Requires pyorbital module
 
    
     [X]          ([C -S 0])[X]
     [Y]     = inv([S  C 0])[Y]
     [Z]ecef      ([0  0 1])[Z]eci


    Empirically:
     [X]       [ C  S  0][X]
     [Y]     = [-S  C  0][Y]
     [Z]ecef   [ 0  0  1][Z]eci



     C and S are cos() and sin() of gmst (Greenwich Meridian Sideral Time)

    References
    ----------
    http://ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
    Inspired from satellite-js (https://github.com/shashwatak/satellite-js)
    
    
    Note
    ----
    Quick mode of the reverse fct, can be improved
    """
    # XYZ and utc_time must have the same shape
    #if not xyz.shape[:-1] == utc_times.shape:
    #    raise ValueError("shape mismatch for XYZ and utc_times (got {} and {})".format(xyz.shape[:-1],utc_times.shape))

    #    gmst = -1 * astronomy.gmst(utc_times) 
    # EDIT 180430 : Why this -1 ??? removed because wrong ! ...
    gmst = 1 * astronomy.gmst(utc_times)
    
    ecef = xyz.copy()
    ecef[:,0] = + xyz[:,0]*np.cos(gmst) + xyz[:,1]*np.sin(gmst)
    ecef[:,1] = - xyz[:,0]*np.sin(gmst) + xyz[:,1]*np.cos(gmst)
    
    return ecef
