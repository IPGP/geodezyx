# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 13:55:33 2019

@author: psakic
"""

import copy
import datetime as dt

#### Import the logger
import logging
import os

########## BEGIN IMPORT ##########
#### External modules
from collections import Counter

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from matplotlib.widgets import MultiCursor

#### geodeZYX modules
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import reffram
from geodezyx import stats
from geodezyx import time_series
from geodezyx import utils

log = logging.getLogger("geodezyx")
log.setLevel(logging.INFO)

##########  END IMPORT  ##########


class Point:

    def __init__(
        self,
        A=0.0,
        B=0.0,
        C=0.0,
        T=0.0,
        initype="XYZ",
        sA=0.0,
        sB=0.0,
        sC=0.0,
        name="noname",
        anex=None,
    ):
        """
        Initialize a Point Object by importing the coordinate component

        Parameters
        ----------
        A : float
            X, F (latitude), E (depends on initype).
        B : float
            Y, L (longitude), N (depends on initype).
        C : float
            Z, H (hight), U (depends on initype).
        T : float, optional
            Epoch of the measure. The default is 0..
        initype : str, optional
            The inital coordinates type. The default is 'XYZ'.
            The others are 'FLH', 'ENU' and 'NED'
        sA : float, optional
            Sigma (standard deviation) of A component. The default is 0.
        sB : float, optional
            Sigma (standard deviation) of B component. The default is 0.
        sC : float, optional
            Sigma (standard deviation) of C component. The default is 0.
        name : str, optional
            Flexible name for the Point identification. The default is 'noname'.
        anex : dict, optional
            Additional data. The default is None. See Note

        Notes
        -----

        A dictionary called anex is also initialized to allow a
        versatile storage of a variety of data.

        Examples of dictionary keys:
            - RMS: average RMS (for gipsy)
            - sdAB, sdBC, sdAC: the variances between A,B,C (for rtklib)
            - sdXY, sdXZ, sdYZ: the variances between XYZ (pbo.pos)
            - Vx, Vy, Vz, sVx, sVy, sVz: velocity of the point (EPOS coordinates)
        """

        self.Tset(T)
        self.ENUset()
        self.name = name
        if not anex:
            self.anex = dict()
        else:
            self.anex = anex

        # le dico "anex" permet de stocker de manière polyvalente des données diverses
        # On trouvera (LISTE SE DEVANT ETRE LA PLUS EXHAUSTIVE POSSIBLE )
        #
        # RMS : moyenne RMS (pour les gipsy bosser)
        # sdAB , sdBC , sdAC : les variances entre A,B,C (pour les rtklib)
        # sdXY , sdXZ , sdYZ : les variances entre XYZ (pbo.pos)
        # Vx , Vy , Vz , sVx , sVy , sVz : velocity of the point (EPOS coordinates)

        if initype == "XYZ":
            self.XYZset(A, B, C, sA, sB, sC)
        elif initype == "FLH":  # On travaille en degres decimaux
            self.FLHset(A, B, C, sA, sB, sC)

        elif initype == "ENU":
            self.ENUset(A, B, C, sA, sB, sC)

        elif initype == "NED":
            self.NEDset(A, B, C, sA, sB, sC)

        elif initype == "UTM":
            self.UTMset(A, B, C, sA, sB, sC)
        else:
            log.error("wrong initype")

    def __call__(self):

        if self.initype == "XYZ":
            return self.X, self.Y, self.Z, self.Tdt, self.T
        elif self.initype == "FLH":
            return self.F, self.L, self.H, self.Tdt, self.T
        elif self.initype == "ENU":
            return self.E, self.N, self.U, self.Tdt, self.T
        elif self.initype == "NED":
            return self.N, self.E, self.D, self.Tdt, self.T
        elif self.initype == "UTM":
            return self.Eutm, self.Nutm, self.Uutm, self.Tdt, self.T
        else:
            log.error("wrong initype")

    def __repr__(self):
        if not hasattr(self, "X"):
            return "{},{},{},{},{}".format(self.E, self.N, self.U, self.Tdt, self.T)
        else:
            return "{},{},{},{},{}".format(self.X, self.Y, self.Z, self.Tdt, self.T)

    def XYZset(self, X=0.0, Y=0.0, Z=0.0, sX=0.0, sY=0.0, sZ=0.0):
        """
        Set Cartesian XYZ coordinates for the Point.

        Parameters
        ----------
        X : float, optional
            X Cartesian coordinate. The default is 0.
        Y : float, optional
            Y Cartesian coordinate. The default is 0.
        Z : float, optional
            Z Cartesian coordinate. The default is 0.
        sX : float, optional
            Sigma (standard deviation) of X. The default is 0.
        sY : float, optional
            Sigma (standard deviation) of Y. The default is 0.
        sZ : float, optional
            Sigma (standard deviation) of Z. The default is 0.

        Returns
        -------
        None
        """
        self.X = X
        self.Y = Y
        self.Z = Z
        self.sX = sX
        self.sY = sY
        self.sZ = sZ

        self.initype = "XYZ"
        self.F, self.L, self.H = conv.xyz2geo(self.X, self.Y, self.Z)

    def FLHset(self, F=0.0, L=0.0, H=0.0, sF=0.0, sL=0.0, sH=0.0):
        """
        Set geodetic coordinates (F=latitude, L=longitude, H=height) for the Point.

        Parameters
        ----------
        F : float, optional
            Latitude in decimal degrees. The default is 0.
        L : float, optional
            Longitude in decimal degrees. The default is 0.
        H : float, optional
            Height (altitude) in meters. The default is 0.
        sF : float, optional
            Sigma (standard deviation) of F. The default is 0.
        sL : float, optional
            Sigma (standard deviation) of L. The default is 0.
        sH : float, optional
            Sigma (standard deviation) of H. The default is 0.

        Returns
        -------
        None
        """
        self.F = F
        self.L = L
        self.H = H
        self.sF = sF
        self.sL = sL
        self.sH = sH

        self.initype = "FLH"
        self.X, self.Y, self.Z = conv.geo2xyz(self.F, self.L, self.H)
        self.sX, self.sY, self.sZ = conv.sigma_geo2xyz(F, L, H, sF, sL, sH)

    def ENUset(self, E=np.nan, N=np.nan, U=np.nan, sE=np.nan, sN=np.nan, sU=np.nan):
        """
        Set East-North-Up local topocentric coordinates for the Point.

        Parameters
        ----------
        E : float, optional
            East component in meters. The default is NaN.
        N : float, optional
            North component in meters. The default is NaN.
        U : float, optional
            Up component in meters. The default is NaN.
        sE : float, optional
            Sigma (standard deviation) of E. The default is NaN.
        sN : float, optional
            Sigma (standard deviation) of N. The default is NaN.
        sU : float, optional
            Sigma (standard deviation) of U. The default is NaN.

        Returns
        -------
        None
        """
        self.E = E
        self.N = N
        self.U = U
        self.sE = sE
        self.sN = sN
        self.sU = sU

        self.initype = "ENU"
        self._flh_computed = False

    def NEDset(self, N=np.nan, E=np.nan, D=np.nan, sN=np.nan, sE=np.nan, sD=np.nan):
        """
        Set North-East-Down local topocentric coordinates for the Point.

        Parameters
        ----------
        N : float, optional
            North component in meters. The default is NaN.
        E : float, optional
            East component in meters. The default is NaN.
        D : float, optional
            Down component in meters. The default is NaN.
        sN : float, optional
            Sigma (standard deviation) of N. The default is NaN.
        sE : float, optional
            Sigma (standard deviation) of E. The default is NaN.
        sD : float, optional
            Sigma (standard deviation) of D. The default is NaN.

        Returns
        -------
        None
        """
        self.N = N
        self.E = E
        self.D = D
        self.sN = sN
        self.sE = sE
        self.sD = sD

        self.initype = "NED"
        self._flh_computed = False

    def UTMset(
        self,
        Eutm=np.nan,
        Nutm=np.nan,
        Uutm=np.nan,
        sEutm=np.nan,
        sNutm=np.nan,
        sUutm=np.nan,
    ):
        """
        Set UTM projected coordinates for the Point.

        Parameters
        ----------
        Eutm : float, optional
            UTM Easting coordinate in meters. The default is NaN.
        Nutm : float, optional
            UTM Northing coordinate in meters. The default is NaN.
        Uutm : float, optional
            Height (altitude) in meters. The default is NaN.
        sEutm : float, optional
            Sigma (standard deviation) of Eutm. The default is NaN.
        sNutm : float, optional
            Sigma (standard deviation) of Nutm. The default is NaN.
        sUutm : float, optional
            Sigma (standard deviation) of Uutm. The default is NaN.

        Returns
        -------
        None
        """
        self.Eutm = Eutm
        self.Nutm = Nutm
        self.Uutm = Uutm
        self.sEutm = sEutm
        self.sNutm = sNutm
        self.sUutm = sUutm

        self.initype = "UTM"
        self._flh_computed = False

    def add_offset(self, dA, dB, dC, coortype="ENU"):
        """
        Add an offset to the Point coordinates.

        Parameters
        ----------
        dA : float
            Offset for the A component (depends on coortype).
        dB : float
            Offset for the B component (depends on coortype).
        dC : float
            Offset for the C component (depends on coortype).
        coortype : str, optional
            Coordinate type ('ENU', 'XYZ', 'FLH', 'NED', 'UTM'). The default is 'ENU'.

        Returns
        -------
        None
        """
        temp = time_series.add_offset_point(self, dA, dB, dC, coortype=coortype)
        self.__dict__ = temp.__dict__

    def Tset(self, T=0):
        """
        Set the time/epoch of the Point.

        Parameters
        ----------
        T : float or datetime.datetime, optional
            Time value. If 0, defaults to current time. Can be POSIX timestamp
            or datetime object. The default is 0.

        Returns
        -------
        None

        Notes
        -----
        The time is stored in both POSIX timestamp format (self.T) and
        datetime object format (self.Tdt).
        """
        if T == 0:
            T = dt.datetime.now()

        if type(T) == dt.datetime:
            self.Tdt = T
            self.T = conv.dt2posix(T)

        else:
            self.T = float(T)
            self.Tdt = conv.posix2dt(float(T))

    def ENUcalc_pt(self, refENU):
        """
        Calculate East-North-Up (ENU) coordinates relative to a reference point.

        Parameters
        ----------
        refENU : Point
            Reference point for ENU coordinate calculation.

        Returns
        -------
        None

        Notes
        -----
        Updates the E, N, U attributes and calculates sigmas if available.
        """

        dX = self.X - refENU.X
        dY = self.Y - refENU.Y
        dZ = self.Z - refENU.Z

        Etmp, Ntmp, Utmp = conv.xyz2enu(
            self.X, self.Y, self.Z, refENU.X, refENU.Y, refENU.Z
        )

        self.E, self.N, self.U = Etmp[0], Ntmp[0], Utmp[0]

        if self.initype == "FLH" and hasattr(self, "s_f"):
            if not np.isnan(self.sF):
                self.sE, self.sN, self.sU = conv.sigma_geo2enu(
                    self.F, self.L, self.H, self.sF, self.sL, self.sH
                )
        elif self.initype == "XYZ" and hasattr(self, "sX"):
            if np.isnan(self.sX):

                return

            if "sdXY" in self.anex.keys():
                sXY = self.anex["sdXY"]
                sXZ = self.anex["sdXZ"]
                sYZ = self.anex["sdYZ"]
            else:
                sXY, sXZ, sYZ = 0, 0, 0

            self.sE, self.sN, self.sU = conv.sigma_xyz2enu(
                self.X,
                self.Y,
                self.Z,
                self.sX,
                self.sY,
                self.sZ,
                s_xy=sXY,
                s_yz=sYZ,
                s_xz=sXZ,
            )

    def UTMcalc_pt(self, ellips="wgs84"):
        """
        Calculate UTM projected coordinates.

        Parameters
        ----------
        ellips : str, optional
            Ellipsoid model to use. The default is 'wgs84'.

        Returns
        -------
        None

        Notes
        -----
        Updates the Eutm, Nutm, Uutm attributes.
        """
        self.Eutm, self.Nutm, _ = conv.utm_geo2xy(self.F, self.L)
        self.Uutm = self.H

    def keysanex(self):
        """
        Get the keys of the anex (annex data) dictionary.

        Returns
        -------
        list
            List of keys in the anex dictionary.
        """
        return list(self.anex.keys())

    def helmert_trans(self, params="itrf2008_2_etrf2000", invert=False):
        """
        Apply Helmert transformation to the Point coordinates.

        Parameters
        ----------
        params : str, optional
            Transformation parameter set name. The default is 'itrf2008_2_etrf2000'.
        invert : bool, optional
            If True, apply inverse transformation. The default is False.

        Returns
        -------
        None

        Notes
        -----
        Updates the Point coordinates in-place based on the Helmert transformation.
        """
        Xb = reffram.helmert_trans(np.array([self.X, self.Y, self.Z]), params, invert)
        self.XYZset(*Xb)
        return None

    def velocity_trans(self, vx, vy, vz, epoc_init="auto", epoc_end="auto"):
        """
        Apply velocity-based coordinate transformation over time.

        Parameters
        ----------
        vx : float
            Velocity in X direction (m/year).
        vy : float
            Velocity in Y direction (m/year).
        vz : float
            Velocity in Z direction (m/year).
        epoc_init : float or 'auto', optional
            Initial epoch in decimal years. If 'auto', uses the point's epoch.
            The default is 'auto'.
        epoc_end : float or 'auto', optional
            Final epoch in decimal years. If 'auto', uses the point's epoch.
            The default is 'auto'.

        Returns
        -------
        None

        Notes
        -----
        At least one of epoc_init or epoc_end must not be 'auto'.
        Updates the Point coordinates based on the velocity transformation.
        """

        if epoc_init == "auto" and epoc_end == "auto":
            log.error("epoc_init == auto and epoc_end == auto")
            return None

        tdt = conv.posix2dt(self.T)
        yeardec = conv.dt2year_decimal(tdt)

        if epoc_init == "auto":
            epoc_init = yeardec

        if epoc_end == "auto":
            epoc_end = yeardec

        Xb = reffram.itrf_speed_calc(
            self.X, self.Y, self.Z, epoc_init, vx, vy, vz, epoc_end
        )
        self.XYZset(*Xb)
        return None


class TimeSeriePoint:

    def __init__(self, stat="STAT"):
        """
        Initialize a TimeSeriePoint object

        Parameters
        ----------
        stat : str, optional
            station 4-char. code. The default is 'STAT'.

        """

        self.pts = []
        self.i_nomi = 0

        self.meta_set(stat=stat)

        self.boolENU = False
        self.boolUTM = False

        self.bool_interp_uptodate = False
        self.bool_discont = False
        self.bool_discont_manu = False

        self.anex = dict()

    def __repr__(self):
        """
        Return string representation of the TimeSeriePoint object.

        Returns
        -------
        str
            String containing station code, number of points, date range,
            number of days, and data coverage percentage.
        """
        if self.pts == []:
            raise Exception("ERR: TimeSeriePoint is empty ...")

        start = self.startdate()
        end = self.enddate()
        nbday = int((end - start).days + 1.0)

        ratio = self.nbpts * 100.0 / nbday

        return "{} {} {} {} {} {} {} {:5.2f}{}".format(
            self.stat, self.nbpts, "points", start, end, nbday, "nb days", ratio, "%"
        )

    def __getitem__(self, i):
        """
        Get a Point from the TimeSerie by index.

        Parameters
        ----------
        i : int
            Index of the Point.

        Returns
        -------
        Point
            Point object at index i.
        """
        return self.pts[i]

    @property
    def nbpts(self):
        """
        Method to have the length of the TimeSerie

        Returns
        -------
        int
            Length of the TimeSerie.

        """
        return len(self.pts)

    def meta_set(self, path="", stat="STAT", name=""):
        """
        Set metadata about the TimeSerie.

        Parameters
        ----------
        path : str, optional
            File path. The default is ''.
        stat : str, optional
            Station 4-character code. The default is 'STAT'.
        name : str, optional
            Free name for the TimeSerie (e.g., experiment, period, software name).
            The default is ''.

        Returns
        -------
        None
        """

        self.path = path
        self.stat = stat

        bn = os.path.basename(path)
        dn = os.path.dirname(path)

        if name == "":
            if bn == "tdp_final":
                self.name = os.path.basename(dn)
            else:
                self.name = bn
        else:
            self.name = name
        self.interval_nominal()

    def del_data(self):
        """
        Method to purge the data in the TimeSeriePoint

        Returns
        -------
        None.

        """
        self.pts = []

        self.bool_interp_uptodate = False
        self.bool_discont = False

    def readfile(self, filein):
        """
        Read time series data from a file.

        Parameters
        ----------
        filein : str
            Path of the file to read.

        Returns
        -------
        None

        Notes
        -----
        Should be used with care. Replaces all internal data.
        """
        self.__dict__ = files_rw.read_all_points(filein).__dict__

        self.interp_set()

    def add_point(self, point_inp):
        """
        Add a Point to the TimeSerie.

        Parameters
        ----------
        point_inp : Point
            Point object to add to the time series.

        Returns
        -------
        None

        Notes
        -----
        Marks interpolation as outdated and should be recalculated if needed.
        """
        self.pts.append(point_inp)
        # this line is discontiued, because now nbpts is a property
        # self.nbpts = len(self.pts)

        self.bool_interp_uptodate = False

    def aleapt(self):
        """
        Get a random Point from the TimeSeries.

        Returns
        -------
        Point
            Randomly selected Point object.
        """
        ipt = np.random.randint(self.nbpts)
        log.info("random point selected #%s", str(ipt))

        return self.pts[ipt]

    def startdate(self):
        """
        Get the first epoch of the data in the TimeSerie.

        Returns
        -------
        datetime.datetime
            First timestamp in the time series.
        """
        self.sort()
        return conv.posix2dt(self.pts[0].T)

    def enddate(self):
        """
        Get the last epoch of the data in the TimeSerie.

        Returns
        -------
        datetime.datetime
            Last timestamp in the time series.
        """
        self.sort()
        return conv.posix2dt(self.pts[-1].T)

    def len_period(self, output_seconds=False):
        """
        Get the period length of the TimeSerie.

        Parameters
        ----------
        output_seconds : bool, optional
            If True, return the result in seconds. If False, return as timedelta.
            The default is False.

        Returns
        -------
        int or datetime.timedelta
            Period length (as timedelta if output_seconds=False, as int seconds otherwise).
        """
        delta = self.enddate() - self.startdate()

        if output_seconds:
            return delta.days * 86400 + delta.seconds
        else:
            return delta

    def interval_nominal(self):
        """
        Get the nominal interval between consecutive epochs.

        Returns
        -------
        float
            Nominal interval between epochs in seconds (rounded to 1 decimal place).
        """

        if len(self.pts) < 2:
            self.i_nomi = 0
        else:
            Ttemp = np.sort([pt.T for pt in self.pts])
            self.i_nomi = np.round(np.min(np.diff(Ttemp)), 1)

        return self.i_nomi

    def from_list(self, T, A, B, C, coortype="XYZ", sA=[], sB=[], sC=[]):
        """
        Load data from lists into the TimeSerie.

        Parameters
        ----------
        T : list of float
            Times (POSIX timestamps).
        A : list of float
            First component (X, F latitude, or E depending on coortype).
        B : list of float
            Second component (Y, L longitude, or N depending on coortype).
        C : list of float
            Third component (Z, H height, or U depending on coortype).
        coortype : str, optional
            Coordinate type ('XYZ', 'FLH', 'ENU', 'NED', 'UTM'). The default is 'XYZ'.
        sA : list of float, optional
            Sigma (standard deviation) of A component. The default is [].
        sB : list of float, optional
            Sigma (standard deviation) of B component. The default is [].
        sC : list of float, optional
            Sigma (standard deviation) of C component. The default is [].

        Returns
        -------
        None
        """

        if not sA:
            sA = np.zeros(len(A))
        if not sB:
            sB = np.zeros(len(B))
        if not sC:
            sC = np.zeros(len(C))

        for t, a, b, c, sa, sb, sc in zip(T, A, B, C, sA, sB, sC):
            point = Point(a, b, c, t, coortype, sa, sb, sc)

            self.add_point(point)

        if coortype == "ENU":
            self.boolENU = True

        if coortype == "UTM":
            self.boolUTM = True

        self.sort()

        return None

    def to_list(self, coortype="XYZ", specific_output=None, time_as_datetime=False):
        """
        Export the TimeSerie as lists of numpy arrays.

        Parameters
        ----------
        coortype : str, optional
            Coordinate type to export ('XYZ', 'FLH', 'ENU', 'UTM').
            The default is 'XYZ'.
        specific_output : int, optional
            If specified, return only one array (index 0-6).
            The default is None (returns all 7 arrays).
        time_as_datetime : bool, optional
            If True, return time as datetime objects. If False, return POSIX timestamps.
            The default is False.

        Returns
        -------
        tuple of 7 numpy arrays
            (A, B, C, T, sA, sB, sC) where:
            - A, B, C depend on coortype (e.g., X, Y, Z for XYZ)
            - T is time (datetime or POSIX depending on time_as_datetime)
            - sA, sB, sC are standard deviations
        """

        if coortype == "XYZ":
            a_lbl, b_lbl, c_lbl = "X", "Y", "Z"
            sa_lbl, sb_lbl, sc_lbl = "sX", "sY", "sZ"

        elif coortype == "FLH":
            a_lbl, b_lbl, c_lbl = "F", "L", "H"
            sa_lbl, sb_lbl, sc_lbl = "sF", "sL", "sH"

        elif coortype == "ENU":
            if not self.boolENU:
                log.warning("no ENU coord. for " + self.name)
                return None

            a_lbl, b_lbl, c_lbl = "E", "N", "U"
            sa_lbl, sb_lbl, sc_lbl = "sE", "sN", "sU"

        elif coortype == "UTM":
            if not self.boolUTM:
                log.warning("no UTM coord. for " + self.name)
                return None

            a_lbl, b_lbl, c_lbl = "Eutm", "Nutm", "Uutm"
            sa_lbl, sb_lbl, sc_lbl = "sEutm", "sNutm", "sUutm"

        else:
            log.error("coortype does not exist")
            raise Exception

        if self.nbpts == 0:
            log.error(self.name + " the timeserie is empty")

        a = np.asarray([getattr(pt, a_lbl) for pt in self.pts])
        b = np.asarray([getattr(pt, b_lbl) for pt in self.pts])
        c = np.asarray([getattr(pt, c_lbl) for pt in self.pts])
        t = np.asarray([pt.T for pt in self.pts])

        if hasattr(self.pts[0], sa_lbl):
            sa = np.asarray([getattr(pt, sa_lbl) for pt in self.pts])
            sb = np.asarray([getattr(pt, sb_lbl) for pt in self.pts])
            sc = np.asarray([getattr(pt, sc_lbl) for pt in self.pts])
        else:
            sa = np.asarray([np.nan] * len(self.pts))
            sb = np.asarray([np.nan] * len(self.pts))
            sc = np.asarray([np.nan] * len(self.pts))
        # il faut squeezer les vecteurs parce que des fois on se retrouve
        # avec des matrices

        if time_as_datetime:
            tout = conv.posix2dt(t)
        else:
            tout = t

        sq = np.squeeze
        outtup = (sq(a), sq(b), sq(c), sq(tout), sq(sa), sq(sb), sq(sc))
        if specific_output is None:
            return outtup
        elif type(specific_output) is int:
            return outtup[specific_output]
        else:
            log.info("this mode must be implemented ;)")
            log.info("use an int as index instead")
            return outtup

    def to_dataframe(self, coortype="XYZ", anex_key_list=None):
        """
        Export the TimeSerie as a pandas DataFrame.

        Parameters
        ----------
        coortype : str or iterable of str, optional
            Coordinate type(s) to export ('XYZ', 'FLH', 'ENU', 'UTM').
            Can be a single string or tuple of strings for multiple coordinates.
            The default is 'XYZ'.
        anex_key_list : list of str, optional
            List of keys from the point's anex (annex) dictionary to add as columns.
            The default is None.

        Returns
        -------
        df : pandas.DataFrame
            DataFrame with columns for epoch, time, coordinates, and uncertainties.
        """

        if not utils.is_iterable(coortype):
            coortype = (coortype,)

        col_stk = tuple()
        col_name_stk = []

        for icoty, coty in enumerate(coortype):
            a, b, c, t, sa, sb, sc = self.to_list(coty)

            if coty == "UTM":
                cotycolnam = ["e_utm", "n_utm", "u_utm"]
            else:
                cotycolnam = coty

            if icoty == 0:
                tdt = conv.posix2dt(t)
                col_stk = col_stk + (tdt, t, a, b, c, sa, sb, sc)
                col_name_stk = (
                    ["epoch", "t"]
                    + [e.lower() for e in cotycolnam]
                    + ["s" + e.lower() for e in cotycolnam]
                )
            else:
                col_stk = col_stk + (a, b, c, sa, sb, sc)
                col_name_stk = [e.lower() for e in coty] + [
                    "s" + e.lower() for e in coty
                ]

            if anex_key_list:
                for key in anex_key_list:
                    val_list = [
                        pt.anex[key] if key in pt.anex.keys() else np.nan
                        for pt in self.pts
                    ]
                    col_stk = col_stk + (val_list,)
                    col_name_stk = col_name_stk + [key]

        big = np.column_stack(col_stk)
        df = pd.DataFrame(big)
        df.columns = col_name_stk
        df = df.infer_objects()

        return df


    def sort(self):
        """
        Sort points in the TimeSerie by time (in-place).

        Returns
        -------
        None

        Notes
        -----
        Modifies the internal list of points in chronological order.
        """
        self.pts.sort(key=lambda x: x.T)

    def plot(
        self,
        coortype="ENU",
        diapt=1.5,
        alpha=0.8,
        fig=1,
        errbar=True,
        symbol=".",
        errbar_width=1,
        ylim=None,
        legend_loc="best",
        legend_ncol=1,
    ):
        """
        Plot data in a TimeSerie Object

        Parameters
        ----------
        coortype : str, optional
            The coordinates type. The default is 'ENU'.
        diapt : float, optional
            Point diameter. The default is 1.
        alpha : float, optional
            Alpha (transparency) of points.
            The default is 0.8.
        fig : int or Figure object, optional
            Figure ID where the data will be plotted
            can accept a int (id of a Figure)
            OR the figure Object itself.
            The default is 1.
        errbar : bool, optional
            Plot the error bars. The default is True.
        symbol : str, optional
            symbol. The default is '.'.
        errbar_width : float, optional
            coefficient for the error bar size. The default is 1.
        ylim : tuple, optional
            Y-axis limits. The default is None.
        legend_loc : str, optional
            Location of the legend. The default is 'best'.
        legend_ncol : int, optional
            Number of columns in the legend. The default is 1.

        Returns
        -------
        figobj : Figure object
            figure object of the plot.
        axes : array of Axes objects
            array of the 3 axes objects of the plot.

        """

        try:
            A, B, C, T, sA, sB, sC = self.to_list(coortype=coortype)
        except TypeError as tyer:
            log.error("unable to get coordinates")
            log.info("TRICK : check if the given coortype is in the timeserie")
            raise tyer

        coor_tup = (A, B, C, T, sA, sB, sC)

        return time_series.plot_timeseries(
            coor_tup,
            coortype=coortype,
            diapt=diapt,
            alpha=alpha,
            fig=fig,
            errbar=errbar,
            symbol=symbol,
            errbar_width=errbar_width,
            ylim=ylim,
            legend_loc=legend_loc,
            legend_ncol=legend_ncol,
            name=self.name,
            stat=self.stat,
        )

    def plot_discont(self, fig=1):
        """
        Plot discontinuities of the TimeSerie on existing plot.

        Parameters
        ----------
        fig : int or matplotlib.figure.Figure, optional
            Figure ID or Figure object where discontinuities will be plotted.
            The default is 1.

        Returns
        -------
        None
        """

        if not self.bool_discont:
            log.warning("no discontinuities in TimeSerie")
            return None

        if type(fig) is int:
            figobj = plt.figure(fig)
        elif type(fig) is plt.Figure:
            figobj = fig

        for ax in figobj.axes:
            stats.plot_vertical_bar_ax(self.discont, ax, "r")

        if self.bool_discont_manu:
            for ax in figobj.axes:
                stats.plot_vertical_bar_ax(self.discont_manu, ax, "g")
        return None

    def discont_manu_click(self, fig=1):
        """
        Interactively record manual discontinuities by clicking on plot.

        Parameters
        ----------
        fig : int or matplotlib.figure.Figure, optional
            Figure ID or Figure object for interaction. The default is 1.

        Returns
        -------
        tuple
            (multi, cid) cursor objects that must be stored as global variables.

        Notes
        -----
        Manual discontinuities are recorded in both the "main" discont list
        and a separate discont_manu list for identification.

        Use SPACE key to record a discontinuity at the cursor position.

        Important: cursor objects must be stored as global variables:

            multi, cid = tsout.discont_manu_click()

        See Also
        --------
        point_n_click_plot : More complete interactive plotting method.
        """

        if type(fig) is int:
            figobj = plt.figure(fig)
        elif type(fig) is plt.Figure:
            figobj = fig

        if not self.bool_discont:
            self.discont = []

        if not self.bool_discont_manu:
            self.discont_manu = []

        log.info("press SPACE to record a manual discontinuity")

        def onclick_discont(event):
            ix, iy = (
                matplotlib.dates.num2date(event.xdata).replace(tzinfo=None),
                event.ydata,
            )

            log.info("discontinuity recorded : %s", ix)
            for ax in figobj.axes:
                stats.plot_vertical_bar_ax([ix], ax, "g")
            #

            self.bool_discont = True
            self.bool_discont_manu = True

            self.discont.append(ix)
            self.discont = sorted(self.discont)

            self.discont_manu.append(ix)
            self.discont_manu = sorted(self.discont_manu)

            # figobj.show()
            plt.draw()

            return None

        multi = MultiCursor(figobj.canvas, figobj.axes, color="k", lw=1)
        cid = figobj.canvas.mpl_connect("key_press_event", onclick_discont)

        return multi, cid

    def initype(self):
        """
        Get the most common coordinate type in the TimeSerie.

        Returns
        -------
        str
            The coordinate type that appears most frequently ('XYZ', 'FLH', 'ENU', etc.).
        """

    def ENUcalc(self, refENU):
        """
        Calculate ENU coordinates relative to a reference point or time series.

        Parameters
        ----------
        refENU : Point or TimeSeriePoint
            Reference point or time series for ENU calculation.

        Returns
        -------
        None

        Notes
        -----
        If refENU is a TimeSeriePoint, interpolation is performed to get
        reference coordinates at each measurement epoch.
        """
        if refENU.__class__.__name__ == "Point":
            self.refENU = refENU
            [pt.ENUcalc_pt(refENU) for pt in self.pts]
            self.boolENU = True
            self.bool_interp_uptodate = False
            self.interp_set()
        elif refENU.__class__.__name__ == "TimeSeriePoint":
            self.refENU = refENU.mean_posi()
            refENU.interp_set()
            [
                pt.ENUcalc_pt(
                    Point(
                        A=refENU.XfT(pt.T),
                        B=refENU.YfT(pt.T),
                        C=refENU.ZfT(pt.T),
                        initype="XYZ",
                        T=pt.T,
                    )
                )
                for pt in self.pts
            ]

            self.boolENU = True
            self.bool_interp_uptodate = False
            self.interp_set()

    def ENUcalc_from_mean_posi(self, mean_type="median"):
        """
        Calculate ENU coordinates relative to the mean/median position of the TimeSerie.

        Parameters
        ----------
        mean_type : str, optional
            Type of mean to use ('median' or 'mean'). The default is 'median'.

        Returns
        -------
        None
        """
        self.ENUcalc(self.mean_posi(mean_type="median"))
        return None

    def ENUcalc_from_first_posi(self):
        """
        Calculate ENU coordinates relative to the first position in the TimeSerie.

        Returns
        -------
        None
        """
        self.ENUcalc(self.pts[0])
        return None

    def from_uniq_point(self, Point, startdate, enddate, pas=1):
        """
        Create a TimeSerie from a single Point with defined time range.

        Parameters
        ----------
        Point : Point
            Reference Point object to replicate at different epochs.
        startdate : float or datetime.datetime
            Start date (POSIX timestamp or datetime object).
        enddate : float or datetime.datetime
            End date (POSIX timestamp or datetime object).
        pas : float, optional
            Time step in seconds. The default is 1.

        Returns
        -------
        None

        Notes
        -----
        Creates a time series by replicating the Point at regular intervals
        from startdate to enddate.
        """

        if type(startdate) == dt.datetime:
            startdate = conv.dt2posix(startdate)
        if type(enddate) == dt.datetime:
            enddate = conv.dt2posix(enddate)

        N = int(np.round(enddate - startdate / pas))
        # datelist = np.arange(startdate,enddate,pas)

        Point.Tset(startdate)

        for i in range(N):
            Point.Tset(Point.T + pas)
            self.add_point(copy.copy(Point))

        self.i_nomi = pas

    def UTMcalc(self):
        """
        Calculate UTM projected coordinates for all points in the TimeSerie.

        Returns
        -------
        None

        Notes
        -----
        Updates the UTM coordinates (Eutm, Nutm, Uutm) for all points.
        Sets boolUTM to True.
        """
        self.boolUTM = True
        [pt.UTMcalc_pt() for pt in self.pts]

    def time_win(self, windows, mode="keep"):
        """
        Apply a time window filter to the TimeSerie.

        Parameters
        ----------
        windows : list of tuple or list of datetime
            Time windows to process. Format depends on implementation.
        mode : str, optional
            Filter mode ('keep' or 'remove'). The default is 'keep'.

        Returns
        -------
        None

        Notes
        -----
        This operation replaces the internal data. Be cautious when applying
        this method as it modifies the TimeSerie in place.
        """
        self.__dict__ = time_series.time_win(self, windows, mode).__dict__

    def interp_set(self, interptype="slinear"):
        """
        Set up coordinate interpolators for the TimeSerie.

        Parameters
        ----------
        interptype : str, optional
            Interpolation type. The default is 'slinear' (linear).

        Returns
        -------
        None

        Notes
        -----
        Creates interpolation functions for each available coordinate type
        (ENU, XYZ, FLH, UTM). Interpolators are stored as attributes
        (EfT, NfT, UfT, etc.) for later use.
        """
        if (not hasattr(self.pts[0], "E")) or np.isnan(self.pts[0].E) == True:
            log.warning("no ENU for " + self.name)
        else:
            E, N, U, T, _, _, _ = self.to_list("ENU")
            self.EfT = scipy.interpolate.interp1d(
                T, E, bounds_error=False, kind=interptype
            )
            self.NfT = scipy.interpolate.interp1d(
                T, N, bounds_error=False, kind=interptype
            )
            self.UfT = scipy.interpolate.interp1d(
                T, U, bounds_error=False, kind=interptype
            )

        if (not hasattr(self.pts[0], "X")) or np.isnan(self.pts[0].X) == True:
            log.warning("no XYZ for " + self.name)
        else:

            X, Y, Z, T, _, _, _ = self.to_list("XYZ")

            self.XfT = scipy.interpolate.interp1d(
                T, X, bounds_error=False, kind=interptype
            )
            self.YfT = scipy.interpolate.interp1d(
                T, Y, bounds_error=False, kind=interptype
            )
            self.ZfT = scipy.interpolate.interp1d(
                T, Z, bounds_error=False, kind=interptype
            )

        if (not hasattr(self.pts[0], "F")) or np.isnan(self.pts[0].L) == True:
            log.warning("no FLH for " + self.name)
        else:
            F, L, H, T, _, _, _ = self.to_list("FLH")

            self.FfT = scipy.interpolate.interp1d(
                T, F, bounds_error=False, kind=interptype
            )
            self.LfT = scipy.interpolate.interp1d(
                T, L, bounds_error=False, kind=interptype
            )
            self.HfT = scipy.interpolate.interp1d(
                T, H, bounds_error=False, kind=interptype
            )

        if (not hasattr(self.pts[0], "Eutm")) or np.isnan(self.pts[0].Eutm) == True:
            # log.warning("no UTM for " + self.name)
            pass
        else:
            Eutm, Nutm, Uutm, T, _, _, _ = self.to_list("UTM")

            self.EutmfT = scipy.interpolate.interp1d(
                T, Eutm, bounds_error=False, kind=interptype
            )
            self.NutmfT = scipy.interpolate.interp1d(
                T, Nutm, bounds_error=False, kind=interptype
            )
            self.UutmfT = scipy.interpolate.interp1d(
                T, Uutm, bounds_error=False, kind=interptype
            )

        self.bool_interp_uptodate = True

    def interp_get(self, T, coortype="ENU"):
        """
        Get interpolated coordinates at given times.

        Parameters
        ----------
        T : float or list of float
            Time(s) in POSIX format where interpolation is desired.
        coortype : str, optional
            Coordinate type ('ENU', 'XYZ', 'FLH', 'UTM'). The default is 'ENU'.

        Returns
        -------
        TimeSeriePoint
            New TimeSerie with interpolated points at requested times.

        Notes
        -----
        Interpolators must be set up first using interp_set().
        If interpolators are outdated, they are automatically recalculated.
        """

        if self.bool_interp_uptodate == False:
            log.warning("interp obsolete, recalcul auto")
            self.interp_set()

        tsout = copy.copy(self)
        tsout.del_data()

        if not utils.is_iterable(T):
            T = np.array([T])

        if coortype == "ENU":
            A = self.EfT(T)
            B = self.NfT(T)
            C = self.UfT(T)

        if coortype == "XYZ":
            A = self.XfT(T)
            B = self.YfT(T)
            C = self.ZfT(T)

        if coortype == "FLH":
            A = self.FfT(T)
            B = self.LfT(T)
            C = self.HfT(T)

        if coortype == "UTM":
            A = self.EutmfT(T)
            B = self.NutmfT(T)
            C = self.UutmfT(T)

        else:
            log.error("coortype does not exist")
            raise Exception

        for i in range(len(T)):
            tsout.add_point(Point(A[i], B[i], C[i], T=T[i], initype=coortype))

        return tsout

    def set_discont(self, indiscont):
        """
        Set the list of discontinuities in the TimeSerie.

        Parameters
        ----------
        indiscont : list of datetime or float
            List of discontinuity times (datetime objects or POSIX timestamps).

        Returns
        -------
        None

        Notes
        -----
        Discontinuities are typically used for visualization and
        can be detected or manually set.
        """
        self.discont = indiscont
        self.bool_discont = True

    def mean_posi(self, coortype="XYZ", outtype="point", mean_type="median"):
        """
        Calculate the mean position of the TimeSerie.

        Parameters
        ----------
        coortype : str, optional
            Coordinate type for calculation ('XYZ', 'FLH', 'ENU', 'UTM').
            The default is 'XYZ'.
        outtype : str, optional
            Output format ('point' for Point object, 'tuple' for coordinates).
            The default is 'point'.
        mean_type : str, optional
            Type of mean ('mean' or 'median'). The default is 'median'.

        Returns
        -------
        Point or tuple
            Mean position as Point object or tuple of coordinates (A, B, C, T).
        """

        # special case where only one point
        if self.nbpts == 1:
            return copy.copy(self)

        A, B, C, T, sA, sB, sC = self.to_list(coortype=coortype)

        if mean_type == "mean":
            Aout = np.nanmean(A)
            Bout = np.nanmean(B)
            Cout = np.nanmean(C)
        elif mean_type == "median":
            Aout = np.nanmedian(A)
            Bout = np.nanmedian(B)
            Cout = np.nanmedian(C)
        else:
            log.error("mean_type does not exist")
            raise Exception

        Tout = (np.max(T) - np.min(T)) / 2

        if outtype == "point":
            out = Point(Aout, Bout, Cout, Tout)
        elif outtype == "tuple":
            out = (Aout, Bout, Cout, Tout)
        else:
            out = (Aout, Bout, Cout, Tout)

        return out

    def add_offset(self, dA, dB, dC, coortype="ENU"):
        """
        Add an offset to all points in the TimeSerie.

        Parameters
        ----------
        dA : float
            Offset for the A component.
        dB : float
            Offset for the B component.
        dC : float
            Offset for the C component.
        coortype : str, optional
            Coordinate type for offset ('ENU', 'XYZ', etc.). The default is 'ENU'.

        Returns
        -------
        None
        """
        for pt in self.pts:
            pt.add_offset(dA, dB, dC, coortype="ENU")

    def decimate(self, dec):
        """
        Decimate the TimeSerie by keeping 1 out of every dec points.

        Parameters
        ----------
        dec : int
            Decimation factor (keep 1/dec points).

        Returns
        -------
        None

        Notes
        -----
        Modifies the TimeSerie in place, removing points to reduce data density.
        """
        self.__dict__ = time_series.decimate_cleaner(self, dec).__dict__
        # decimate_cleaner(self,dec,True)

    def find_point(self, tin, tol=0.001, stop_when_found=True):
        """
        Find a Point by timestamp with tolerance.

        Parameters
        ----------
        tin : float or datetime.datetime
            Target timestamp (POSIX or datetime).
        tol : float, optional
            Tolerance in seconds. The default is 0.001.
        stop_when_found : bool, optional
            If True, stop at first match. If False, return all matches.
            The default is True.

        Returns
        -------
        tuple
            If stop_when_found=True: (Point, int) - Point and its index
            If stop_when_found=False: (list of Point, list of int) - Points and indices
        """

        find = False

        if type(tin) is dt.datetime:
            tin = conv.dt2posix(tin)

        tmin = self.pts[0].T
        tmax = self.pts[-1].T

        if not (tmin <= tin <= tmax):
            log.warning("t !€ [tmin , tmax]")

        if stop_when_found:
            for i, p in enumerate(self.pts):
                if np.abs(p.T - tin) < tol:
                    find = True
                    break
            if find:
                return p, i
            else:
                return Point(), np.nan

        else:
            pts_stk = []
            i_stk = []
            for i, p in enumerate(self.pts):
                if np.abs(p.T - tin) < tol:
                    pts_stk.append(p)
                    i_stk.append(i)
            return pts_stk, i_stk

    def rm_duplicat_pts(self, coortype="XYZ"):
        """
        Remove duplicate points from the time series.

        Parameters
        ----------
        coortype : str, optional
            The coordinate type to consider for duplication check. The default is "XYZ".

        Returns
        -------
        None
        """
        t = self.to_dataframe(coortype)["t"]

        dup_bool = t.duplicated()

        if dup_bool.sum() > 0:
            log.warning(
                "%s duplicated point(s) removed for %s", dup_bool.sum(), self.name
            )

        self.pts = list(pd.Series(self.pts)[np.logical_not(dup_bool)])
