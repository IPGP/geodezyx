# -*- coding: utf-8 -*-
"""
Created on Fri Apr 3 2026

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


#  ______                      _                      _        _    _____ _
# |  ____|                    (_)                    | |      | |  / ____| |
# | |__  __  ___ __   ___ _ __ _ _ __ ___   ___ _ __ | |_ __ _| | | |    | | __ _ ___ ___  ___  ___
# |  __| \ \/ / '_ \ / _ \ '__| | '_ ` _ \ / _ \ '_ \| __/ _` | | | |    | |/ _` / __/ __|/ _ \/ __|
# | |____ >  <| |_) |  __/ |  | | | | | | |  __/ | | | || (_| | | | |____| | (_| \__ \__ \  __/\__ \
# |______/_/\_\ .__/ \___|_|  |_|_| |_| |_|\___|_| |_|\__\__,_|_|  \_____|_|\__,_|___/___/\___||___/
#             | |
#             |_|


class Attitude:
    """Attitude (Roll, Pitch, Yaw) representation."""

    def __init__(
        self, R=0, P=0, Y=0, T=0, sR=0, sP=0, sY=0, devID="NULL", angtype="deg"
    ):
        """
        Initialize an Attitude object.

        Parameters
        ----------
        R : float, optional
            Roll angle. The default is 0.
        P : float, optional
            Pitch angle. The default is 0.
        Y : float, optional
            Yaw angle. The default is 0.
        T : float or datetime.datetime, optional
            Epoch time. The default is 0 (current time).
        sR : float, optional
            Standard deviation of Roll. The default is 0.
        sP : float, optional
            Standard deviation of Pitch. The default is 0.
        sY : float, optional
            Standard deviation of Yaw. The default is 0.
        devID : str, optional
            Device identifier. The default is 'NULL'.
        angtype : str, optional
            Angle type ('deg' for degrees or 'rad' for radians). The default is 'deg'.

        Raises
        ------
        Exception
            If angtype is not 'deg' or 'rad'.
        """

        self.Tset(T)
        self.devID = devID

        if angtype == "deg":
            self.R = R
            self.P = P
            self.Y = Y

        elif angtype == "rad":
            self.R = np.rad2deg(R)
            self.P = np.rad2deg(P)
            self.Y = np.rad2deg(Y)
        else:
            raise Exception("Mauvais angtype")

        self.Qcalc()

    def __call__(self):
        """
        Return the Attitude as a tuple of Roll, Pitch, Yaw, and datetime.

        Returns
        -------
        tuple
            (Roll, Pitch, Yaw, datetime)
        """
        return self.R, self.P, self.Y, self.Tdt

    def __repr__(self):
        """
        Return string representation of the Attitude.

        Returns
        -------
        str
            String in format "R,P,Y,datetime"
        """
        return "{},{},{},{}".format(self.R, self.P, self.Y, self.Tdt)

    def Tset(self, T=0):
        """
        Set the time/epoch of the Attitude.

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

    def RPYget(self):
        """
        Get Roll, Pitch, Yaw angles.

        Returns
        -------
        tuple
            (Roll, Pitch, Yaw) in degrees.
        """
        return self.R, self.P, self.Y

    def RPYset(self, R=0, P=0, Y=0, sR=0, sP=0, sY=0):
        """
        Set Roll, Pitch, Yaw angles and their standard deviations.

        Parameters
        ----------
        R : float, optional
            Roll angle in degrees. The default is 0.
        P : float, optional
            Pitch angle in degrees. The default is 0.
        Y : float, optional
            Yaw angle in degrees. The default is 0.
        sR : float, optional
            Standard deviation of Roll. The default is 0.
        sP : float, optional
            Standard deviation of Pitch. The default is 0.
        sY : float, optional
            Standard deviation of Yaw. The default is 0.

        Returns
        -------
        None
        """

    def Qcalc(self):
        """
        Calculate quaternion from Roll, Pitch, Yaw angles.

        Returns
        -------
        None

        Notes
        -----
        Stores the quaternion in self.Q.
        """
        self.Q = conv.quaternion(self.R, self.P, self.Y, "deg")
        return None


class TimeSerieObs(object):
    """
    Time series observation class for homogeneous observation data.

    Unlike TimeSeriePoint, objects contain only one type of observation data
    in a single coordinate form. When reading from files with multiple devices,
    the read functions return a list of TimeSerieObs objects.

    Attributes
    ----------
    obs : list
        List of observation objects.
    typeobs : str
        Type of observations (e.g., 'RPY' for Roll-Pitch-Yaw).
    nbobs : int
        Number of observations.
    """

    def __init__(self, typeobs="NULL", filepath=""):
        """
        Initialize a TimeSerieObs object.

        Parameters
        ----------
        typeobs : str, optional
            Type of observations ('RPY' or other). The default is 'NULL'.
        filepath : str, optional
            File path for metadata. The default is ''.

        Returns
        -------
        None
        """

        self.obs = []
        self.nbobs = 0
        self.i_nomi = 0
        self.typeobs = typeobs
        # ASM Attitude
        self.bool_interp_uptodate = False

        self.meta_set(filepath)

    def meta_set(self, path="", devID="NULL", name=""):
        """
        Set metadata for the TimeSerieObs.

        Parameters
        ----------
        path : str, optional
            File path. The default is ''.
        devID : str, optional
            Device identifier. The default is 'NULL'.
        name : str, optional
            Observation series name. The default is ''.

        Returns
        -------
        None
        """
        self.path = path
        self.devID = devID

        bn = os.path.basename(path)
        dn = os.path.dirname(path)

        if name == "":
            if bn == "tdp_final":
                self.name = os.path.basename(dn)
            else:
                self.name = bn

        self.interval_nominal()

    def del_data(self):
        """
        Delete all observations from the TimeSerieObs.

        Returns
        -------
        None
        """

    def readfile(self, filein, indtab=0):
        """
        Read observation data from a file.

        Parameters
        ----------
        filein : str
            Path to the file to read.
        indtab : int, optional
            Index of the device/table to read. The default is 0.

        Returns
        -------
        None

        Notes
        -----
        Replaces all internal data with data from file.
        """

    def add_obs(self, inObs):
        """
        Add an observation to the TimeSerieObs.

        Parameters
        ----------
        inObs : Attitude or observation object
            Observation to add.

        Returns
        -------
        None

        Notes
        -----
        Marks interpolation as outdated.
        """

    def aleaobs(self):
        """
        Get a random observation from the TimeSerieObs.

        Returns
        -------
        observation object
            Randomly selected observation.
        """

    def interval_nominal(self):
        """
        Get the nominal interval between consecutive observations.

        Returns
        -------
        float
            Nominal interval in seconds (rounded to 1 decimal place).
        """

    def timewin(self, windows, mode="keep"):
        """
        Apply a time window filter to the TimeSerieObs.

        Parameters
        ----------
        windows : list
            Time windows to process.
        mode : str, optional
            Filter mode ('keep' or 'remove'). The default is 'keep'.

        Returns
        -------
        None
        """
        self.__dict__ = time_series.time_win(self, windows, mode).__dict__

    def startdate(self):
        """
        Get the first epoch of observations.

        Returns
        -------
        float
            First POSIX timestamp.
        """
        return self.obs[0].T

    def enddate(self):
        """
        Get the last epoch of observations.

        Returns
        -------
        float
            Last POSIX timestamp.
        """
        return self.obs[-1].T

    def to_list(self):
        """
        Export observations as lists of numpy arrays.

        Returns
        -------
        tuple of 7 numpy arrays
            (A, B, C, T, sA, sB, sC) where:
            - A, B, C are the main observation components (e.g., R, P, Y for RPY)
            - T is time (POSIX timestamps)
            - sA, sB, sC are standard deviations
        """

        if self.typeobs == "RPY":
            A, B, C = "R", "p", "Y"
            sA, sB, sC = "sR", "sP", "sY"
        else:
            log.error("typeobs not recognized")
            raise Exception

        A = np.asarray([getattr(o, A) for o in self.obs])
        B = np.asarray([getattr(o, B) for o in self.obs])
        C = np.asarray([getattr(o, C) for o in self.obs])
        T = np.asarray([o.T for o in self.obs])

        if hasattr(self, sA):
            sA = np.asarray([getattr(o, sA) for o in self.obs])
            sB = np.asarray([getattr(o, sB) for o in self.obs])
            sC = np.asarray([getattr(o, sC) for o in self.obs])

        else:
            sA = np.asarray([np.nan] * len(self.obs))
            sB = np.asarray([np.nan] * len(self.obs))
            sC = np.asarray([np.nan] * len(self.obs))

        return A, B, C, T, sA, sB, sC

    def interp_set(self, interptype="slinear"):
        """
        Set up observation interpolators.

        Parameters
        ----------
        interptype : str, optional
            Interpolation type (e.g., 'slinear'). The default is 'slinear'.

        Returns
        -------
        None

        Notes
        -----
        Creates interpolation functions for each observation component.
        For RPY observations, creates RfT, PfT, YfT interpolators.
        """

        if self.typeobs == "RPY":

            R, P, Y, T, _, _, _ = self.to_list()

            self.RfT = scipy.interpolate.interp1d(
                T, R, bounds_error=False, kind=interptype
            )
            self.PfT = scipy.interpolate.interp1d(
                T, P, bounds_error=False, kind=interptype
            )
            self.YfT = scipy.interpolate.interp1d_ang(
                T, Y, bounds_error=False, kind=interptype
            )

        self.bool_interp_uptodate = True

    def interp_get(self, T):
        """
        Get interpolated observations at given times.

        Parameters
        ----------
        T : float or list of float
            Time(s) in POSIX format where interpolation is desired.

        Returns
        -------
        TimeSerieObs
            New TimeSerieObs with interpolated observations at requested times.

        Notes
        -----
        Interpolators must be set up first using interp_set().
        """

        if self.bool_interp_uptodate == False:
            log.warning("interp obsolete, recalcul auto")
            self.interp_set()

        tsout = copy.copy(self)
        tsout.del_data()

        if not utils.is_iterable(T):
            T = np.array([T])

        if self.typeobs == "NULL":
            log.error("no typeobs defined (NULL)")
            return 0

        if self.typeobs == "RPY":
            A = self.RfT(T)
            B = self.PfT(T)
            C = self.YfT(T)

            for i in range(len(T)):
                tsout.add_obs(Attitude(A[i], B[i], C[i], T=T[i]))

        return tsout

    def plot(self, diapt=10, alpha=0.8, fig=1, new_style=True):
        """
        Plot observations as multi-panel figure.

        Parameters
        ----------
        diapt : float, optional
            Point marker size. The default is 10.
        alpha : float, optional
            Alpha (transparency) of points. The default is 0.8.
        fig : int or matplotlib.figure.Figure, optional
            Figure ID or Figure object. The default is 1.
        new_style : bool, optional
            If True, use 410 layout. If False, use 220 layout. The default is True.

        Returns
        -------
        None

        Notes
        -----
        Creates a 4-panel plot with phase space and time series for each component.
        """

        if self.typeobs == "RPY":
            listtitle = ["", "Roll", "Pitch", "Yaw"]
        else:
            listtitle = ["", "", "", ""]

        namest = 0
        namend = 10
        Tdt = conv.posix2dt(T)

        if new_style:
            styleint = 410
        else:
            styleint = 220

        plt.figure(fig)
        plt.subplot(styleint + 1)
        plt.axis("equal")
        plt.plot(
            A,
            B,
            ".",
            label=str(self.devID)[namest:namend],
            markersize=diapt,
            alpha=alpha,
        )
        plt.legend()
        plt.title(listtitle[0])

        plt.subplot(styleint + 2)
        plt.plot(
            Tdt,
            A,
            ".",
            label=str(self.devID)[namest:namend],
            markersize=diapt,
            alpha=alpha,
        )
        plt.legend()
        plt.title(listtitle[1])

        plt.subplot(styleint + 3)
        plt.plot(
            Tdt,
            B,
            ".",
            label=str(self.devID)[namest:namend],
            markersize=diapt,
            alpha=alpha,
        )
        plt.legend()
        plt.title(listtitle[2])

        plt.subplot(styleint + 4)
        plt.plot(
            Tdt,
            C,
            ".",
            label=str(self.devID)[namest:namend],
            markersize=diapt,
            alpha=alpha,
        )
        plt.legend()
        plt.title(listtitle[3])

    #  _____       _                      _   _             _____  _       _
    # |_   _|     | |                    | | (_)           |  __ \| |     | |
    #   | |  _ __ | |_ ___ _ __ __ _  ___| |_ ___   _____  | |__) | | ___ | |_
    #   | | | '_ \| __/ _ \ '__/ _` |/ __| __| \ \ / / _ \ |  ___/| |/ _ \| __|
    #  _| |_| | | | ||  __/ | | (_| | (__| |_| |\ v /  __/ | |    | | (_) | |_
    # |_____|_| |_|\__\___|_|  \__,_|\___|\__|_|_\_/ \___| |_|    |_|\___/ \__|
    # |  __ \    (_)     | |     ___     / ____| (_)    | |
    # | |__) |__  _ _ __ | |_   ( _ )   | |    | |_  ___| | __
    # |  ___/ _ \| | '_ \| __|  / _ \/\ | |    | | |/ __| |/ /
    # | |  | (_) | | | | | |_  | (_>  < | |____| | | (__|   <
    # |_|   \___/|_|_| |_|\__|  \___/\/  \_____|_|_|\___|_|\_\

    """
    USAGE :

        Then :

            PnC = point_n_click_plot()

            multi , cid = PnC(fig=6)

            PnC.selectedX

        i.e. :

            Create an object point_n_click_plot (here it is PnC in the exemple below)

            Call the object like a function with the id of the plot figure

            Make your selection using SPACE key

            Get your results in a list called PnC.selectedX




    """


class point_n_click_plot:
    """
    Interactive point selection tool for plots.

    Allows interactive selection of points on a plot by clicking, useful for
    identifying offsets, discontinuities, or other features of interest.

    Attributes
    ----------
    selectedX : list
        List of selected X-values (timestamps if Xdata_are_time=True).
    ver_bar_stk : list
        Stack of vertical bar objects for visualization.

    Examples
    --------
    Basic usage:

    >>> PnC = point_n_click_plot()
    >>> multi, cid = PnC(fig=1, Xdata_are_time=True)
    >>> selected_times = PnC.selectedX

    Notes
    -----
    - Press SPACE to record a point
    - Press R to remove the last recorded point
    - Cursor objects (multi, cid) must be stored as global variables
      to keep them in scope during interactive use

    See Also
    --------
    TimeSeriePoint.discont_manu_click : Similar functionality for discontinuities
    """

    def __init__(self):
        """
        Initialize the point_n_click_plot object.

        Returns
        -------
        None
        """
        self.selectedX = []
        self.ver_bar_stk = []

    def __call__(self, fig=1, Xdata_are_time=True):
        """
        Activate interactive point selection on a figure.

        Parameters
        ----------
        fig : int or matplotlib.figure.Figure, optional
            Figure ID or Figure object for interaction. The default is 1.
        Xdata_are_time : bool, optional
            If True, X-data are treated as timestamps (datetime).
            If False, X-data are treated as regular float values.
            The default is True.

        Returns
        -------
        tuple
            (multi, cid) cursor objects that must be stored as global variables.

        Notes
        -----
        - Press SPACE to record the current X-value
        - Press R to remove the previously recorded value

        Important: cursor objects must be stored as global variables:

            multi, cid = PnC(fig=1)
        """

        if type(fig) is int:
            figobj = plt.figure(fig)
        elif type(fig) is matplotlib.figure.Figure:
            figobj = fig

        log.info(
            "press SPACE to record a X-value, \n       press R to Remove the previously recorded one"
        )

        def onclick_discont(event):

            if event.key == " ":
                if Xdata_are_time:
                    ix, iy = (
                        matplotlib.dates.num2date(event.xdata).replace(tzinfo=None),
                        event.ydata,
                    )
                else:
                    ix, iy = event.xdata, event.ydata

                log.info("X value recorded : ", ix)

                for ax in figobj.axes:
                    out_bar_list = stats.plot_vertical_bar_ax(
                        [ix], ax, "b", linewidth=1
                    )

                    self.ver_bar_stk.append(out_bar_list[0])

                self.selectedX.append(ix)

                plt.draw()

            elif event.key in ("r", "R") and len(self.selectedX) > 0:
                last = self.selectedX[-1]

                self.selectedX.remove(last)

                last_bars = self.ver_bar_stk[-len(figobj.axes) :]

                for bar in last_bars:
                    self.ver_bar_stk.remove(bar)
                    bar.remove()

                log.info("value removed : ", last)

                plt.draw()

            return None

        multi = MultiCursor(figobj.canvas, figobj.axes, color="k", lw=1)
        cid = figobj.canvas.mpl_connect("key_press_event", onclick_discont)

        return multi, cid
