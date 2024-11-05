# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 13:55:33 2019

@author: psakicki
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

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

##########  END IMPORT  ##########

class Point():
      
    def __init__(self,A=0.,B=0.,C=0.,T=0.,initype='XYZ',
                 sA=0.,sB=0.,sC=0.,name='noname',anex=None):
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
        sA : TYPE, optional
            sigma of A component. The default is 0.
        sB : TYPE, optional
            sigma of B component. The default is 0.
        sC : TYPE, optional
            sigma of C component. The default is 0.
        name : str, optional
            Flexible name for the Point identification. The default is 'noname'.
        anex : dict, optional
            Additional data. The default is None. See Note
        
        Note
        ----
        
        A dictionary called anex is also initialized to allow a 
        versatile storage of a variety of data
        
        Exemple of dictionary keys 
        RMS: average RMS (for gipsy)
        sdAB , sdBC , sdAC : the variances between A,B,C (for rtklib)
        sdXY , sdXZ , sdYZ : the variances between XYZ (pbo.pos)
        Vx , Vy , Vz , sVx , sVy , sVz : velocity of the point (EPOS coordinates)
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

        if initype == 'XYZ':
            self.XYZset(A,B,C,sA,sB,sC)
        elif initype == 'FLH': # On travaille en degres decimaux
            self.FLHset(A,B,C,sA,sB,sC)

        elif initype == 'ENU':
            self.ENUset(A,B,C,sA,sB,sC)

        elif initype == 'NED':
            self.NEDset(A,B,C,sA,sB,sC)
            
        elif initype == 'UTM':
            self.UTMset(A,B,C,sA,sB,sC)
        else:
            log.error("wrong initype")

    def __call__(self):

        if self.initype == 'XYZ':
            return self.X,self.Y,self.Z,self.Tdt,self.T
        elif self.initype == 'FLH':
            return self.F,self.L,self.H,self.Tdt,self.T
        elif self.initype == 'ENU':
            return self.E,self.N,self.U,self.Tdt,self.T
        elif self.initype == 'NED':
            return self.N,self.E,self.D,self.Tdt,self.T
        elif self.initype == 'UTM':
            return self.Eutm,self.Nutm,self.Uutm,self.Tdt,self.T
        else:
            log.error("wrong initype")

    def __repr__(self):
        if (not hasattr(self,'X')):
            return "{},{},{},{},{}".format(self.E,self.N,self.U,self.Tdt,self.T)
        else:
            return "{},{},{},{},{}".format(self.X,self.Y,self.Z,self.Tdt,self.T)

    def XYZset(self,X=0,Y=0,Z=0,sX=0,sY=0,sZ=0):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.sX = sX
        self.sY = sY
        self.sZ = sZ

        self.initype = 'XYZ'
        self.F,self.L,self.H = conv.XYZ2GEO(self.X,self.Y,self.Z)

    def FLHset(self,F=0,L=0,H=0,sF=0,sL=0,sH=0):
        self.F = F
        self.L = L
        self.H = H
        self.sF = sF
        self.sL = sL
        self.sH = sH

        self.initype = 'FLH'
        self.X,self.Y,self.Z = conv.GEO2XYZ(self.F,self.L,self.H)
        self.sX,self.sY,self.sZ = conv.sFLH2sXYZ(F,L,H,sF,sL,sH)
        

    def ENUset(self,E=np.nan,N=np.nan,U=np.nan,
               sE=np.nan,sN=np.nan,sU=np.nan):
        self.E = E
        self.N = N
        self.U = U
        self.sE = sE
        self.sN = sN
        self.sU = sU

        self.initype = 'ENU'

    def NEDset(self,N=np.nan,E=np.nan,D=np.nan,
               sN=np.nan,sE=np.nan,sD=np.nan):
        self.N = N
        self.E = E
        self.D = D
        self.sN = sN
        self.sE = sE
        self.sD = sD

        self.initype = 'NED'


    def UTMset(self,Eutm=np.nan,Nutm=np.nan,Uutm=np.nan,
               sEutm=np.nan,sNutm=np.nan,sUutm=np.nan):
        self.Eutm = Eutm
        self.Nutm = Nutm
        self.Uutm = Uutm
        self.sEutm = sEutm
        self.sNutm = sNutm
        self.sUutm = sUutm

        self.initype = 'UTM'
        

    def add_offset(self,dA,dB,dC):
        log.warning("add_offset as method are hazardous ...")
        temp = time_series.add_offset_point(self,dA,dB,dC)
        self.__dict__ = temp.__dict__

    def Tset(self,T=0):
        if T == 0:
            T = dt.datetime.now()

        if type(T) == dt.datetime:
            self.Tdt = T
            self.T = conv.dt2posix(T)

        else:
            self.T = float(T)
            self.Tdt = conv.posix2dt(float(T))

    def ENUcalc_pt(self,refENU ):

        self.refENU = refENU

        dX =  self.X - refENU.X
        dY =  self.Y - refENU.Y
        dZ =  self.Z - refENU.Z

        Etmp,Ntmp,Utmp = conv.XYZ2ENU(dX,dY,dZ,refENU.F,refENU.L)
        self.E,self.N,self.U = Etmp[0],Ntmp[0],Utmp[0]

        if self.initype == 'FLH' and hasattr(self,'sF'):
            if not np.isnan(self.sF):
                self.sE,self.sN,self.sU = conv.sFLH2sENU(self.F,self.L,self.H,
                                                         self.sF,self.sL,self.sH)
        elif self.initype == 'XYZ' and hasattr(self,'sX'):
            if not np.isnan(self.sX):
                self.sE,self.sN,self.sU = conv.sXYZ2sENU(self.X,self.Y,self.Z,
                                                         self.sX,self.sY,self.sZ)


    def UTMcalc_pt(self,ellips="wgs84"):
        self.Eutm , self.Nutm , _ = conv.utm_geo2xy(self.F,self.L)
        self.Uutm = self.H
        

    def keysanex(self):
        return list(self.anex.keys())


    def helmert_trans(self,params='itrf2008_2_etrf2000',invert=False):
        Xb = reffram.helmert_trans(np.array([self.X,self.Y,self.Z]),params,invert)
        self.XYZset(*Xb)
        return None

    def velocity_trans(self, vx, vy, vz, epoc_init = 'auto' , epoc_end = 'auto'):
        """
        auto == epoc of the measures
        """

        if epoc_init == 'auto' and epoc_end == 'auto':
            log.error('epoc_init == auto and epoc_end == auto')
            return None

        tdt = conv.posix2dt(self.T)
        yeardec = conv.dt2year_decimal(tdt)

        if epoc_init  == 'auto':
            epoc_init = yeardec

        if epoc_end  == "auto":
            epoc_end = yeardec

        Xb = reffram.itrf_speed_calc(self.X,self.Y,self.Z, epoc_init ,
                                  vx, vy, vz, epoc_end)
        self.XYZset(*Xb)
        return None

class TimeSeriePoint:
    
    def __init__(self,stat='STAT'):
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
        Representation of a TimeSeriePoint object
        """
        if self.pts == []:
            raise Exception('ERR: TimeSeriePoint is empty ...')

        start = self.startdate()
        end = self.enddate()
        nbday = int((end - start).days + 1.)

        ratio = self.nbpts * 100. / nbday

        return '{} {} {} {} {} {} {} {:5.2f}{}'.format(self.stat,self.nbpts,'points',
                start,end,nbday,"nb days", ratio,"%")

    def __getitem__(self,i):
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


    def meta_set(self,path='',stat='STAT',name=''):
        """
        Set meta data about the TimeSerie

        Parameters
        ----------
        path : str, optional
            file path. The default is ''.
        stat : str, optional
            station 4-char. code. The default is 'STAT'.
        name : str, optional
            free name of for the TS, 
            like the experience, the periode , the software ...
            The default is ''.

        Returns
        -------
        None.
        """

        self.path = path
        self.stat = stat

        bn = os.path.basename(path)
        dn = os.path.dirname(path)

        if name == '' :
            if bn == 'tdp_final':
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


    def readfile(self,filein):
        """
        Method to read the data form a file
        Should be used with care

        Parameters
        ----------
        filein : str
            path of the file.

        Returns
        -------
        None.

        """
        self.__dict__ = files_rw.read_all_points(filein).__dict__

        self.interp_set()

    def add_point(self,inPoint):
        """
        Method to add a Point in the TimeSerie Object

        Parameters
        ----------
        inPoint : Point Object

        Returns
        -------
        None.

        """
        self.pts.append(inPoint)
        # this line is discontiued, because now nbpts is a property
        #self.nbpts = len(self.pts)

        self.bool_interp_uptodate = False

    def aleapt(self):
        """
        Method to get a random Point in the TimeSeries

        Returns
        -------
        Point Object

        """
        ipt = np.random.randint(self.nbpts)
        log.info("random point selected #%s",str(ipt))

        return self.pts[ipt]

    def startdate(self):
        """
        Method to get the first epoch of the data in the TimeSerie

        Returns
        -------
        DateTime

        """
        self.sort()
        return conv.posix2dt(self.pts[0].T)

    def enddate(self):
        """
        Method to get the last epoch of the data in the TimeSerie

        Returns
        -------
        DateTime

        """
        self.sort()
        return conv.posix2dt(self.pts[-1].T)
    
    
    def len_period(self,output_seconds=False):
        """
        Method to get the period length

        Returns
        -------
        timedelta or 
        """    
        delta = self.enddate() - self.startdate()
        
        if output_seconds:
            return delta.days*86400 + delta.seconds
        else:
            return delta

    def interval_nominal(self):
        """
        Method to get the nominal internal between two epochs.

        Returns
        -------
        float
            interval nominal.

        """

        if len(self.pts) < 2:
            self.i_nomi = 0
        else:
            Ttemp = np.sort([pt.T for pt in self.pts])
            self.i_nomi = np.round(np.min(np.diff(Ttemp)),1)

        return self.i_nomi


    def from_list(self,T,A,B,C,coortype='XYZ',sA=[],sB=[],sC=[]):
        """
        Method to load data from lists to the TimeSerie

        Parameters
        ----------
        T : float
            Time.
        A : list of float
            X, F (latitude), E..
        B : list of float
            Y, L (longitude), N.
        C : list of float
            Z, H (hight), U.
        coortype : str, optional
            The coordinates type. The default is 'XYZ'.
        sA : list of float, optional
            sigma of A component. The default is [].
        sB : list of float, optional
            sigma of B component. The default is [].
        sC : list of float, optional
            sigma of C component. The default is [].

        Returns
        -------
        None.

        """

        if not sA:
            sA = np.zeros(len(A))
        if not sB:
            sB = np.zeros(len(B))
        if not sC:
            sC = np.zeros(len(C))

        for t,a,b,c,sa,sb,sc in zip(T,A,B,C,sA,sB,sC):
            point = Point(a,b,c,t,coortype,sa,sb,sc)

            self.add_point(point)

        if coortype == "ENU":
            self.boolENU = True

        if coortype == "UTM":
            self.boolUTM = True

        self.sort()

        return None


    def to_list(self,coortype='XYZ',specific_output=None,
                time_as_datetime=False):
        """
        Export the TimeSerie Object as Lists (Numpy Arrays)

        Parameters
        ----------
        coortype : str, optional
        The coordinates type exported to the list.
        The default is 'XYZ'.
        specific_output : int, optional
            ask for a specific list, ranges between 0 and 6.
            The default is None.
        time_as_datetime : bool, optional
            if True the Time list is exported in datetime
            if False the Time list is exported in Posix time

        Returns
        -------
        A,B,C,T,sA,sB,sC : lists
            A = X, F (latitude), E.
            B = Y, L (longitude), N.
            C = Z, H (hight), U.
            T = Time
            sA = sigma of A component
            sB = sigma of B component
            sC = sigma of C component
        """


        if coortype == 'XYZ':
            A,B,C = 'X','Y','Z'
            sA,sB,sC = 'sX','sY','sZ'

        elif coortype == 'FLH':
            A,B,C = 'F','L','H'
            sA,sB,sC = 'sF','sL','sH'

        elif coortype == 'ENU':
            if self.boolENU == False:
                log.warning("no ENU coord. for " + self.name)
                return None

            A,B,C = 'E','N','U'
            sA,sB,sC = 'sE','sN','sU'

        elif coortype == 'UTM':
            if self.boolUTM == False:
                log.warning("no UTM coord. for " + self.name)
                return None

            A,B,C = 'Eutm','Nutm','Uutm'
            sA,sB,sC = 'sEutm','sNutm','sUutm'

        else:
            log.error("coortype does not exist")

        if self.nbpts == 0:
            log.error(self.name + " the timeserie is empty")

        A = np.asarray([ getattr(pt,A) for pt in self.pts ])
        B = np.asarray([ getattr(pt,B) for pt in self.pts ])
        C = np.asarray([ getattr(pt,C) for pt in self.pts ])
        T = np.asarray([ pt.T for pt in self.pts ])

        if hasattr(self.pts[0],sA):
            sA = np.asarray([ getattr(pt,sA) for pt in self.pts ])
            sB = np.asarray([ getattr(pt,sB) for pt in self.pts ])
            sC = np.asarray([ getattr(pt,sC) for pt in self.pts ])
        else:
            sA = np.asarray([ np.nan ] * len(self.pts))
            sB = np.asarray([ np.nan ] * len(self.pts))
            sC = np.asarray([ np.nan ] * len(self.pts))
        # il faut squeezer les vecteurs parce que des fois on se retrouve
        # avec des matrices
        
        if time_as_datetime:
            Tout = conv.posix2dt(T)
        else:
            Tout = T
        
        sq = np.squeeze
        outtup = (sq(A) , sq(B) , sq(C) , sq(Tout) , sq(sA) , sq(sB) , sq(sC))
        if specific_output == None:
            return outtup
        elif type(specific_output) is int:
            return outtup[specific_output]
        else:
            log.info("this mode must be implemented ;)")
            log.info("use an int as index instead")
            return outtup
        
    def to_dataframe(self,coortype ='XYZ'):
        """
        Export the TimeSerie Object as DataFrame
        
        Parameters
        ----------
        coortype : str or iterable of str.
            The coordinates type exported to the DataFrame.
            'XYZ', 'FLH', 'ENU', 'NED'
            can be also an iterable like ('XYZ','FLH')
            The default is 'XYZ'.

        Returns
        -------
        DF : DataFrame
            output DataFrame.
        """
        
        if not utils.is_iterable(coortype):
            coortype = (coortype,)
        
        ColStk = tuple()
        ColNameStk = []
        
        for icoty , coty in enumerate(coortype):
            A,B,C,T,sA,sB,sC = self.to_list(coty)
            
            if coty == "UTM":
                cotycolnam = ["Eutm","Nutm","Uutm"]
            else:
                cotycolnam = coty
                
            if icoty == 0:
                Tdt = conv.posix2dt(T)
                ColStk = ColStk + (Tdt,T,A,B,C,sA,sB,sC)  
                ColNameStk = ["Tdt","T"] + [e for e in cotycolnam] + ["s" + e for e in cotycolnam]
            else:
                ColStk = ColStk + (A,B,C,sA,sB,sC)
                ColNameStk = [e for e in coty] + ["s" + e for e in coty]
                
        BIG = np.column_stack(ColStk)
        DF = pd.DataFrame(BIG)
        DF.columns = ColNameStk
        DF = DF.infer_objects()
        
        return DF

    def sort(self):
        """
        Internal method to sort the point in the TimeSerie Object

        Returns
        -------
        None.

        """
        self.pts.sort(key=lambda x: x.T)

    def plot(self,coortype='ENU',
             diapt=2,
             alpha=0.8,
             fig=1,
             errbar=True,
             new_style=True,
             symbol = '.',
             errbar_width=1,
             ylim=None):
        """
        Plot data in a TimeSerie Object

        Parameters
        ----------
        coortype : str, optional
            The coordinates type. The default is 'ENU'.
        diapt : float, optional
            Point diamaeter. The default is 2.
        alpha : float, optional
            Alpha (transparency) of points. The default is 0.8.
        fig : int or Figure object, optional
            Figure ID where the data will be plotted
            can accept a int (id of a Figure)
            OR the figure Object itself.
            The default is 1.
        errbar : bool, optional
            Plot the error bars. The default is True.
        new_style : bool, optional
            Plot in a new style.
            The old style is only kept for legacy
            The default is True.
        symbol : str, optional
            symbol. The default is '.'.
        errbar_width : TYPE, optional
            coefficient for the error bar size. The default is 1.

        Returns
        -------
        The matplotlib Figure object.

        """
        
        log.setLevel(logging.INFO)

        if new_style:
            styleint = 310
        else:
            styleint = 220

        try:
            A, B, C , T , sA, sB, sC = self.to_list(coortype=coortype)
        except TypeError as tyer:
            log.error("unable to get coordinates")
            log.info("TRICK : check if the given coortype is in the timeserie")
            raise tyer

        if coortype == 'ENU':
            Atitle = 'East'
            Btitle = 'North'
            Ctitle = 'Up'
            ABtitle = 'East North'
            yylabel = 'displacement (m)'

        elif coortype == 'XYZ':
            Atitle = 'X'
            Btitle = 'Y'
            Ctitle = 'Z'
            yylabel = 'displacement (m)'
            ABtitle = 'X Y (sans signification)'

        elif coortype == 'FLH':
            Atitle = 'Phi'
            Btitle = 'Lambda'
            Ctitle = 'Haut'
            yylabel = 'displacement (m)'
            ABtitle = 'Phi Lambda (sans signification)'


        elif coortype == 'UTM':
            Atitle = 'East (UTM)'
            Btitle = 'North (UTM)'
            Ctitle = 'Up'
            yylabel = 'displacement (m)'
            ABtitle = 'East North (UTM)'

        else:
            Atitle = 'A'
            Btitle = 'B'
            Ctitle = 'C'
            yylabel = 'displacement (??)'
            ABtitle = 'A & B'

        log.info("plot : %s, pts : %s", self.nbpts, self.stat)

        namest=0
        namend=10

        Tdt = conv.posix2dt(T)

        if type(fig) is int:
            figobj = plt.figure(fig)
        elif type(fig) is plt.Figure:
            figobj = fig

        figobj.suptitle(self.stat)


        if self.name:
            name4plot = self.name[namest:namend]
        else:
            name4plot = self.stat

        plt.subplot(styleint+1)
        
        if errbar:
            plt.errorbar(Tdt,A,sA, fmt=symbol,label=name4plot,
                         markersize=diapt, alpha=alpha,ecolor='xkcd:light grey',
                         elinewidth=errbar_width)
        else:
            plt.plot(Tdt,A,symbol,label=name4plot,
                     markersize=diapt, alpha=alpha)
        try:
            plt.legend()
        except:
            pass
        #plt.xlabel('Date')
        plt.ylabel(yylabel)
        plt.title(Btitle)
        plt.title(Atitle)
                
#        ax = plt.gca()

#        if coortype == 'ENU':
##            refstr = 'ref XYZ = ' + utils.join_improved(',',self.refENU.X ,
##                                                          self.refENU.Y ,
##                                                          self.refENU.Z)
        plt.subplot(styleint+2)
        if errbar:
            plt.errorbar(Tdt,B,sB, fmt=symbol,label=name4plot,
                         markersize=diapt, alpha=alpha,ecolor='xkcd:light grey',
                         elinewidth=errbar_width)
        else:
            plt.plot(Tdt,B,symbol,label=name4plot,
                     markersize=diapt, alpha=alpha)
        try:
            plt.legend()
        except:
            pass
        #plt.xlabel('Date')
        plt.ylabel(yylabel)
        plt.title(Btitle)

        plt.subplot(styleint+3)
        if errbar:
            plt.errorbar(Tdt,C,sC, fmt=symbol,label=name4plot,
                         markersize=diapt, alpha=alpha,ecolor='xkcd:light grey',
                         elinewidth=errbar_width)
        else:
            plt.plot(Tdt,C,symbol,label=name4plot,
                     markersize=diapt, alpha=alpha)
        try:
            plt.legend()
        except:
            pass
        
        if ylim:
            [a.set_ylim(ylim) for a in figobj.axes]
        
        plt.xlabel('Date')
        plt.ylabel(yylabel)
        plt.title(Ctitle)
        figobj.autofmt_xdate()
        figobj.set_size_inches(8.27,11.69)
        figobj.tight_layout()
        plt.subplots_adjust(top=0.93)

        if not new_style:
            plt.subplot(styleint+4)
            plt.axis('equal')
            try:
                plt.legend()
            except:
                pass
            plt.plot(A,B,'.',label=name4plot,
                     markersize=diapt, alpha=alpha)
            plt.xlabel(Atitle + ' ' + yylabel)
            plt.ylabel(Btitle + ' ' + yylabel)
            plt.title(ABtitle)

        return figobj


    def plot_discont(self,fig=1):
        """
        Plot discontinuties of a TimeSerie Object contained in discont list

        Parameters
        ----------
        fig : int or Figure object, optional
            Figure ID where the data will be plotted
            can accept a int (id of a Figure)
            OR the figure Object itself.
            The default is 1.

        Returns
        -------
        None.

        """

        if not self.bool_discont:
            log.warning("no discontinuities in TimeSerie")
            return None

        if type(fig) is int:
            figobj = plt.figure(fig)
        elif type(fig) is plt.Figure:
            figobj = fig

        for ax in figobj.axes:
            stats.plot_vertical_bar_ax(self.discont,ax,"r")

        if self.bool_discont_manu:
            for ax in figobj.axes:
                stats.plot_vertical_bar_ax(self.discont_manu,ax,"g")

#        figobj.axes[1]
#        stats.plot_vertical_bar(self.discont)
#
#        figobj.axes[2]
#        stats.plot_vertical_bar(self.discont)


    def discont_manu_click(self,fig=1):
        """
        manual discontinuities are both recorded in the "main" discont list
        and in a new discont_manu list,
        thus the manual discontinuites can be identified

        IMPORTANT : cursor objects (multi , cid)
                    must be stored as global variables like this :
                    multi , cid = tsout.discont_manu_click()


        NOTE : This method was created before point_n_click_plot():
            this other one is more complete
            both has to be merged ASAP !!!!!

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
            ix, iy = matplotlib.dates.num2date(event.xdata).replace(tzinfo=None), event.ydata

            log.info("discontinuity recorded : %s" , ix)
            for ax in figobj.axes:
                stats.plot_vertical_bar_ax([ix],ax,"g")
#

            self.bool_discont      = True
            self.bool_discont_manu = True

            self.discont.append(ix)
            self.discont = sorted(self.discont)

            self.discont_manu.append(ix)
            self.discont_manu = sorted(self.discont_manu)

            #figobj.show()
            plt.draw()

            return None

        multi = MultiCursor(figobj.canvas, figobj.axes , color='k', lw=1)
        cid   = figobj.canvas.mpl_connect('key_press_event', onclick_discont)

        return multi , cid

    def initype(self):
        L = [ pt.initype for pt in self.pts ]
        return Counter(L).most_common(1)[0][0]

    def ENUcalc(self, refENU ):
        """
        Method to determine the ENU components based on a reference point

        Parameters
        ----------
        refENU : Point Object or TimeSeriePoint Object
            Reference point.

        Returns
        -------
        None.

        """
        if refENU.__class__.__name__ == 'Point':
            self.refENU = refENU
            [ pt.ENUcalc_pt(refENU) for pt in self.pts ]
            self.boolENU = True
            self.bool_interp_uptodate = False
            self.interp_set()
        elif refENU.__class__.__name__ == 'TimeSeriePoint':
            self.refENU = refENU.mean_posi()
            refENU.interp_set()
            [ pt.ENUcalc_pt(Point(A=refENU.XfT(pt.T),
                                  B=refENU.YfT(pt.T),
                                  C=refENU.ZfT(pt.T),
                                  initype='XYZ',
                                  T=pt.T)) for pt in self.pts ]

            self.boolENU = True
            self.bool_interp_uptodate = False
            self.interp_set()

    def ENUcalc_from_mean_posi(self,mean_type="median"):
        """
        Method to determine the ENU components based directly 
        on the mean/median position

        Returns
        -------
        None.

        """
        self.ENUcalc(self.mean_posi(mean_type="median"))
        return None
    
    
    def ENUcalc_from_first_posi(self):
        """
        Method to determine the ENU components based directly
        on the mean/median position

        Returns
        -------
        None.

        """
        self.ENUcalc(self.pts[0])
        return None

    def from_uniq_point(self,Point,startdate,enddate,pas=1):
        self.del_data()

        if type(startdate) == dt.datetime:
            startdate = conv.dt2posix(startdate)
        if type(enddate) == dt.datetime:
            enddate = conv.dt2posix(enddate)

        N = int(np.round( enddate - startdate / pas ))
        #datelist = np.arange(startdate,enddate,pas)

        Point.Tset(startdate)

        for i in range(N):
            Point.Tset(Point.T + pas)
            self.add_point(copy.copy(Point))

        self.i_nomi = pas
        
        
    def UTMcalc(self):
        """
        Method to determine the UTM E and N projected coordinates

        Parameters
        ----------

        Returns
        -------
        None.

        """
        self.boolUTM = True
        [ pt.UTMcalc_pt() for pt in self.pts ]
        

    def time_win(self,windows,mode='keep'):
        '''IL EST TRES DANGEREUX DE L'APPLIQUER UN FENETRAGE A SOI MEME'''
        self.__dict__ = time_series.time_win(self,windows,mode).__dict__


    def interp_set(self,interptype = 'slinear'):
        """
        Method to set the coordinate interpolators

        Parameters
        ----------
        interptype : TYPE, optional
            Interpolation type. The default is 'slinear'.

        Returns
        -------
        None.

        """
        if (not hasattr(self.pts[0],'E')) or np.isnan(self.pts[0].E) == True:
            log.warning("no ENU for " + self.name)
        else:
            E,N,U,T,_,_,_ = self.to_list('ENU')
            self.EfT = scipy.interpolate.interp1d(T,E,bounds_error=False,kind=interptype)
            self.NfT = scipy.interpolate.interp1d(T,N,bounds_error=False,kind=interptype)
            self.UfT = scipy.interpolate.interp1d(T,U,bounds_error=False,kind=interptype)

        if (not hasattr(self.pts[0],'X')) or np.isnan(self.pts[0].X) == True:
            log.warning("no XYZ for " + self.name)
        else:

            X,Y,Z,T,_,_,_ = self.to_list('XYZ')

            self.XfT = scipy.interpolate.interp1d(T,X,bounds_error=False,kind=interptype)
            self.YfT = scipy.interpolate.interp1d(T,Y,bounds_error=False,kind=interptype)
            self.ZfT = scipy.interpolate.interp1d(T,Z,bounds_error=False,kind=interptype)


        if (not hasattr(self.pts[0],'F')) or np.isnan(self.pts[0].L) == True:
            log.warning("no FLH for " + self.name)
        else:
            F,L,H,T,_,_,_ = self.to_list('FLH')

            self.FfT = scipy.interpolate.interp1d(T,F,bounds_error=False,kind=interptype)
            self.LfT = scipy.interpolate.interp1d(T,L,bounds_error=False,kind=interptype)
            self.HfT = scipy.interpolate.interp1d(T,H,bounds_error=False,kind=interptype)

        if (not hasattr(self.pts[0],'Eutm')) or np.isnan(self.pts[0].Eutm) == True:
            log.warning("no UTM for " + self.name)
        else:
            Eutm,Nutm,Uutm,T,_,_,_ = self.to_list('UTM')

            self.EutmfT = scipy.interpolate.interp1d(T,Eutm,bounds_error=False,kind=interptype)
            self.NutmfT = scipy.interpolate.interp1d(T,Nutm,bounds_error=False,kind=interptype)
            self.UutmfT = scipy.interpolate.interp1d(T,Uutm,bounds_error=False,kind=interptype)

        self.bool_interp_uptodate = True

    def interp_get(self,T,coortype='ENU'):
        """
        Method to get the coordinate interpolators

        Parameters
        ----------
        T : float or list of float
            Time (IN POSIX Time) where the interpolation is wished.
        coortype : str, optional
            The coordinates type. The default is 'ENU'.

        Returns
        -------
        tsout : 
            DESCRIPTION.

        """

        if self.bool_interp_uptodate == False:
            log.warning("interp obsolete, recalcul auto")
            self.interp_set()

        tsout = copy.copy(self)
        tsout.del_data()

        if not utils.is_iterable(T):
            T = np.array([T])

        if coortype == 'ENU':
            A = self.EfT(T)
            B = self.NfT(T)
            C = self.UfT(T)

        if coortype == 'XYZ':
            A = self.XfT(T)
            B = self.YfT(T)
            C = self.ZfT(T)

        if coortype == 'FLH':
            A = self.FfT(T)
            B = self.LfT(T)
            C = self.HfT(T)

        if coortype == 'UTM':
            A = self.EutmfT(T)
            B = self.NutmfT(T)
            C = self.UutmfT(T)

        for i in range(len(T)):
            tsout.add_point(Point(A[i],B[i],C[i],T=T[i],initype=coortype))

        return tsout

    def set_discont(self,indiscont):
        """
        Method to set the discontinuties list

        Parameters
        ----------
        indiscont : list of time
            Discontinuities in the TimeSerie.

        Returns
        -------
        None.

        """
        self.discont = indiscont
        self.bool_discont = True

    def mean_posi(self,coortype='XYZ',outtype='point',mean_type='median'):
        """
        Method to determine the mean position of the TimeSerie

        Parameters
        ----------
        coortype : TYPE, optional
            The coordinates type. The default is 'XYZ'.
        outtype : TYPE, optional
            'point' or 'tuple'. The default is 'point'.
        mean_type : TYPE, optional
            'mean' or 'median'. The default is 'median'.

        Returns
        -------
        Point or coordinates tuple

        """

        #special case where only one point
        if self.nbpts == 1:
            return copy.copy(self)

        A,B,C,T,sA,sB,sC = self.to_list(coortype=coortype)


        if mean_type == 'mean':
            Aout = np.nanmean(A)
            Bout = np.nanmean(B)
            Cout = np.nanmean(C)
        elif mean_type == 'median':
            Aout = np.nanmedian(A)
            Bout = np.nanmedian(B)
            Cout = np.nanmedian(C)

        Tout = (np.max(T) - np.min(T)) / 2

        if outtype == "point":
            out = Point(Aout,Bout,Cout,Tout)
        elif outtype == "tuple":
            out = (Aout,Bout,Cout,Tout)
        else:
            out = (Aout,Bout,Cout,Tout)

        return out

    def add_offset(self,dA,dB,dC):
        """
        NOTE 160415 : add_offset as method are hazardous ...
        use fct add_offset_ts instead
        """
        log.warning("add_offset as method are hazardous ...")
        for pt in self.pts:
            pt.add_offset(dA,dB,dC)

    def decimate(self,dec):
        """
        Method to decimate a TimeSerie

        Parameters
        ----------
        dec : int
            keep 1/dec point in the TimeSerie.

        Returns
        -------
        None.

        """
        self.__dict__ = time_series.decimate_cleaner(self,dec).__dict__
        #decimate_cleaner(self,dec,True)

    def find_point(self,tin,tol=0.001,stop_when_found=True):
        """
        Method to find a specific point according to its timestamp

        Parameters
        ----------
        tin : float or datetime
            timestamp of the researched point.
        tol : float, optional
            tolerence of the research. The default is 0.001.
        stop_when_found : bool, optional
            Stop the research when a point is found. The default is True.

        Returns
        -------
        Point Object
            Point Found.
        int or list of int
            index of the point.

        """

        find = False

        if type(tin) is dt.datetime:
            tin = conv.dt2posix(tin)

        tmin = self.pts[0].T
        tmax  = self.pts[-1].T

        if not ( tmin <= tin <= tmax ):
            log.warning("t !€ [tmin , tmax]")

        if stop_when_found:
            for i,p in enumerate(self.pts):
                if np.abs(p.T - tin) < tol:
                    find = True
                    break
            if find:
                return p,i
            else:
                return Point(),np.nan

        else:
            pts_stk = []
            i_stk = []
            for i,p in enumerate(self.pts):
                if np.abs(p.T - tin) < tol:
                    pts_stk.append(p)
                    i_stk.append(i)
            return pts_stk  , i_stk
        
    def remove_duplicate_pts(self,coortype="XYZ"):
        T = self.to_dataframe(coortype)["T"] 
        
        dup_bool = T.duplicated()
        
        if dup_bool.sum() > 0:
            log.warn("%s duplicated point(s) removed for %s",
                     dup_bool.sum(), self.name)
            
        self.pts = list(pd.Series(self.pts)[np.logical_not(dup_bool)])
    
    
    


 #  ______                      _                      _        _    _____ _                         
 # |  ____|                    (_)                    | |      | |  / ____| |                        
 # | |__  __  ___ __   ___ _ __ _ _ __ ___   ___ _ __ | |_ __ _| | | |    | | __ _ ___ ___  ___  ___ 
 # |  __| \ \/ / '_ \ / _ \ '__| | '_ ` _ \ / _ \ '_ \| __/ _` | | | |    | |/ _` / __/ __|/ _ \/ __|
 # | |____ >  <| |_) |  __/ |  | | | | | | |  __/ | | | || (_| | | | |____| | (_| \__ \__ \  __/\__ \
 # |______/_/\_\ .__/ \___|_|  |_|_| |_| |_|\___|_| |_|\__\__,_|_|  \_____|_|\__,_|___/___/\___||___/
 #             | |                                                                                   
 #             |_|                                                                                   


    
    
class Attitude:
    def __init__(self,R=0,P=0,Y=0,T=0,sR=0,sP=0,sY=0,devID='NULL',angtype='deg'):

        self.Tset(T)
        self.devID = devID

        if angtype == 'deg':
            self.R = R
            self.P = P
            self.Y = Y

        elif angtype == 'rad' :
            self.R = np.rad2deg(R)
            self.P = np.rad2deg(P)
            self.Y = np.rad2deg(Y)
        else:
            raise Exception("Mauvais angtype")

        self.Qcalc()

    def __call__(self):
        return self.R,self.P,self.Y,self.Tdt

    def __repr__(self):
        return "{},{},{},{}".format(self.R,self.P,self.Y,self.Tdt)


    def Tset(self,T=0):
        if T == 0:
            T = dt.datetime.now()

        if type(T) == dt.datetime:
            self.Tdt = T
            self.T = conv.dt2posix(T)

        else:
            self.T = float(T)
            self.Tdt = conv.posix2dt(float(T))

    def RPYget(self):
        return self.R,self.P,self.Y

    def RPYset(self,R=0,P=0,Y=0,sR=0,sP=0,sY=0):
        self.R  = R
        self.P  = P
        self.Y  = Y
        self.sR  = sR
        self.sP  = sP
        self.sY  = sY

    def Qcalc(self):
        self.Q = conv.quaternion(self.R , self.P , self.Y , 'deg')
        return None


    
class TimeSerieObs(object):

    ''' LES DIFFERENCES AVEC TSPOINT
        * Les objets ne contiennent qu'un type de données sous une seul forme
        (a la difference d'un point qui peut exister sous plusieurs formes)
        * Dans un fichier en input, il peut y avoir plusieurs "devices"
          => les fonctions de lectures produisent donc obligatoirement des listes
          de TS (le cas échéant une liste à 1 élt)
          => la methode readfile() nécessite donc l'indice de la device'''

    def __init__(self,typeobs='NULL',filepath =''):

        self.obs = []
        self.nbobs = 0
        self.i_nomi = 0
        self.typeobs = typeobs
        # ASM Attitude
        self.bool_interp_uptodate = False

        self.meta_set(filepath)

    def meta_set(self,path='',devID='NULL',name=''):
        self.path = path
        self.devID = devID

        bn = os.path.basename(path)
        dn = os.path.dirname(path)

        if name == '' :
            if bn == 'tdp_final':
                self.name = os.path.basename(dn)
            else:
                self.name = bn

        self.interval_nominal()

    def del_data(self):
        self.obs = []
        self.nbobs = 0

    def readfile(self,filein,indtab=0):
        log.info("selection device %s" , indtab)
        temp = files_rw.read_all_obs(filein)[indtab]
        self.__dict__ = temp.__dict__

        self.interp_set()

    def add_obs(self,inObs):
        self.obs.append(inObs)
        self.nbobs = len(self.obs)

        self.bool_interp_uptodate = False

    def aleaobs(self):
        iobs = randrange(self.nbobs)
        log.info("observation no " + str(iobs))

        log.info(self.obs[iobs])

        return self.obs[iobs]


    def interval_nominal(self):
        if len(self.obs) < 2:
            self.i_nomi = 0
        else:
            Ttemp = np.sort([o.T for o in self.obs])
            self.i_nomi = np.round(np.min(np.diff(Ttemp)),1)

        return self.i_nomi

    def timewin(self,windows,mode='keep'):
        '''IL EST TRES DANGEREUX DE L'APPLIQUER UN FENETRAGE A SOI MEME'''
        self.__dict__ = time_series.time_win(self,windows,mode).__dict__

    def startdate(self):
        return self.obs[0].T

    def enddate(self):
        return self.obs[-1].T

    def to_list(self):
        if self.typeobs == 'NULL':
            log.error("pas de typeobs defini (NULL)")
            return 0

        if self.typeobs == 'RPY':
            A,B,C = 'R','P','Y'
            sA,sB,sC = 'sR','sP','sY'

        A = np.asarray([ getattr(o,A) for o in self.obs ])
        B = np.asarray([ getattr(o,B) for o in self.obs ])
        C = np.asarray([ getattr(o,C) for o in self.obs ])
        T = np.asarray([ o.T for o in self.obs ])

        if hasattr(self,sA):
            sA = np.asarray([ getattr(o,sA) for o in self.obs ])
            sB = np.asarray([ getattr(o,sB) for o in self.obs ])
            sC = np.asarray([ getattr(o,sC) for o in self.obs ])

        else:
            sA = np.asarray([ np.nan ] * len(self.obs))
            sB = np.asarray([ np.nan ] * len(self.obs))
            sC = np.asarray([ np.nan ] * len(self.obs))

        return A, B, C , T , sA , sB , sC



    def interp_set(self,interptype='slinear'):

        if self.typeobs == 'NULL':
            log.error("no typeobs defined (NULL)")
            return 0

        if self.typeobs == 'RPY':

            R,P,Y,T,_,_,_ = self.to_list()

            self.RfT = scipy.interpolate.interp1d(T,R,bounds_error=False,kind=interptype)
            self.PfT = scipy.interpolate.interp1d(T,P,bounds_error=False,kind=interptype)
            self.YfT = scipy.interpolate.interp1d_ang(T,Y,bounds_error=False,kind=interptype)

        self.bool_interp_uptodate = True

    def interp_get(self,T):

        if self.bool_interp_uptodate == False:
            log.warning("interp obsolete, recalcul auto")
            self.interp_set()

        tsout = copy.copy(self)
        tsout.del_data()

        if not utils.is_iterable(T):
            T = np.array([T])


        if self.typeobs == 'NULL':
            log.error("no typeobs defined (NULL)")
            return 0

        if self.typeobs == 'RPY':
            A = self.RfT(T)
            B = self.PfT(T)
            C = self.YfT(T)

            for i in range(len(T)):
                    tsout.add_obs(Attitude(A[i],B[i],C[i],T=T[i]))

        return tsout

    def plot(self,diapt=10,alpha=0.8,fig=1,new_style=True):

        A, B, C , T , sA, sB, sC = self.to_list()

        log.info("plot : %s, pts : %s", self.nbpts, self.stat)

        if self.typeobs == 'RPY':
            listtitle = ['','Roll','Pitch','Yaw']
        else:
            listtitle = ['','','','']

        namest=0
        namend=10
        Tdt = conv.posix2dt(T)

        if new_style:
            styleint = 410
        else:
            styleint = 220

        plt.figure(fig)
        plt.subplot(styleint + 1)
        plt.axis('equal')
        plt.plot(A,B,'.',label=str(self.devID)[namest:namend], markersize=diapt, alpha=alpha)
        plt.legend()
        plt.title(listtitle[0])

        plt.subplot(styleint + 2)
        plt.plot(Tdt,A,'.',label=str(self.devID)[namest:namend], markersize=diapt, alpha=alpha)
        plt.legend()
        plt.title(listtitle[1])

        plt.subplot(styleint + 3)
        plt.plot(Tdt,B,'.',label=str(self.devID)[namest:namend], markersize=diapt, alpha=alpha)
        plt.legend()
        plt.title(listtitle[2])

        plt.subplot(styleint + 4)
        plt.plot(Tdt,C,'.',label=str(self.devID)[namest:namend], markersize=diapt, alpha=alpha)
        plt.legend()
        plt.title(listtitle[3])


#  _____       _                      _   _             _____  _       _
# |_   _|     | |                    | | (_)           |  __ \| |     | |
#   | |  _ __ | |_ ___ _ __ __ _  ___| |_ ___   _____  | |__) | | ___ | |_
#   | | | '_ \| __/ _ \ '__/ _` |/ __| __| \ \ / / _ \ |  ___/| |/ _ \| __|
#  _| |_| | | | ||  __/ | | (_| | (__| |_| |\ V /  __/ | |    | | (_) | |_
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


class point_n_click_plot():
    """
    This method allow to do "point and click" on a plot, to localize offsets 
    for instance
    
    Usage
    -----
    Data have to be ploted already in a figure
    
    .. code-block:: python
    
        PnC = point_n_click_plot()
        multi , cid = PnC(fig=1,Xdata_are_time=True)
        PnC.selectedX

    i.e.
    
    Create an object point_n_click_plot (here it is PnC in the exemple below) \n
    Call the object like a function with as 1st argument the id of the plot figure or the plot figure itself \n
    Make your selection using the SPACE key \n
    Get your results in a list called PnC.selectedX   
    
    Important
    ---------
    
    cursor objects (i.e. multi & cid)
    must be stored as global variables when you call the method
    like this :

    .. code-block:: python
    
        multi , cid = PnC(fig=1)
    
    """
    



    def __init__(self):
        self.selectedX = []
        self.ver_bar_stk = []


    def __call__(self,fig=1,Xdata_are_time=True):
        """
        IMPORTANT : cursor objects (i.e. multi & cid)
                    must be stored as global variables when you call the method
                    like this :
                    multi , cid = PnC(fig=1)
        """

        if type(fig) is int:
            figobj = plt.figure(fig)
        elif type(fig) is matplotlib.figure.Figure:
            figobj = fig

        log.info("press SPACE to record a X-value, \n       press R to Remove the previously recorded one")


        def onclick_discont(event):

            if event.key == ' ':
                if Xdata_are_time:
                    ix, iy = matplotlib.dates.num2date(event.xdata).replace(tzinfo=None),event.ydata
                else:
                    ix, iy = event.xdata , event.ydata

                log.info("X value recorded : " , ix)

                for ax in figobj.axes:
                    out_bar_list = stats.plot_vertical_bar_ax([ix],ax,"b",
                                                             linewidth=1)

                    self.ver_bar_stk.append(out_bar_list[0])

                self.selectedX.append(ix)

                plt.draw()

            elif event.key in ('r','R') and len(self.selectedX) > 0:
                last = self.selectedX[-1]

                self.selectedX.remove(last)

                last_bars = self.ver_bar_stk[-len(figobj.axes):]
                
                for bar in last_bars:
                    self.ver_bar_stk.remove(bar)    
                    bar.remove()
                    
                log.info("value removed : " , last )

                plt.draw()

            return None

        multi = MultiCursor(figobj.canvas, figobj.axes , color='k', lw=1)
        cid   = figobj.canvas.mpl_connect('key_press_event', onclick_discont)

        return multi , cid
