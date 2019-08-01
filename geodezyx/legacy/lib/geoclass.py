
 # -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 10:40:16 2014

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

import datetime as dt
import numpy as np
import dateutil.parser
import re
import glob
import scipy
from scipy.interpolate import interp1d
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.figure import Figure
import csv
import pandas as pd
from collections import Counter
import copy
import warnings
import time
from random import randrange
import operator
from natsort import natsorted, ns
import tabulate
import geo_files_converter_lib as gfc
import geo_trop as gtro
from matplotlib.widgets import MultiCursor
import matplotlib
import linecache
import math

import genefun as genefun
import softs_runner
import geodetik as geok


#import combi_mgex as cmg



dUTCGPS = 16

class Point():
    def __init__(self,A=0.,B=0.,C=0.,T=0.,initype='XYZ',
                 sA=0.,sB=0.,sC=0.,name='noname'):

        self.Tset(T)
        self.ENUset()
        self.name = name
        self.anex = dict()

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
        else:
            print("ERR : init point : mauvais initype")

    def __call__(self):

        if self.initype == 'XYZ':
            return self.X,self.Y,self.Z,self.Tdt,self.T
        elif self.initype == 'FLH':
            return self.F,self.L,self.H,self.Tdt,self.T
        elif self.initype == 'ENU':
            return self.E,self.N,self.U,self.Tdt,self.T
        elif self.initype == 'NED':
            return self.N,self.E,self.D,self.Tdt,self.T
        else:
            print("ERR : call point : mauvais initype")

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
        self.F,self.L,self.H = geok.XYZ2GEO(self.X,self.Y,self.Z)

    def FLHset(self,F=0,L=0,H=0,sF=0,sL=0,sH=0):
        self.F = F
        self.L = L
        self.H = H
        self.sF = sF
        self.sL = sL
        self.sH = sH

        self.initype = 'FLH'
        self.X,self.Y,self.Z = geok.GEO2XYZ(self.F,self.L,self.H)
        self.sX,self.sY,self.sZ = geok.sFLH2sXYZ(F,L,H,sF,sL,sH)

    def ENUset(self,E=np.nan,N=np.nan,U=np.nan,sE=np.nan,sN=np.nan,sU=np.nan):
        self.E = E
        self.N = N
        self.U = U
        self.sE = sE
        self.sN = sN
        self.sU = sU

        self.initype = 'ENU'

    def NEDset(self,N=np.nan,E=np.nan,D=np.nan,sN=np.nan,sE=np.nan,sD=np.nan):
        self.N = N
        self.E = E
        self.D = D
        self.sN = sN
        self.sE = sE
        self.sD = sD

        self.initype = 'NED'

    def add_offset(self,dA,dB,dC):
        print("NOTE 160415 : add_offset as method are hazardous ...")
        temp = add_offset_point(self,dA,dB,dC)
        self.__dict__ = temp.__dict__

    def Tset(self,T=0):
        if T == 0:
            T = dt.datetime.now()

        if type(T) == dt.datetime:
            self.Tdt = T
            self.T = geok.dt2posix(T)

        else:
            self.T = float(T)
            self.Tdt = geok.posix2dt(float(T))

    def ENUcalc_pt(self,refENU ):

        self.refENU = refENU

        dX =  self.X - refENU.X
        dY =  self.Y - refENU.Y
        dZ =  self.Z - refENU.Z

        self.E,self.N,self.U = geok.XYZ2ENU(dX,dY,dZ,refENU.F,refENU.L)

        if self.initype == 'FLH' and hasattr(self,'sF'):
            if not np.isnan(self.sF):
                self.sE,self.sN,self.sU = geok.sFLH2sENU(self.F,self.L,self.H,self.sF,self.sL,self.sH)
        elif self.initype == 'XYZ' and hasattr(self,'sX'):
            if not np.isnan(self.sX):
                self.sE,self.sN,self.sU = geok.sXYZ2sENU(self.X,self.Y,self.Z,self.sX,self.sY,self.sZ)

    def keysanex(self):
        return list(self.anex.keys())


    def helmert_trans(self,params='itrf2008_2_etrf2000',invert=False):
        Xb = geok.helmert_trans(np.array([self.X,self.Y,self.Z]),params,invert)
        self.XYZset(*Xb)
        return None

    def velocity_trans(self, vx, vy, vz, epoc_init = 'auto' , epoc_end = 'auto'):
        """
        auto == epoc of the measures
        """

        if epoc_init == 'auto' and epoc_end == 'auto':
            print('ERR : epoc_init == auto and epoc_end == auto')
            return None

        tdt = geok.posix2dt(self.T)
        yeardec = geok.dt2year_decimal(tdt)

        if epoc_init  == 'auto':
            epoc_init = yeardec

        if epoc_end  == "auto":
            epoc_end = yeardec

        Xb = geok.itrf_speed_calc(self.X,self.Y,self.Z, epoc_init ,
                                  vx, vy, vz, epoc_end)
        self.XYZset(*Xb)
        return None



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
            self.T = geok.dt2posix(T)

        else:
            self.T = float(T)
            self.Tdt = geok.posix2dt(float(T))

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
        self.Q = geok.quaternion(self.R , self.P , self.Y , 'deg')
        return None


#class Evenement:
#    ''' Nature : '''
#
#    def __init__(self,T,Nature='event_NULL'):
#
#        if type(T) == dt.datetime:
#            self.Tdt = T
#            self.T = geok.dt2posix(T)
#
#        else:
#            self.T = float(T)
#            self.Tdt = geok.posix2dt(float(T))
#
#        allowed_nature = ['event_NULL','chgt_rec','chgt_ant']
#
#        self.Nature = Nature
#
#        if Nature in allowed_nature:
#            self.good_nature = True
#        else:
#            self.good_nature = False

#
#    def __call__(self):
#        return self.T, self.Nature
#
#    def __repr__(self):
#        return "{},{}".format(self.Tdt,self.Nature)
#
#

class TimeSeriePoint:
    def __init__(self,stat='STAT'):

        self.pts = []
        self.i_nomi = 0

        self.meta_set(stat=stat)

        self.boolENU = False
        self.bool_interp_uptodate = False
        self.bool_discont = False
        self.bool_discont_manu = False

        self.anex = dict()

    def __repr__(self):
        if self.pts == []:
            raise Exception('TS is empty ...')

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
        return len(self.pts)

    def meta_set(self,path='',stat='STAT',name=''):
        """
        stat:
            the 4char. code of the station
        name :
            is a free name for the TS, like the experience, the periode , the software ...
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
        self.pts = []

        self.bool_interp_uptodate = False
        self.bool_discont = False


    def readfile(self,filein):
        self.__dict__ = read_all_points(filein).__dict__

        self.interp_set()

    def add_point(self,inPoint):
        self.pts.append(inPoint)
        # this line is discontiued, because now nbpts is a property
        #self.nbpts = len(self.pts)

        self.bool_interp_uptodate = False

    def aleapt(self):
        ipt = randrange(self.nbpts)
        print("point no " + str(ipt))

        return self.pts[ipt]

    def startdate(self):
        self.sort()
        return geok.posix2dt(self.pts[0].T)

    def enddate(self):
        self.sort()
        return geok.posix2dt(self.pts[-1].T)

    def interval_nominal(self):

        if len(self.pts) < 2:
            self.i_nomi = 0
        else:
            Ttemp = np.sort([pt.T for pt in self.pts])
            self.i_nomi = np.round(np.min(np.diff(Ttemp)),1)

        return self.i_nomi


    def from_list(self,T,A,B,C,coortype='XYZ',sA=[],sB=[],sC=[]):
        """
        T in datetime
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

        self.sort()

        return None


    def to_list(self,coortype = 'XYZ',specific_output=''):

        if coortype == 'XYZ':
            A,B,C = 'X','Y','Z'
            sA,sB,sC = 'sX','sY','sZ'

        elif coortype == 'FLH':
            A,B,C = 'F','L','H'
            sA,sB,sC = 'sF','sL','sH'

        elif coortype == 'ENU':

            if self.boolENU == False:
                print("WARN : to_list : pas de coord. ENU pour " + self.name)
                return 0

            A,B,C = 'E','N','U'
            sA,sB,sC = 'sE','sN','sU'

        else:
            print("ERR : to_list : le coortype n'existe pas")

        if self.nbpts == 0:
            print("ERR : to_list : " + self.name + " the timeserie is empty")

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
        sq = np.squeeze
        outtup = (sq(A) , sq(B) , sq(C) , sq(T) , sq(sA) , sq(sB) , sq(sC))
        if specific_output == '':
            return outtup
        elif type(specific_output) is int:
            return outtup[specific_output]
        else:
            print("INFO : this mode must be implemented ;)")
            print("use an int as index instead")
            return outtup



    def sort(self):
        self.pts.sort(key=lambda x: x.T)

    def plot(self,coortype='ENU',diapt=2,alpha=0.8,fig=1,
             errbar=True,new_style=True,symbol = '.',errbar_width=1):
        """ fig can accept a int (id of a Figure)
         OR the figure Object itself """

        if new_style:
            styleint = 310
        else:
            styleint = 220

        try:
            A, B, C , T , sA, sB, sC = self.to_list(coortype=coortype)
        except TypeError as tyer:
            print("ERR   : unable to get coordinates")
            print("TRICK : check if the given coortype is in the timeserie")
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

        else:
            Atitle = 'A'
            Btitle = 'B'
            Ctitle = 'C'
            yylabel = 'displacement (??)'
            ABtitle = 'A & B'

        print("plot :", self.nbpts , "pts, ", self.stat)

        namest=0
        namend=10

        Tdt = geok.posix2dt(T)

        if type(fig) is int:
            figobj = plt.figure(fig)
        elif type(fig) is Figure:
            figobj = fig

        figobj.suptitle(self.stat)


        if self.name:
            name4plot = self.name[namest:namend]
        else:
            name4plot = self.stat


        plt.subplot(styleint+1)
        if errbar:
            plt.errorbar(Tdt,A,sA, fmt=symbol,label=name4plot,
                         markersize=diapt, alpha=alpha,ecolor='tab:cyan',
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

        ax = plt.gca()

#        if coortype == 'ENU':
##            refstr = 'ref XYZ = ' + genefun.join_improved(',',self.refENU.X ,
##                                                          self.refENU.Y ,
##                                                          self.refENU.Z)
        plt.subplot(styleint+2)
        if errbar:
            plt.errorbar(Tdt,B,sB, fmt=symbol,label=name4plot,
                         markersize=diapt, alpha=alpha,ecolor='tab:cyan',
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
                         markersize=diapt, alpha=alpha,ecolor='tab:cyan',
                         elinewidth=errbar_width)
        else:
            plt.plot(Tdt,C,symbol,label=name4plot,
                     markersize=diapt, alpha=alpha)
        try:
            plt.legend()
        except:
            pass
        plt.xlabel('Date')
        plt.ylabel(yylabel)
        plt.title(Ctitle)
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

        if not self.bool_discont:
            print("WARN : pas de discontinuité dans la TS")
            return 0

        if type(fig) is int:
            figobj = plt.figure(fig)
        elif type(fig) is Figure:
            figobj = fig

        for ax in figobj.axes:
            geok.plot_vertical_bar_ax(self.discont,ax,"r")

        if self.bool_discont_manu:
            for ax in figobj.axes:
                geok.plot_vertical_bar_ax(self.discont_manu,ax,"g")

#        figobj.axes[1]
#        geok.plot_vertical_bar(self.discont)
#
#        figobj.axes[2]
#        geok.plot_vertical_bar(self.discont)


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
        elif type(fig) is Figure:
            figobj = fig

        if not self.bool_discont:
            self.discont = []

        if not self.bool_discont_manu:
            self.discont_manu = []

        print("INFO : press SPACE to record a manual discontinuity")


        def onclick_discont(event):
            ix, iy = matplotlib.dates.num2date(event.xdata).replace(tzinfo=None), event.ydata

            print("INFO : discontinuity recorded : " , ix)
            for ax in figobj.axes:
                geok.plot_vertical_bar_ax([ix],ax,"g")
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

    def ENUcalc_from_mean_posi(self):
        self.ENUcalc(self.mean_posi())
        return None

    def from_uniq_point(self,Point,startdate,enddate,pas=1):
        self.del_data()

        if type(startdate) == dt.datetime:
            startdate = geok.dt2posix(startdate)
        if type(enddate) == dt.datetime:
            enddate = geok.dt2posix(enddate)

        N = int(np.round( enddate - startdate / pas ))
        #datelist = np.arange(startdate,enddate,pas)

        Point.Tset(startdate)

        for i in range(N):
            Point.Tset(Point.T + pas)
            self.add_point(copy.copy(Point))

        self.i_nomi = pas

    def timewin(self,windows,mode='keep'):
        '''IL EST TRES DANGEREUX DE L'APPLIQUER UN FENETRAGE A SOI MEME'''
        self.__dict__ = time_win(self,windows,mode).__dict__


    def interp_set(self,interptype = 'slinear'):
        if (not hasattr(self.pts[0],'E')) or np.isnan(self.pts[0].E) == True:
            print("WARN : interp_set : pas de ENU pour " + self.name)
        else:
            E,N,U,T,_,_,_ = self.to_list('ENU')

            self.EfT = interp1d(T,E,bounds_error=False,kind=interptype)
            self.NfT = interp1d(T,N,bounds_error=False,kind=interptype)
            self.UfT = interp1d(T,U,bounds_error=False,kind=interptype)

        if (not hasattr(self.pts[0],'X')) or np.isnan(self.pts[0].X) == True:
            print("WARN : interp_set : pas de XYZ pour " + self.name)
        else:

            X,Y,Z,T,_,_,_ = self.to_list('XYZ')

            self.XfT = interp1d(T,X,bounds_error=False,kind=interptype)
            self.YfT = interp1d(T,Y,bounds_error=False,kind=interptype)
            self.ZfT = interp1d(T,Z,bounds_error=False,kind=interptype)


        if (not hasattr(self.pts[0],'F')) or np.isnan(self.pts[0].L) == True:
            print("WARN : interp_set : pas de FLH pour " + self.name)
        else:
            F,L,H,T,_,_,_ = self.to_list('FLH')

            self.FfT = interp1d(T,F,bounds_error=False,kind=interptype)
            self.LfT = interp1d(T,L,bounds_error=False,kind=interptype)
            self.HfT = interp1d(T,H,bounds_error=False,kind=interptype)

        self.bool_interp_uptodate = True

    def interp_get(self,T,coortype='ENU'):

        if self.bool_interp_uptodate == False:
            print("WARN : interp obsolete, recalcul auto")
            self.interp_set()

        tsout = copy.copy(self)
        tsout.del_data()

        if not genefun.is_iterable(T):
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

        for i in range(len(T)):
            tsout.add_point(Point(A[i],B[i],C[i],T=T[i],initype=coortype))

        return tsout

    def set_discont(self,indiscont):
        self.discont = indiscont
        self.bool_discont = True

    def mean_posi(self,coortype='XYZ',outtype='point',meanormed='mean'):

        #special case where only one point
        if self.nbpts == 1:
            return copy.copy(self)

        A,B,C,T,sA,sB,sC = self.to_list(coortype=coortype)


        if meanormed == 'mean':
            Aout = np.nanmean(A)
            Bout = np.nanmean(B)
            Cout = np.nanmean(C)
        elif meanormed == 'med':
            Aout = np.median(A)
            Bout = np.median(B)
            Cout = np.median(C)

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
        print("NOTE 160415 : add_offset as method are hazardous ...")
        for pt in self.pts:
            pt.add_offset(dA,dB,dC)

    def decimate(self,dec):
        self.__dict__ = decimate_cleaner(self,dec).__dict__
        #decimate_cleaner(self,dec,True)

    def find_point(self,tin,tol=0.001,stop_when_found=True):

        find = False

        if type(tin) is dt.datetime:
            tin = geok.dt2posix(tin)

        tmin = self.pts[0].T
        tmax  = self.pts[-1].T

        if not ( tmin <= tin <= tmax ):
            print("WARN : find_point : t !€ [tmin , tmax]")

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
        print("selection device " , indtab)
        temp = read_all_obs(filein)[indtab]
        self.__dict__ = temp.__dict__

        self.interp_set()

    def add_obs(self,inObs):
        self.obs.append(inObs)
        self.nbobs = len(self.obs)

        self.bool_interp_uptodate = False

    def aleaobs(self):
        iobs = randrange(self.nbobs)
        print("observation no " + str(iobs))

        print(self.obs[iobs])

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
        self.__dict__ = time_win(self,windows,mode).__dict__

    def startdate(self):
        return self.obs[0].T

    def enddate(self):
        return self.obs[-1].T

    def to_list(self):
        if self.typeobs == 'NULL':
            print("ERR : pas de typeobs defini (NULL)")
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
            print("ERR : pas de typeobs defini (NULL)")
            return 0

        if self.typeobs == 'RPY':

            R,P,Y,T,_,_,_ = self.to_list()

            self.RfT = interp1d(T,R,bounds_error=False,kind=interptype)
            self.PfT = interp1d(T,P,bounds_error=False,kind=interptype)
            self.YfT = geok.interp1d_ang(T,Y,bounds_error=False,kind=interptype)

        self.bool_interp_uptodate = True

    def interp_get(self,T):

        if self.bool_interp_uptodate == False:
            print("WARN : interp obsolete, recalcul auto")
            self.interp_set()

        tsout = copy.copy(self)
        tsout.del_data()

        if not genefun.is_iterable(T):
            T = np.array([T])


        if self.typeobs == 'NULL':
            print("ERR : pas de typeobs defini (NULL)")
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

        print("plot :", self.nbobs , "obs., ", self.name)


        if self.typeobs == 'RPY':
            listtitle = ['','Roll','Pitch','Yaw']
        else:
            listtitle = ['','','','']

        namest=0
        namend=10
        Tdt = geok.posix2dt(T)

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

# ================================== FONCTION DE LECTURE ==================================

def read_all_points(filein):
    """selectionne automatiquement le type de fichier brut en entrée
       INPUT  : chemin du fichier brut de POINTS
       OUTPUT : Une TimeSeriePoint """

    firstline = open(filein).readline()

    if re.compile('RTKLIB').search(firstline):
        tsout = read_rtklib(filein)

    elif re.compile('Kinematic Processing').search(firstline):
        tsout = read_gipsy_apps(filein)

    elif re.compile('tdp.llh').search(filein):
        tsout = read_gipsy_bosser(filein)

    elif re.compile('STA').search(firstline) or re.compile('tdp').search(filein) or re.compile('TRPAZ').search(firstline):
        tsout = read_tdp(filein)

    elif re.compile('YY  MM DD HR MIN').search(firstline):
        tsout = read_track(filein)

    elif re.compile('Latitude').search(firstline):
        tsout = read_sonardyne_posi(filein)

    elif re.compile('Heading').search(firstline):
        tsout = read_sonardyne_attitude(filein)

    elif re.compile('\*\*\* warning').search(firstline):
        tsout = read_gins(filein,'kine')

    elif re.compile('#GINS_VERSION').search(firstline):
        tsout = read_gins_solution(filein)

    elif re.compile('OCTANS_ATTITUDE').search(firstline):
        tsout = read_qinsy(filein,2014,0o4,0o4)

    elif re.compile('PBO Station Position Time Series').search(firstline):
        tsout= read_pbo_pos(filein)

    elif re.compile('latitude_degre_decimal').search(firstline):
        tsout= read_nrcan_csv(filein)

    elif re.compile('--------------------------------------------------').search(firstline): # NDLR : Best parser ever
        tsout= read_nrcan_pos(filein)

    else:
        print("ERR : read_all : pas de motif valide pour lect. auto")
        print(filein)
        print(firstline)

    return tsout

def read_all_obs(filein):

    """selectionne automatiquement le type de fichier brut en entrée
       INPUT  : chemin du fichier brut de observations génériques
       OUTPUT : Une LISTE de TimeSerieObs """

    firstline = open(filein).readline()

    if re.compile('Heading').search(firstline):
        tsout = read_sonardyne_attitude(filein)

    else:
        print("ERR : read_all : pas de motif valide pour lect. auto")

    return tsout

def read_rtklib(filein):

    """lit un fichier de type RTKLIB
       INPUT  : chemin du fichier brut
       OUTPUT : Une TimeSeriePoint """

    tsout = TimeSeriePoint()

    for line in open(filein):

        if re.compile('e-baseline').search(line):
            initype='ENU'
        elif re.compile('x-ecef').search(line):
            initype='XYZ'
        elif 'latitude(deg)' in line:
            initype='FLH'

        if re.compile('%  UTC').search(line):
            dUTCGPS = 16
            print('WARN : MAYBE A WRONG LEAP SECOND !!!!')
        elif re.compile('%  GPST').search(line):
            dUTCGPS = 0

        if line[0] == '%':
            continue

        fields = line.split()
        date1 = re.findall(r"[\w']+",fields[0] + ':' + fields[1])
        date2 = tuple([int(d) for d in date1])

        T = (dt.datetime(date2[0],date2[1],date2[2],date2[3],date2[4],date2[5],date2[6]*1000) + dt.timedelta(seconds=dUTCGPS))
        A = (float(fields[2]))
        B = (float(fields[3]))
        C = (float(fields[4]))
        if initype == 'XYZ':
            sA = (float(fields[7]))
            sB = (float(fields[8]))
            sC = (float(fields[9]))
        elif initype == 'FLH':
            sA = (float(fields[7]))
            sB = (float(fields[8]))
            sC = (float(fields[9]))
            sA,sB,sC = geok.sENU2sFLH(A,B,C,sA,sB,sC)


        point = Point(A,B,C,T,initype,sA,sB,sC)
        point.anex['sdAB'] = float(fields[10])
        point.anex['sdBC'] = float(fields[11])
        point.anex['sdAC'] = float(fields[12])

        tsout.add_point(point)

    tsout.meta_set(filein)

    if initype == 'ENU':
        tsout.boolENU = True


    try:
        inpfilis = genefun.grep(filein,'inp file')
        tsout.anex['rover'] = os.path.basename(inpfilis[0].split()[-1])[0:4].upper()
        tsout.anex['base']  = os.path.basename(inpfilis[1].split()[-1])[0:4].upper()
    except:
        pass

    return tsout

def read_tdp(filein):

    print('TDPclassic')

    X,Y,Z = 0,0,0
    Tx , Ty , Tz, T = 111,222,333,0
    sX,sY,sZ = 0,0,0

    tsout = TimeSeriePoint()

    for line in open(filein):

        fields = line.split()

        if fields[4] == 'STA' and fields[5] == 'Z':
            Tz = geok.tgipsy2dt(fields[0])
            Z = (float(fields[2])* 1000)
            sZ = (float(fields[3])* 1000)

        if fields[4] == 'STA' and fields[5] == 'Y':
            Ty = geok.tgipsy2dt(fields[0])
            Y = (float(fields[2])* 1000)
            sY = (float(fields[3])* 1000)

        if fields[4] == 'STA' and fields[5] == 'X':
            Tx = geok.tgipsy2dt(fields[0])
            X = (float(fields[2])* 1000)
            sX = (float(fields[3])* 1000)
            STAT = fields[6]

        if  Tx == Ty == Tz :
            T = Tx
            point = Point(X,Y,Z,T,'XYZ',sX,sY,sZ)
            tsout.add_point(point)

            Tx = 111
            Ty = 222
            Tz = 333

    tsout.meta_set(filein,stat=STAT)

    return tsout


def read_gipsy_tdp(filein):
    """
    wrapper de read_tdp pour un nom plus explicite (et qui evitera de recoder la fct ...)
    """
    return read_tdp(filein)


def read_gipsy_tdp_list(filelistin):
    tslist = []
    for fil in filelistin:
        ts = read_gipsy_tdp(fil)
        tslist.append(ts)

    tsout = merge_ts(tslist)

    return tsout

def read_gipsy_bosser(filein):
    print('TDPBOSSER')
    F,L,H = 0,0,0
    T = 0
    sF,sL,sH = 0,0,0

    tsout = TimeSeriePoint()

    for line in open(filein):

        f = line.split()

        y = int(f[0])
        doydec = float(f[1])
        doy = int(doydec)

        T = geok.doy2dt(y,doy) + dt.timedelta(days=(doydec - doy))

        F = np.rad2deg(float(f[2]))
        L = np.rad2deg(float(f[3]))
        H = float(f[4])
        RMS = float(f[8])

        point = Point(F,L,H,T,'FLH')
        point.anex['RMS'] = RMS
        tsout.add_point(point)


    tsout.meta_set(filein,stat='STAT')

    return tsout

def read_gipsy_apps(filein):
    """ WARN : not optimized for > 1Hz !!! """
    print("GIPSY APPS")

    tsout = TimeSeriePoint()

    for l in open(filein):
        if l[0] == '#' or ('Kinematic Processing' in l):
            continue
        f = l.split()

        date_lis = [int(float(e)) for e in f[1].split(':')]

        T = dt.datetime(*date_lis)

        X  = float(f[2])
        sX = float(f[3])
        Y  = float(f[4])
        sY = float(f[5])
        Z  = float(f[6])
        sZ = float(f[7])

        point = Point(X,Y,Z,T,'XYZ',sX,sY,sZ)
        tsout.add_point(point)

    tsout.meta_set(filein)
    return tsout

#def read_gipsy_apps_tdp(filein):
#    """
#    double emploi avec read_tdp ... bravo champion
#    """
#
#    tsout = TimeSeriePoint()
#
#    TX = dt.datetime(1980,1,1)
#    TY = dt.datetime(1980,1,2)
#    TZ = dt.datetime(1980,1,3)
#
#    X , Y , Z = 0,0,0
#
#    for line in open(filein):
#        f = line.split()
#        stat = f[-1][-4:]
#        if not f[-3] == 'STA' and f[-4] in ('X','Y','Z'):
#            continue
#        elif f[-3] == 'STA' and f[-2] == 'X':
#            TX = geok.tgipsy2dt(float(f[0]))
#            X  = float(f[2])
#            sX = float(f[3])
#        elif f[-3] == 'STA' and f[-2] == 'Y':
#            TY = geok.tgipsy2dt(float(f[0]))
#            Y  = float(f[2])
#            sY = float(f[3])
#        elif f[-3] == 'STA' and f[-2] == 'Z':
#            TZ = geok.tgipsy2dt(float(f[0]))
#            Z  = float(f[2])
#            sZ = float(f[3])
#
#        if TX == TY == TZ:
#            point = Point(X,Y,Z,TX,'XYZ',sX,sY,sZ,name=stat)
#            tsout.add_point(point)
#
#    tsout.meta_set(filein,stat=stat)
#    return tsout


def read_track(filein):

    tsout = TimeSeriePoint()

    for line in open(filein):

        if re.compile('dNorth').search(line):
            initype='ENU'
        elif re.compile('dX').search(line):
            initype='XYZ'

        if line[0] != ' ':
            continue

        fields = line.split()

        T = (dt.datetime(int(fields[0]),
                         int(fields[1]),
                         int(fields[2]),
                         int(fields[3]),
                         int(fields[4]),
                         int(fields[5].split('.')[0]),
                         int(np.round(float(fields[5].split('.')[1]),4))))
        A = (float(fields[6]))
        B = (float(fields[8]))
        C = (float(fields[10]))

        sA = (float(fields[7]))
        sB = (float(fields[9]))
        sC = (float(fields[11]))

        # On inverse A et B car NEU => ENU
        if initype == 'ENU':
            point = Point(B,A,C,T,initype,sB,sA,sC)
        elif initype == 'XYZ':
            point = Point(A,B,C,T,initype,sA,sB,sC)
        else:
            "ERR : track_read : mauvais initype"

        tsout.add_point(point)

    if initype == 'ENU':
        tsout.boolENU = True

    tsout.meta_set(filein)

    tsout.anex['rover'] = tsout.name.split('.')[-2].upper()

    # recherche de la base
    try:
        if '.LC' in filein:
            stat_n_base_lc_files = glob.glob(filein[:-8] + '*')
            stat_n_base_lc_files.remove(filein)
            statbase = os.path.basename(stat_n_base_lc_files[0])[-7:-3]
            tsout.anex['base'] = statbase
        elif '.L1+L2' in filein:
            stat_n_base_l1l2_files = glob.glob(filein[:-11] + '*')
            stat_n_base_l1l2_files.remove(filein)
            statbase = os.path.basename(stat_n_base_l1l2_files[0])[-10:-6]
            tsout.anex['base'] = statbase
    except:
        print('WARN : unable to find the base for the TRACK experience')
        print(filein)
        pass

    return tsout

def read_gins_solution(filein,mode="cinematic"):
    """
    mode = cinematic : retrun a TimeSerie
    mode = static : retrun a point
    """

    F = open(filein)

    Pts_list_tmp = []

    for l in F:
        f = l.split()

        if 'STATION_NAME' in l:
            if len(f) > 1:
                namestat = f[1]
            else:
                namestat = f[-1]


        if l[0] == '#':
            continue

        Traw = float(f[2])

        if 'XYZ_SOL' in l:
            coordstype = 'XYZ'
            X =  float(f[4])
            Y =  float(f[6])
            Z =  float(f[8])

            sX =  np.sqrt(float(f[9]))
            sY =  np.sqrt(float(f[10]))
            sZ =  np.sqrt(float(f[11]))

            if "T24:00:00.000" in f[1]: # Manage the special case if we are at the border bw 2 days
                Txyz = geok.string_date2dt(f[1][:10])  + dt.timedelta(days=+1) + dt.timedelta(seconds=-19)
            else:
                Txyz = geok.string_date2dt(f[1]) + dt.timedelta(seconds=-19)

            point = Point(X,Y,Z,Txyz,coordstype,sX,sY,sZ,name=namestat)

            point.anex['sdXY'] = float(f[10])
            point.anex['sdXZ'] = float(f[11])
            point.anex['sdYZ'] = float(f[12])

        elif 'FLH_SOL' in l:
            coordstype = 'FLH'
            F =  float(f[4])
            L =  float(f[6])
            H =  float(f[8])

            sF =  np.rad2deg(np.sqrt(float(f[9])))
            sL =  np.rad2deg(np.sqrt(float(f[10])))
            sH =  np.sqrt(float(f[11]))

            point.F , point.L , point.H = F,L,H
            point.sF , point.sL , point.sH = sF,sL,sH

            #point.X , point.Y , point.Z = X,Y,Z

            point.anex['sdFL'] = float(f[10])
            point.anex['sdFH'] = float(f[11])
            point.anex['sdLH'] = float(f[12])

            Pts_list_tmp.append(point)

    #### End of reading, export
    if not Pts_list_tmp:
        print("WARN : no point found in :")
        print(filein)
        print("returns None")
        return None

    if mode == "cinematic":
        tsout = TimeSeriePoint()
        for point in Pts_list_tmp:
            tsout.add_point(point)
        tsout.meta_set(filein,namestat)
        return tsout

    elif mode == "static":
        pt_out = Pts_list_tmp[0]
        return pt_out


def read_gins_solution_multi(filein_list,return_dict = True):

    filein_list  = sorted(filein_list)
    Points_list  = []
    statname_stk = []

    for fil in filein_list:
        point_daily   = read_gins_solution(fil,mode="static")
        if not point_daily:
            continue
        Points_list.append(point_daily)
        statname_stk.append(point_daily.name)

    statname_uniq = sorted(list(set(statname_stk)))

    ts_dict = dict()

    for point in Points_list:
        if not point.name in ts_dict.keys():
            ts_dict[point.name] = TimeSeriePoint(stat=point.name)
        ts_dict[point.name].add_point(point)

    if return_dict:
        return ts_dict
    else:
        ts_list = []
        for k , val in ts_dict.items():
            ts_list.append(val)
        return ts_list


def read_epos_sta_kinematics(filein):
    """
    read an EPOS kinematic solutions
    """

    F = open(filein)
    Lines_4_DF_stk = []
    for l in F:
        fields = l.split()
        if l[0] != "K" and l[0] != "U" and l[0] != "X":
            continue
        if l[0] == "K" or l[0] == "U" or l[0] == "X":
            namstat = fields[2]
            numstat = int(fields[1])
            MJD_epo = float(fields[3])
            numobs = int(fields[4])

            X = float(fields[6])
            Y = float(fields[7])
            Z = float(fields[8])
            sX = float(fields[10])
            sY = float(fields[11])
            sZ = float(fields[12])

            N = float(fields[14])
            E = float(fields[15])
            U = float(fields[16])
            sN = float(fields[18])
            sE = float(fields[19])
            sU = float(fields[20])

            tup_4_df = (namstat,numstat,MJD_epo,numobs,X,Y,Z,sX,sY,sZ,
                        N,E,U,sN,sE,sU)
            Lines_4_DF_stk.append(tup_4_df)

    columns = ("site","site_num",
                   "MJD_epo","numobs",
                   "x","y","z","sx","sy","sz",
                   "N","E","U","sN","sE","sU")

    DFout = pd.DataFrame(Lines_4_DF_stk,
                     columns=columns)
    return DFout

def read_epos_sta_coords_mono(filein,return_df=True):
    """
    read an EPOS's YYYY_DDD_XX_sta_coordinates coordinates files
    and return a list of Points objects

    ... TBC ...
    """
    F = open(filein)

    Points_list_stk = []
    Lines_4_DF_stk = []

    for l in F:
        fields = l.split()
        if l[0] != " ":
            continue
        if "SITE" in fields[0]:
            namestat = fields[8]
            numstat  = int(fields[2])
            tecto_plate = fields[4]
            MJD_ref  = int(fields[5])
            MJD_strt = int(fields[6])
            MJD_end  = int(fields[7])
            MJD_mid = np.mean([MJD_strt , MJD_end])
            T = geok.MJD2dt(MJD_mid)

        if "POS_VEL:XYZ" in fields[0]:
            X  = float(fields[4])
            Y  = float(fields[5])
            Z  = float(fields[6])
            Vx = float(fields[7])
            Vy = float(fields[8])
            Vz = float(fields[9])

        if "SIG_PV_XYZ" in fields[0]:
            sX  = float(fields[4].replace("D","E"))
            sY  = float(fields[5].replace("D","E"))
            sZ  = float(fields[6].replace("D","E"))
            sVx = float(fields[7])
            sVy = float(fields[8])
            sVz = float(fields[9])

            #### Last useful line for the point, store it
            point = Point(X,Y,Z,T,"XYZ",sX,sY,sZ,name=namestat)
            point.anex["Vx"] = sVx
            point.anex["Vy"] = sVy
            point.anex["Vz"] = sVz
            Points_list_stk.append(point)

            #### And store for the DataFrame
            tup_4_DF = (namestat,numstat,tecto_plate,
                        MJD_ref,MJD_strt,MJD_end,
                        X,Y,Z,sX,sY,sZ,
                        Vx,Vy,Vz,sVx,sVy,sVz)



            Lines_4_DF_stk.append(tup_4_DF)


        columns = ("site","site_num","tecto_plate",
                   "MJD_ref","MJD_start","MJD_end",
                   "x","y","z","sx","sy","sz",
                   "Vx","Vy","Vz","sVx","sVy","sVz")

        DFout = pd.DataFrame(Lines_4_DF_stk,
                     columns=columns)


    if return_df:
        return DFout
    else:
        return Points_list_stk

def read_epos_sta_coords_multi(filein_list,return_dict = True):

    filein_list  = sorted(filein_list)
    Points_list  = []
    statname_stk = []

    for fil in filein_list:
        Points_daily_list = read_epos_sta_coords_mono(fil,return_df=False)
        Points_list   = Points_list + Points_daily_list
        statname_stk  = statname_stk + [e.name for e in Points_daily_list]

    statname_uniq = sorted(list(set(statname_stk)))

    ts_dict = dict()

    for point in Points_list:
        if not point.name in ts_dict.keys():
            ts_dict[point.name] = TimeSeriePoint(stat=point.name)
        ts_dict[point.name].add_point(point)

    if return_dict:
        return ts_dict
    else:
        ts_list = []
        for k , val in ts_dict.items():
            ts_list.append(val)
        return ts_list


def read_epos_slv_times(p):
    L = genefun.extract_text_between_elements_2(p,"\+sum_times/estimates","\-sum_times/estimates")


    Lgood = []
    for l in L[1:-1]:
        if "EPOCHE" in l:
            cur_epoc_line = l
            cur_epoc_f = cur_epoc_line.split()
            cur_epoc   = geok.MJD2dt(int(cur_epoc_f[1])) +  dt.timedelta(seconds=int(86400*float(cur_epoc_f[2])))

        if re.match("^   [0-9]{4}.*",l):
            Lgood.append([cur_epoc] + [float(e) for e in l.split()])


    DF = pd.DataFrame(Lgood,columns=["epoch","stat","offset","offset_sig"])

    DF["stat"] = DF["stat"].astype('int')

    return DF



def write_epos_sta_coords(DF_in,file_out):

    DF_work = DF_in.sort_values(["site","MJD_start"])


    Stat_lines_blk_stk = []

    generic_header = """+info
 FLATTENING                  298.2550
 MAJOR_AXIS              6378140.0000
 REFERENCE_FRAME                IGS14
 NUMBER_OF_STATIONS             {:5d}
 REF_MJD                        {:5d}
-info
"""

    generic_header = generic_header.format(len(DF_work),
                                           genefun.most_common(DF_work["MJD_ref"]))

    Stat_lines_blk_stk.append(generic_header)

    Stat_lines_blk_stk.append("+station_coordinates")

    for site in DF_work["site"].unique():

        Stat_lines_blk_stk.append("*------------------------- ---- ----- -beg- -end- -**- ------------------------------------------------\n*")

        for i_l ,(_ , l) in enumerate(DF_work[DF_work["site"] == site].iterrows()):

            line_site_fmt = " SITE            m {:4d}  {:1d} {:} {:5d} {:5d} {:5d} {:}   A  0      LOG_CAR       LOG_CAR"
            line_posi_fmt = " POS_VEL:XYZ     m {:4d}  {:1d} {:+15.4f} {:+15.4f} {:+15.4f}      {:+6.4f} {:+6.4f} {:+6.4f}"
            line_velo_fmt = " SIG_PV_XYZ      m {:4d}  {:1d} {:+15.4f} {:+15.4f} {:+15.4f}      {:+6.4f} {:+6.4f} {:+6.4f}"

            line_site = line_site_fmt.format(l["site_num"],i_l,l["tecto_plate"].upper(), l["MJD_ref"],l["MJD_start"],l["MJD_end"],l["site"])
            line_posi = line_posi_fmt.format(l["site_num"],i_l,l["x"],l["y"], l["z"],l["sx"],l["sy"],l["sz"])
            line_velo = line_velo_fmt.format(l["site_num"],i_l,l["Vx"],l["Vy"], l["Vz"],l["sVx"],l["sVy"],l["sVz"])

            Stat_lines_blk_stk.append(line_site)
            Stat_lines_blk_stk.append(line_posi)
            Stat_lines_blk_stk.append(line_velo)
            Stat_lines_blk_stk.append("*")


    Stat_lines_blk_stk.append("-station_coordinates")


    final_str = "\n".join(Stat_lines_blk_stk)


    with open(file_out,"w+") as f:
        f.write(final_str)

    return final_str


def prn_int_2_prn_str(prn_int,full_out=False):
    """
    for read_combi_sum_full

    if full_out : return e.g. "G04","G",4
    """
    const = "X"

    prn_int = int(prn_int)

    prn_int_out = prn_int

    if prn_int >= 400:
        prn_int_out = prn_int - 400
        const = "J"
    elif 300 <= prn_int < 400:
        prn_int_out = prn_int - 300
        const = "C"
    elif 200 <= prn_int < 300:
        prn_int_out = prn_int - 200
        const = "E"
    elif 100 <= prn_int < 200:
        prn_int_out = prn_int - 100
        const = "R"
    else:
        prn_int_out = prn_int
        const = "G"

    prn_str = const + str(prn_int_out).zfill(2)

    if not full_out:
        return prn_str
    else:
        return prn_str , const , prn_int_out

def read_combi_sum_full(sum_full_file,RMS_lines_output=True,
                        set_PRN_as_index=True):
    Vals_stk = []

    for l in open(sum_full_file):

        F = l.split()
        if "|" in F:
            F.remove("|")

        ### Find date line
        if "MJD:" in l:
            date_line  = l

        ### Skip useless lines
        if not "|" in l or "------" in l:
            continue

        ### Find AC list
        if "PRN" in l:
            ACs_list  = F
            ACs_list.append("RMS_sat")
            ACs_list.append("PRN_str")
            ACs_list.append("CONST")


        elif F[0].isnumeric():
            Fout = [float(f) for f in F]
            Fout[0] = int(Fout[0])
            #Add the PRN string and the constellation
            Fout.append(prn_int_2_prn_str(int(Fout[0])))
            Fout.append(Fout[-1][0])

            Vals_stk.append(Fout)

        elif "RMS" in F[0] and RMS_lines_output:
            Fout = [float(f) for f in F[1:]]
            Fout.append(np.nan)
            Fout.insert(0,F[0])
            #Add FAKE the PRN string and the constellation
            Fout.append(F[0])
            Fout.append(None)

            Vals_stk.append(Fout)


    DF = pd.DataFrame(Vals_stk,columns=ACs_list)

    ### Date management
    mjd = float(date_line.split("MJD:")[1].split()[0])
    date_dt = geok.MJD2dt(mjd)

    DF.date_mjd = mjd
    DF.date_dt  = date_dt
    DF.date_gps = genefun.join_improved("",geok.dt2gpstime(date_dt))

    if set_PRN_as_index:
        DF.set_index("PRN_str",inplace=True)

    return DF


def read_combi_sum_exclu(sum_file,return_as_df=True,
                         use_intuitive_bool = True):


    t_dt = geok.sp3name2dt(sum_file)

    with open(sum_file) as f:
        cont = f.readlines()


    excluded_dic = dict()
    useful_ssection = False
    useful_ssection_k = 0
    for l in cont:
        f = l.split()
        if "---|---" in l and useful_ssection_k < 2:
            useful_ssection = not useful_ssection
            useful_ssection_k +=1
            continue

        if not useful_ssection:
            continue

        prn_raw = f[0]

        if "X" in prn_raw:
            exclu = True
        else:
            exclu = False

        if use_intuitive_bool:
            exclu = not exclu

        prn_int = int(f[0].replace("X","").split()[0])

        prn_good = prn_int_2_prn_str(prn_int)

        excluded_dic[prn_good] = exclu

    if return_as_df:
        return pd.DataFrame(excluded_dic,index=[t_dt])
    else:
        return excluded_dic


def read_combi_clk_rms(sum_file,return_as_df=True,
                       clk_ref_cen_gal = "com"):
    """
    based on : read_good_clk_rms_one
    """

    strt = " RESULTS OF FINAL WEIGHTED COMBINATION"
    end  = " CLK_REF_CEN_GAL: " + clk_ref_cen_gal

    L = genefun.extract_text_between_elements_2(sum_file,strt,end)

    L = L[:-2]

    Lres = [e for e in L if re.search("^ [a-z]{3} \|",e)]

    Lres_splited = [e.split() for e in Lres]

    filnam = os.path.basename(sum_file)
    if "log" in filnam:
        week = int(filnam[4:8])
        dow  = int(filnam[9])
        tdt = geok.gpstime2dt(week,dow)
    elif "cls" in filnam:
        week = int(filnam[3:7])
        dow  = int(filnam[7])
        tdt = geok.gpstime2dt(week,dow)

    rms_dict = dict()

    for e in Lres_splited:
        try:
            rms_dict[e[0]] = int(e[-4])
        except:
            print("WARN : ", e[-4],"not handeled")
            print("replaced with NaN")
            rms_dict[e[0]] = np.nan

    if return_as_df:
        return pd.DataFrame(rms_dict,index=[tdt])
    else:
        return tdt,rms_dict


def read_combi_clk_rms_full_table(path_in):
    """
    recommended for .out file
    """
    strt = "RMS \(ps\) OF AC CLOCK COMPARED TO COMBINATION"
    end  = "---+---"

    Lines = genefun.extract_text_between_elements_2(path_in,strt,end,nth_occur_elt_end=1)

    Lines_good = []

    for l in Lines[1:]:
        if "---+---" in l or "bad" in l:
            continue
        else:
            Lines_good.append(l)


    Lines_good = [e.replace("|","") for e in Lines_good]
    Lines_good = [e.replace("         ","  SAT    ") for e in Lines_good]

    STR = "".join(Lines_good)

    import io

    DF = pd.read_table(io.StringIO(STR),
                      na_values = ["-","X",">>>"] ,
                      delim_whitespace = True ,
                      error_bad_lines=False)

    DF = DF.set_index("SAT")

    return DF


def read_combi_REPORT(Path_list):
    STK = []
    for p in Path_list:
        F = open(p)
        for l in F:
            f = l.split()
            if "epoch" in l:
                epoch = geok.gpstime2dt(int(f[2]),int(f[3]))
            if "orb_flag_x" in l:
                prn_str , const , prn_int = prn_int_2_prn_str((f[2]),True)
                STK.append((epoch,prn_str,const , prn_int,"all"))
            if "orb_excl_sat" in l:
                prn_str , const , prn_int = prn_int_2_prn_str((f[3]),True)
                STK.append((epoch,prn_str,const , prn_int,f[2]))

    DF = pd.DataFrame(STK,columns=("epoch","PRN_str","CONST","PRN","AC"))

    return DF




def read_gins(filein,kineorstatic='kine',flh_in_rad=True,
              force_get_convergence=False,kf_result=False):

    '''
    Static : donne un point
    Kinematic : donne une TS

    force_get_convergence : if there is a bug about
    'COORDONNEES DES STATIONS AJUSTEES EN HAUTE FREQUENCE' field, it is in the listing
    but empty ... so we force the retreive of the 'c o n v e r g e n c e' part
    '''

    if '.prepars' in filein:
        print('WARN :',filein, 'seems to be a prepars file, are you sure of what you are doing ?')

    # Pour le static, il y a des blancs en fin de ligne ...

    if kineorstatic == 'kine':
        #regex = '\[S[PLHXYZ] .*\]$'
        regex = '\[S[PLHXYZ] .*\]$'
        tsout = TimeSeriePoint()
        if kf_result:
            regex = '\[S[PLHXYZ][E ].*\]     $'
    elif kineorstatic == 'static':
        #regex = '\[S[PLHXYZ] .*\]     $'
        regex = '\[S[PLHXYZ][E ].*\]     $'
        tsout = TimeSeriePoint()
    else:
        print("ERR")

    A,B,C = 0,0,0
    Ta , Tb , Tc, T = 111,222,333,0
    sA,sB,sC = 0,0,0

    if kf_result:
        regex = '\[S[PLHXYZ][E ].*\]     $'


    # Specific si 2ble convergence
    grep_conv = genefun.grep(filein,'c o n v e r g e n c e')
    if len(grep_conv) == 2:
        IPPmode = True
        converg_compt = 0
        print("INFO : ", os.path.basename(filein) , 'have 2  c o n v e r g e n c e  fields')
        print("     keeping the last one")
    else:
        IPPmode = False

    # Specific si Ajustement Final
    greped_adj = genefun.grep(filein,'COORDONNEES DES STATIONS AJUSTEES EN HAUTE FREQUENCE')
    if force_get_convergence:
        FinalAdj_mode = False
    elif len(greped_adj) != 0:
        FinalAdj_mode = True
        FinalAdj_found = False
        print("INFO : ", os.path.basename(filein) , 'have a COORD STAT AJ EN HTE FREQ  field')
        print("     keeping this one (and not the convergence)")
    else:
        FinalAdj_mode = False

    fileopened = open(filein,encoding = "ISO-8859-1")

    regex_valid_line_count = 0

    for line in fileopened:

        if re.compile('__Nom__').search(line):
            nextline = next(fileopened).split()
            namestat = nextline[3]
            Xref = float(nextline[4])
            Yref = float(nextline[5])
            Zref = float(nextline[6])

            Fref , Lref , Href = geok.XYZ2GEO(Xref,Yref,Zref)

#        if 'angles en deg' in line:
#            flh_in_rad = False

        # Specific search
        if IPPmode and re.compile('c o n v e r g e n c e').search(line):
            converg_compt = converg_compt + 1
        if FinalAdj_mode and re.compile('COORDONNEES DES STATIONS AJUSTEES EN HAUTE FREQUENCE').search(line):
            FinalAdj_found = True

        # Specific skip
        if FinalAdj_mode and not FinalAdj_found:
            continue
        if IPPmode and converg_compt != 2:
            continue

        if "real" in line:
            rawexectime = line.split()[-1]

        if re.compile(regex).search(line):
            regex_valid_line_count += 1
            fields = line.split()

            if re.compile('[XYZ]').search(line):
                initype = 'XYZ'
                Aref , Bref ,Cref = Xref,Yref,Zref
            elif re.compile('[PLH]').search(line):
                initype = 'FLH'
                Aref , Bref ,Cref = Fref,Lref,Href
            else:
                print("ERR : read gins : wrong initype")

            if (float(fields[2]) == 0):
                continue

            if (fields[0] == 'stations'):
                print('yyyyy')

            # securité pour les lignes du type
            #  ------------------------------------------------------------------------------------------------------
            # stations      175  -0.574353783560234E+07  +/-   0.000000000000000E+00   1  [SX  1212001892701M005071]

            jour = int(line[108:110])
            h    = int(line[110:112])
            m    = int(line[112:114])
            s    = int(line[114:116])
            yy   = int(line[125:127])

            # gestion des années
            if 80 < yy <= 99:
                yy = yy + 1900
            else:
                yy = yy + 2000

            # pour le mois, si > sept (9), alors lettre ...
            mm = line[127]

            if mm == 'O':
                mm = 10
            elif mm == 'N':
                mm = 11
            elif mm == 'D':
                mm = 12
            else:
                mm = int(mm)

            try:
                if h == 24:
                    # cas exceptionnel ou on doit gerer minuit
                    # on retranche l'heure dans l'int en input et on l'ajoute dans le dt
                    Ttemp = (dt.datetime(yy,mm,jour,h-1,m,s) + dt.timedelta(seconds=-19) + dt.timedelta(hours=1))
                else:
                    Ttemp = (dt.datetime(yy,mm,jour,h,m,s) + dt.timedelta(seconds=-19))

                if  line[105] == 'X' or line[105] == 'P':
                    Ta = Ttemp
                    A = (float(fields[3]))
                    sA = (float(fields[4]))

                if  line[105] == 'Y' or line[105] == 'L':
                    Tb = Ttemp
                    B = (float(fields[3]))
                    sB = (float(fields[4]))

                if  line[105] == 'Z' or line[105] == 'H':
                    Tc = Ttemp
                    C = (float(fields[3]))
                    sC = (float(fields[4]))
            except ValueError as err:
                print('ERR :')
                print('yy,mm,jour,h,m,s ',yy,mm,jour,h,m,s)
                raise err

            if  Ta == Tb == Tc :
                T = Ta
                if initype == 'FLH' and not FinalAdj_mode:
                    if flh_in_rad:
                        A = np.rad2deg(A)
                        B = np.rad2deg(B)
                        sA = np.rad2deg(sA)
                        sB = np.rad2deg(sB)
                if kf_result and 0:
                    A = A + Aref
                    B = B + Bref
                    C = C + Cref
                point = Point(A,B,C,T,initype,sA,sB,sC,name=namestat)

                Ta = 111
                Tb = 222
                Tc = 333

                if kineorstatic == 'static':
                    return point
                elif kineorstatic == 'kine':
                    tsout.add_point(point)
                else:
                    print("ERR")

    if  regex_valid_line_count == 0:
        print("WARN : no valid line (with regex check) was found !!!")


    tsout.anex['exec_time'] = rawexectime
    tsout.meta_set(filein,namestat)
    return tsout


def gins_read_time(line):
    jour = int(line[108:110])
    h    = int(line[110:112])
    m    = int(line[112:114])
    s    = int(line[114:116])
    yy   = int(line[125:127])

    # gestion des années
    if 80 < yy <= 99:
        yy = yy + 1900
    else:
        yy = yy + 2000

    # pour le mois, si > sept (9), alors lettre ...
    mm = line[127]

    if mm == 'O':
        mm = 10
    elif mm == 'N':
        mm = 11
    elif mm == 'D':
        mm = 12
    else:
        mm = int(mm)

    try:
        if h == 24:
            # cas exceptionnel ou on doit gerer minuit
            # on retranche l'heure dans l'int en input et on l'ajoute dans le dt
            Ttemp = (dt.datetime(yy,mm,jour,h-1,m,s) +
                     dt.timedelta(seconds=-19) + dt.timedelta(hours=1))
        else:
            Ttemp = (dt.datetime(yy,mm,jour,h,m,s) + dt.timedelta(seconds=-19))
    except:
        Ttemp = dt.datetime(1970,1,1)

    T = Ttemp

    return T


def gins_read_MZB(filein):

    F = open(filein)

    regex = '\[MZB.*\]     $'

    Tstk    = []
    MZBstk  = []
    sMZBstk = []

    for line in F:
        if re.compile('__Nom__').search(line):
            nextline = next(F).split()
            namestat = nextline[3]
            Xref = float(nextline[4])
            Yref = float(nextline[5])
            Zref = float(nextline[6])

        if re.compile(regex).search(line):
            fields = line.split()

            #[MZB  801980015051301 GPS]

            Traw = fields[6][-8:]
            yy = int(Traw[0:2])
            mm = int(Traw[2:4])
            dd = int(Traw[4:6])
            hh = int(Traw[6:8])

            if 80 < yy <= 99:
                yy = yy + 1900
            else:
                yy = yy + 2000

            T = dt.datetime(yy,mm,dd,hh)

            MZB  = float(fields[3])
            sMZB = float(fields[4])

            Tstk.append(T)
            MZBstk.append(MZB)
            sMZBstk.append(sMZB)

    return Tstk , MZBstk , sMZBstk , namestat


def write_ATM_GAMIT(Tstk , MZBstk , sMZBstk ,
                    namestat , file_out):
    Fout = open(file_out,'w+')
    for T,mzb , smzb in zip(Tstk , MZBstk , sMZBstk):
        yy = T.year
        mm = T.month
        dd = T.day
        hh = T.hour
        Line = 'ATM_ZEN X {}  1 {:4} {:2} {:2} {:2}  0  {:6.4f} +-   {:6.4f}    {:6.4f}\n'.format(namestat.upper(),yy,mm,dd,hh,mzb,smzb,mzb)
        Fout.write(Line)
    Fout.close()
    return file_out


def MZB_GINS_2_ATM_GAMIT(listing_in,path_out):
    Tstk , MZBstk , sMZBstk , namestat = gins_read_MZB(listing_in)
    doy,yy = geok.dt2doy_year(Tstk[0],str)
    file_out = os.path.join(path_out,'_'.join(('MZB_GINS_2_ATM_GAMIT',namestat,doy,yy,'.txt')))
    write_ATM_GAMIT(Tstk , MZBstk , sMZBstk ,
                    namestat , file_out)

    return file_out

def read_gins_wrapper(input_list_or_path,flh_in_rad=True):
    """
    pour une liste de paths renvoie une liste de timeseries
    AUTANT DE TIMESERIES QUE DE LISTINGS

    INPUT :
    Une liste de fichiers
    Un path compatible avec glob

    les fichiers sans extensions .gins
    sont automatiquement exclus
    """

    if type(input_list_or_path) is str:
        gins_listings_list = glob.glob(input_list_or_path)
    else:
        gins_listings_list = input_list_or_path

    tslis = []
    for f in gins_listings_list:
        if not '.gins' in f:
            print('WARN : no .gins ext, skipping')
            continue

        ts = read_gins(f,flh_in_rad=flh_in_rad)
        tslis.append(ts)

    return tslis


def convert_sp3_clk_2_GINS_clk(sp3_path_in,
                               clk_gins_out,
                               interpo_30sec = True,
                               return_as_DF = True):
    DF = gcls.read_sp3(sp3_path_in)

    Fout = open(clk_gins_out,"w+")

    def write_GINS_signaletic_elt_clk(dt_in,sv,
                                      signaletic_name="MNG",
                                      add_19sec_to_dt_in = True):
        """
        very beta, only for GPS clk
        """

        if add_19sec_to_dt_in:
            dt_work = dt_in + dt.timedelta(seconds=19)
        else:
            dt_work = dt_in

        jjul = geok.dt2jjulCNES(dt_work)
        sec_in_day = (dt_work - geok.jjulCNES2dt(jjul) ).seconds

        #MNG0000000jjjjjcccccnnnn

        outstr = "[MNG0000000" + str(jjul) + str(sec_in_day).zfill(5) + "GP" + str(sv).zfill(2) + "]"

        return outstr

    c = 299792458

    if interpo_30sec:
        Epoc_work = []
        Sv_work   = []
        Clk_work  = []

        for sv in sorted(DF["sv"].unique()):
            DFsv = DF[DF["sv"] == sv]

            Epoc_inp = np.array(geok.dt2posix(DFsv["epoch"]))
            Clk_inp  = np.array(DFsv["clk"])

            Epoc_interp = np.arange(np.min(Epoc_inp),np.max(Epoc_inp),30)
            Epoc_interp_dt = geok.posix2dt(Epoc_interp)

            I = interpolate.interp1d(Epoc_inp,Clk_inp)

            Clk_interp = I(Epoc_interp)

            Sv_work   = Sv_work   + [sv] * len(Epoc_interp)
            Epoc_work = Epoc_work + list(Epoc_interp_dt)
            Clk_work  = Clk_work  + list(Clk_interp)

    else:
        Epoc_work = DF["epoch"]
        Sv_work   = DF["sv"]
        Clk_work  = DF["clk"]


    DF_work = pd.DataFrame(list(zip(Epoc_work,Sv_work,Clk_work)),
                            columns=("epoch","sv","clk"))

    DF_work.sort_values(["epoch","sv"],inplace=True)

    print("tot",DF_work)

    for epoc , sv , clk in zip(DF_work["epoch"],DF_work["sv"],DF_work["clk"]):
        signaletik = write_GINS_signaletic_elt_clk(epoc,sv,clk)

        str_final = " 0 0 {:}  {:+17.15e} {:+17.15e}\n".format(signaletik,clk* 10**-6 * c,0)

        Fout.write(str_final)

    if not return_as_DF:
        return clk_gins_out
    else:
        return DF_work



#def read_gins_kinematic(filein):
#
#    ''' retourne une TSPoint '''
#
#    tsout = TimeSeriePoint()
#
#    X,Y,Z = 0,0,0
#    Tx , Ty , Tz, T = 111,222,333,0
#    sX,sY,sZ = 0,0,0
#
#    for line in open(filein):
#
#        if re.compile('\[S[XYZ] .*\]$').search(line):
#            fields = line.split()
#
#            if (float(fields[2]) == 0):
#                continue
#
#            if (fields[0] == 'stations'):
#                continue
#            # securité pour les lignes du type
#            #  ------------------------------------------------------------------------------------------------------
#            # stations      175  -0.574353783560234E+07  +/-   0.000000000000000E+00   1  [SX  1212001892701M005071]
#
#            jour = int(line[108:110])
#            h = int(line[110:112])
#            m = int(line[112:114])
#            s = int(line[114:116])
#            yy = int(line[125:127]) + 2000
#            mm = int(line[127])
#
#            if line[105] == 'X':
#                Tx = (dt.datetime(yy,mm,jour,h,m,s) + dt.timedelta(seconds=-19))
#                X = (float(fields[3]))
#                sX = (float(fields[4]))
#
#            if  line[105] == 'Y':
#                Ty = (dt.datetime(yy,mm,jour,h,m,s) + dt.timedelta(seconds=-19))
#                Y = (float(fields[3]))
#                sY = (float(fields[4]))
#
#            if  line[105] == 'Z':
#                Tz = (dt.datetime(yy,mm,jour,h,m,s) + dt.timedelta(seconds=-19))
#                Z = (float(fields[3]))
#                sZ = (float(fields[4]))
#
#
#            if  Tx == Ty == Tz :
#                T = Tx
#                point = Point(X,Y,Z,T,'XYZ',sX,sY,sZ)
#                tsout.add_point(point)
#                Tx = 111
#                Ty = 222
#                Tz = 333
#
#    tsout.meta_set(filein)
#
#    return tsout
#
#
#def read_gins_static_solo(filein):
#
#    X,Y,Z = 0,0,0
#    Tx , Ty , Tz, T = 111,222,333,0
#    sX,sY,sZ = 0,0,0
#    namestat='NULL'
#
#    fileopened = open(filein)
#
#    for line in fileopened:
#
#        if re.compile('__Nom__').search(line):
#            namestat = next(fileopened).split()[3]
#
#        # Pour le static, il y a des blancs en fin de ligne ...
#        if re.compile('\[S[XYZ] .*\]     $').search(line):
#
#            fields = line.split()
#
#            if (float(fields[2]) == 0):
#                continue
#
#            if (fields[0] == 'stations'):
#                continue
#            # securité pour les lignes du type
#            #  ------------------------------------------------------------------------------------------------------
#            # stations      175  -0.574353783560234E+07  +/-   0.000000000000000E+00   1  [SX  1212001892701M005071]
#
#            jour = int(line[108:110])
#            h = int(line[110:112])
#            m = int(line[112:114])
#            s = int(line[114:116])
#            yy = int(line[125:127])
#
#            # gestion des annÃ©es
#            if 80 < yy <= 99:
#                yy = yy + 1900
#            else:
#                yy = yy + 2000
#
#            # pour le mois, si > sept (9), alors lettre ...
#            mm = line[127]
#
#            if mm == 'O':
#                mm = 10
#            elif mm == 'N':
#                mm = 11
#            elif mm == 'D':
#                mm = 12
#            else:
#                mm = int(mm)
#
#            if line[105] == 'X':
#                Tx = (dt.datetime(yy,mm,jour,h,m,s) + dt.timedelta(seconds=-19))
#                X = (float(fields[3]))
#                sX = (float(fields[4]))
#
#            if  line[105] == 'Y':
#                Ty = (dt.datetime(yy,mm,jour,h,m,s) + dt.timedelta(seconds=-19))
#                Y = (float(fields[3]))
#                sY = (float(fields[4]))
#
#            if  line[105] == 'Z':
#                Tz = (dt.datetime(yy,mm,jour,h,m,s) + dt.timedelta(seconds=-19))
#                Z = (float(fields[3]))
#                sZ = (float(fields[4]))
#
#
#            if  Tx == Ty == Tz :
#                T = Tx
#                point = Point(X,Y,Z,T,'XYZ',sX,sY,sZ,name=namestat)
#
#                Tx = 111
#                Ty = 222
#                Tz = 333
#
#                return point

#read_gins_static_solo('/media/DDannex/GINS_AKRIM/listing/listing_NRMD/DIR_sortie_2007_20830_nrmd.140513_060821.gins')

def read_gins_multi_raw_listings(filelistin,kineorstatic='static',flh_in_rad=True):
    """
    traite une liste de listing bruts
    pour obtenir UNE SEULE timeserie
    """

    tsout = TimeSeriePoint()
    refname = 'RIEN'
    if kineorstatic == 'static':
        for filein in filelistin:
            print(filein)
            if not genefun.check_regex(filein,'c o n v e r g e n c e'):
                continue
            pt = read_gins(filein,kineorstatic='static',flh_in_rad=flh_in_rad)
            if refname == 'RIEN':
                refname = pt.name
            if refname != pt.name:
                print("WARN : read_gins_multi : nom de stat. != reference")
            tsout.add_point(pt)
        tsout.meta_set(stat=refname)


    elif kineorstatic == 'kine':
        tsoutlis = []
        for filein in filelistin:
            print(filein)
            if not genefun.check_regex(filein,'c o n v e r g e n c e'):
                continue
            ts = read_gins(filein,kineorstatic='kine',flh_in_rad=flh_in_rad)
            tsoutlis.append(ts)
        tsout = merge_ts(tsoutlis)

    else:
        print("ERR : check kineorstatic keyword")

    tsout.sort()
    return tsout

def read_gins_multi_extracted(filelistin,flh_in_rad=True):
    """
    traite les extractions de listings
    sous la forme par ex. S<HLP>__ddhhiissxxxxxxxxxyym.HOUE
    La liste doit contenir exactement 3 fichiers (pour chacune des composantes)
    """
    tsout = TimeSeriePoint()
    if len(filelistin) != 3:
        print("ERR : read_gins_multi_extracted : listfilein != 3 elts")
        return None
    statnameset = list(set([ f.split('.')[-1] for f in filelistin ]))
    statname = statnameset[0]
    if len(statnameset) != 1:
        print("WARN : read_gins_multi_extracted : len(statnameset) != 1")
    fileopenedlist = [ open(f,'r+') for f in filelistin ]
    coortypelist = [ os.path.basename(f)[1] for f in filelistin ]
    print(coortypelist,filelistin)
    if 'X' in coortypelist:
        initype='XYZ'
        ia = coortypelist.index('X')
        ib = coortypelist.index('Y')
        ic = coortypelist.index('Z')
    elif 'P' in coortypelist:
        initype = 'FLH'
        ia = coortypelist.index('P')
        ib = coortypelist.index('L')
        ic = coortypelist.index('H')

    for lf in zip(*fileopenedlist):
        # OLD MODE : MAVAIS GESTION DE LA DATE
#        t1 = lf[0].split()[0]
#        t2 = lf[1].split()[0]
#        t3 = lf[2].split()[0]
#        if not (t1 == t2 == t3):
#            print "WARN : read_gins_multi_extracted : not (t1 == t2 == t3)"
#
        fields = lf[0].split()
        line = lf[0]
        blok = fields[-1]
        jour = int(blok[0:2])
        h = int(blok[2:4])
        m = int(blok[4:6])
        s = int(blok[6:8])
        yy = int(blok[-4:-2])
        # gestion des années
        if 80 < yy <= 99:
            yy = yy + 1900
        else:
            yy = yy + 2000

        # pour le mois, si > sept (9), alors lettre ...
        mm = blok[-2]

        if mm == 'O':
            mm = 10
        elif mm == 'N':
            mm = 11
        elif mm == 'D':
            mm = 12
        else:
            mm = int(mm)

        t1  =  (dt.datetime(yy,mm,jour,h,m,s) + dt.timedelta(seconds=-19))

        lfiasplit = lf[ia].split()
        lfibsplit = lf[ib].split()
        lficsplit = lf[ic].split()

        A  = float(lfiasplit[4])
        B  = float(lfibsplit[4])
        C  = float(lficsplit[4])
        sA = float(lfiasplit[5])
        sB = float(lfibsplit[5])
        sC = float(lficsplit[5])
        T = t1

        if initype == 'FLH' and flh_in_rad:
            A  = np.rad2deg(A)
            B  = np.rad2deg(B)
            sA = np.rad2deg(sA)
            sB = np.rad2deg(sB)

        pt = Point(A,B,C,T,initype,sA,sB,sC)

        tsout.add_point(pt)
    tsout.meta_set(stat=statname)
    tsout.sort()

    return tsout

def read_gins_double_diff(filein):
    """
    return a list of Point object for a double diff listing
    """

    if genefun.grep(filein , 'c o n v e r g e n c e') == '':
        print('ERR : ' , filein , 'have no convergence, return None')
        return None

    fileopened = open(filein)

    statlist    = []
    rawdatalist = []
    timelist    = []

    regstat   = re.compile(r'[0-9]{5}[A-Z][0-9]{3}  [0-9]{7} .*$')
    regresult = re.compile(r'\[S[PLHXYZ][E ].*\]     $')


    for line in fileopened:
        # finding
        if regstat.search(line):
            stat = line.split()[3]
            statlist.append(stat)

        if regresult.search(line):
            fields = line.split()
            rawdatalist.append([float(e) for e in fields[:-2]])
            timelist.append(gins_read_time(line))


    rawdatatab = np.vstack(rawdatalist)

    Pstk = []

    for istat,stat in enumerate(statlist):
        P = Point(rawdatatab[3*istat,3],
                       rawdatatab[3*istat+1,3],
                       rawdatatab[3*istat+2,3],
                       timelist[3*istat],'XYZ',
                       rawdatatab[3*istat,4],
                       rawdatatab[3*istat,4],
                       rawdatatab[3*istat,4],name=stat)

        Pstk.append(P)

    return Pstk

def read_gins_double_diff_multi(filelistin):
    """
    return a dictionnary with station names as keys and timeseries as values
    """

    Ptsstk = []
    for filein in filelistin:
        Pstk = read_gins_double_diff(filein)
        if Pstk != None:
            Ptsstk = Ptsstk + Pstk

    statlist = set([e.name for e in Ptsstk])

    tsdico = dict()

    for stat in statlist:
        tsdico[stat] = TimeSeriePoint(stat=stat)

    for pt in Ptsstk:
        tsdico[pt.name].add_point(pt)

    for ts in tsdico.values():
        ts.sort()

    return tsdico




def read_jpl_timeseries_solo(latlonrad_files_list):

    tsout = TimeSeriePoint()

    latpath = [f for f in latlonrad_files_list if '.lat' in f][0]
    lonpath = [f for f in latlonrad_files_list if '.lon' in f][0]
    radpath = [f for f in latlonrad_files_list if '.rad' in f][0]

    latfile = open(latpath)
    lonfile = open(lonpath)
    radfile = open(radpath)

    for llat , llon , lrad in zip(latfile,lonfile,radfile):
        flat = [float(e) for e in llat.split()[0:3]]
        flon = [float(e) for e in llon.split()[0:3]]
        frad = [float(e) for e in lrad.split()[0:3]]

        if not ( flat[0] == flon[0] == frad[0]):
            print(flat[0] , flon[0] , frad[0])
            raise Exception('ERR : read_jpl_timeseries_solo : Time dont corresponds !!!')

        statlat = os.path.basename(latpath).split('.')[0]
        statlon = os.path.basename(lonpath).split('.')[0]
        statrad = os.path.basename(radpath).split('.')[0]

        if not ( statlat == statlon == statrad):
            print(statlat , statlon , statrad)
            raise Exception('ERR : read_jpl_timeseries_solo : station name dont corresponds !!!')

        T = geok.year_decimal2dt(flat[0])

        N = flat[1] * 10**-2
        E = flon[1] * 10**-2
        U = frad[1] * 10**-2

        sN = flat[2] * 10**-2
        sE = flon[2] * 10**-2
        sU = frad[2] * 10**-2

        point = Point(E,N,U,T,'ENU',sE,sN,sU)

        tsout.boolENU = True
        tsout.add_point(point)

    tsout.stat = statlat

    return tsout

def read_nevada(filein,input_coords="enu"):
    """
    input_coords="enu" or "xyz"
    """

    tsout = TimeSeriePoint()

    envfile = open(filein)

    if input_coords=="enu":
        for l in envfile:
            f = l.split()

            if "site YYMMMDD" in l:
                continue
            if len(l) == 0:
                continue

            stat = f[0]

            T = geok.year_decimal2dt(float(f[2]))

            N = float(f[10])
            E = float(f[8])
            U = float(f[12])

            sN = float(f[15])
            sE = float(f[14])
            sU = float(f[16])

            point = Point(E,N,U,T,'ENU',sE,sN,sU)

            #tsout.refENU = Point()

            tsout.boolENU = True
            tsout.add_point(point)


    if input_coords=="xyz":
        for l in envfile:
            f = l.split()

            if "site YYMMMDD" in l:
                continue
            if len(l) == 0:
                continue

            stat = f[0]

            T = geok.year_decimal2dt(float(f[2]))

            X = float(f[3])
            Y = float(f[4])
            Z = float(f[5])

            sX = float(f[6])
            sY = float(f[7])
            sZ = float(f[8])

            point = Point(X,Y,Z,T,'XYZ',sX,sY,sZ)

            point.anex['Rxy'] = float(f[9])
            point.anex['Rxz'] = float(f[10])
            point.anex['Ryz'] = float(f[11])

            tsout.add_point(point)

    tsout.stat = stat

    return tsout

def read_IGS_coords(filein,initype='auto'):
    Tstk , Astk, Bstk, Cstk = [] , [] , [] , []
    tsout = TimeSeriePoint()
    for l in open(filein):
        f = l.split()

        if initype == 'auto':
            if "plh" in os.path.basename(filein):
                initype = 'FLH'
            elif "xyz" in os.path.basename(filein):
                initype = 'XYZ'

        T = geok.MJD2dt(float(f[3]))
        A = float(f[6])
        B = float(f[7])
        C = float(f[8])
        sA = float(f[9])
        sB = float(f[10])
        sC = float(f[11])

        initype = "FLH"

        pt = Point(A,B,C,T,initype,sA,sB,sC,f[0])
        tsout.add_point(pt)

    tsout.meta_set(filein,f[0])
    return tsout

def sorting_a_calais_file(openedfile):
    openedfile.seek(0)
    T , A, sA , STAT = [] , [] , [] , []
    for l in openedfile:
        f = l.split()
        T.append(float(f[0]))
        A.append(float(f[1]))
        sA.append(float(f[2]))
        STAT.append(str(f[3]))

    DATA = [ T , A, sA , STAT ]
    DATA2 = genefun.sort_table(DATA,0)

    return DATA2

def read_calais(filelist):
    """ filelistin est une liste de 3 fichier E N & U """

    if filelist == []:
        print('WARN : read_calais : files list empty , exiting ...')
        return None
    filelist.sort()

    fileopenedlist = [open(f) for f in filelist]

    sorted_data_lis = [sorting_a_calais_file(fil) for fil in fileopenedlist]

    statnameset = list(set([ f.split('.')[0] for f in filelist ]))
    statname = os.path.basename(statnameset[0])

    # LOADING ALL DATA IN A BIG MATRIX
    bigT = np.array(sorted(set(np.hstack([np.array(sd[0]) for sd in sorted_data_lis]))))
    bigDATA = np.empty((len(bigT),3))
    bigDATA.fill(np.nan)
    bigDATAsigma = np.empty((len(bigT),3))
    bigDATAsigma.fill(np.nan)

    for i,bt in enumerate(bigT):
        for j,data in enumerate(sorted_data_lis):
            for t , d ,sd, stat in zip(*data):
                if t == bt:
                    bigDATA[i,j] = d
                    bigDATAsigma[i,j] = sd

    DATA = np.hstack((bigDATA / 1000. ,bigDATAsigma / 100.))

    ptslist = []

    # MAKING POINTS
    for i in range(DATA.shape[0]):
        pt = Point(DATA[i,0],DATA[i,1],DATA[i,2],geok.convert_partial_year(bigT[i])
        ,'ENU',DATA[i,3],DATA[i,4],DATA[i,5])
        ptslist.append(pt)

    tsout = TimeSeriePoint()

    for pt in ptslist:
        tsout.add_point(pt)

    #FINDING DISCONT
    # finding discont directly in files

    # Finding the composant with max of data
    lendata = [len(sd[0]) for sd in sorted_data_lis]
    ii = lendata.index(max(lendata))

    discont = []
    T = sorted_data_lis[ii][0]
    STAT = sorted_data_lis[ii][-1]
    for i in range(len(STAT)-1):
        if STAT[i+1] != STAT[i]:
            discont.append(T[i+1])

    tsout.set_discont(discont)
    tsout.meta_set(stat=statname)
    tsout.boolENU = True
    tsout.sort()

    return tsout

def read_renag_synthetic(filein , discont_file_in = None):

    tsout = TimeSeriePoint()

    fil = open(filein)

    for l in fil:

        if l[0] == "#":
            continue

        f = l.split()

        T = geok.year_decimal2dt(float(f[0]))

        N = float(f[1])
        E = float(f[2])
        U = float(f[3])

        sN = float(f[4])
        sE = float(f[5])
        sU = float(f[6])

        point = Point(E,N,U,T,'ENU',sE,sN,sU)

        #tsout.refENU = Point()

        tsout.boolENU = True
        tsout.add_point(point)


    if discont_file_in:
        DiscontInp = open(discont_file_in)

        Discont = []
        for l in DiscontInp:
            if l[0] == "#":
                continue
            f = l.split()
            try:
                Discont.append(geok.doy2dt(int(f[0]),int(f[1])))
            except:
                print("WARN : something went wrong during discont. file reading")
                pass

        Discont = sorted(Discont)
        tsout.set_discont(Discont)

    stat_name = os.path.basename(filein).split(".")[0].split(".")[0]
    print(stat_name)
    tsout.meta_set(path=filein,stat=stat_name,name=stat_name)


    return tsout


def read_jump_file(filein,returned_events=('S','E','D')):
    """
    From a "Jump" File (P. Sakic internal file)
    Return a dictionnairy with events

    Parameters
    ----------
    filein : str
        path of the Jump File

    returned_events : tuple or list
        contains the inital letter of the event type which will be stored in the dico
        (See below)

    Returns
    -------
    jump_dico : dict of dict of datetime
        outputed events, in the form jump_dico["STAT"]["L"] = datetime
        where L is the inital letter of the event type

    Note
    ----
    A jump file contains infos like this :

    >>> STAT S 2000 001
    >>> STAT E 2001 001
    >>> STAT D 2000 06 01

    it can manage YEAR DOY or YEAR MM DD or DECIMAL YEAR

    A non-blank 1st column is a commented line

    After a #, it is a commentary

    event type letters :
        S : Start

        E : End

        D : Discontinuity

    """
    import pytz

    F = open(filein)
    jump_dico = dict()
    for l in F:
        l =  l.split("#")[0]
        f =  l.split()
        print(f)
        ## Skip comment lines
        if (not f) or (l[0] != " ") or ("#" in l):
            continue
        else:
            stat  = f[0]
            event = f[1]
            ## Create the key for the station
            if not f[0] in jump_dico.keys():
                jump_dico[stat] = dict()
                for rtn_evt in returned_events:
                    jump_dico[stat][rtn_evt] = []

            ## Fill the dico
            if event in returned_events:
                if    len(f) == 4: ### DOY
                    date = geok.doy2dt(int(f[2]),int(f[3]))
                elif  len(f) == 5: ### YYYY MM DD
                    date = dt.datetime(*[int(e) for e in f[2:]])
                elif  len(f) == 3: ### DECIMAL YEAR
                    date = geok.year_decimal2dt(float(f[2]))

                date_tz = date.replace(tzinfo=pytz.UTC)
                jump_dico[stat][event].append(date_tz)
    return jump_dico

def read_nav_step1_geodesea(filein):
    M = np.loadtxt(filein)
    tsout = TimeSeriePoint()

    for m in M:
        pt = Point(np.rad2deg(m[1]),np.rad2deg(m[2]),m[3],m[0],initype='FLH')
        tsout.add_point(pt)

    return tsout


def read_nrcan_csv(filein , associated_ps_file = '', statname = ''):
    """
    associated_ps_file is highly recommanded
    because of the time managing

    WARN : Must be avoided b/c of the weak decimal precision of the angles !!!
    """

    if statname == '':
        statname = os.path.basename(filein)[0:4]

    pdcsv = pd.read_csv(filein)

    F     = np.array(pdcsv['latitude_degre_decimal'])
    L     = np.array(pdcsv['longitude_degre_decimal'])
    H     = np.array(pdcsv['hauteur_ellipsoidale_m'])
    heure = np.array(pdcsv['heure_decimal'])
    doy   = np.array(pdcsv['jour_de_l_annee'])
    year  = np.array(pdcsv['annee'])

    T = geok.doy2dt(year,doy,heure)

    if associated_ps_file != '':
        T = []
        for l in open(associated_ps_file):
            if "BWD" in l:
                f = l.split()
                t = geok.date_string_2_dt(f[4] + ' ' + f[5])
                T.append(t)
        if statname == '':
            statname = f[2]

    tsout = TimeSeriePoint()
    for f,l,h,t in zip(F,L,H,T):
        point = Point(f,l,h,t,'FLH',name = statname)
        tsout.add_point(point)

    tsout.meta_set(filein,statname)

    return tsout

def read_nrcan_pos(filein):
    """
    .pos file are more precise than .csv, should be used !
    """
    tsout = TimeSeriePoint()
    start_read = False

    for l in open(filein):
        if l[0:3] == 'DIR':
            start_read = True
            continue
        elif not start_read:
            continue
        else:
            f = l.split()
            lat = (np.abs(float(f[20])) + 1/60. * float(f[21]) + 1/3600. * float(f[22])) * np.sign(float(f[20]))
            lon = (np.abs(float(f[23])) + 1/60. * float(f[24]) + 1/3600. * float(f[25])) * np.sign(float(f[23]))
            sE = float(f[15])
            sN = float(f[16])
            sU = float(f[17])
            h   = float(f[26])
            slat , slon , sh = geok.sENU2sFLH(lat,lon,h,sE,sN,sU)

            t   = geok.date_string_2_dt(f[4] + ' ' + f[5])

            pt = Point(lat,lon,h,t,'FLH',slat,slon,sh,name = f[2])
            tsout.add_point(pt)

    tsout.meta_set(filein,f[2])
    return tsout


def read_qinsy(filein,yy,mm,dd):
    reader = pd.read_csv(open(filein))
    T = [ dateutil.parser.parse(e).replace(year=yy, month=mm, day=dd) + dt.timedelta(seconds=dUTCGPS) for e in list(reader.icol(0))]
    (X,Y,Z) = np.array(geok.GEO2XYZ(reader.icol(12),reader.icol(13),reader.icol(14)))
    #sX,sY,sZ = [] , [] , []
    initype = 'XYZ'
    tsout = TimeSeriePoint()
    for i in range(len(T)):
        point = Point(X[i],Y[i],Z[i],T[i],initype,sA=np.nan,sB=np.nan,sC=np.nan)
        tsout.add_point(point)
    tsout.meta_set(filein)
    return tsout

def read_sonardyne_posi(filein):
    reader = pd.read_csv(open(filein),skip_footer=1)
    T = [ dateutil.parser.parse(e) + dt.timedelta(seconds=dUTCGPS) for e in list(reader['UTCTime'])]
    (L,F,H) = (reader['Longitude'],reader['Latitude'],reader['Altitude'])
    sX,sY,sZ = [] , [] , []
    initype = 'FLH'
    tsout = TimeSeriePoint()
    for i in range(len(T)):
        point = Point(F[i],L[i],H[i],T[i],initype,sA=np.nan,sB=np.nan,sC=np.nan)
        tsout.add_point(point)

    tsout.meta_set(filein)
    return tsout

def read_pbo_pos(filein):
    filobj = open(filein,'r')
    tsout = TimeSeriePoint()
    header = True
    for line in filobj:
        if not header:
            f = line.split()
            f2 = [float(e) for e in f[:-1]]
            t = dt.datetime(int(f[0][0:4]),int(f[0][4:6]),int(f[0][6:]),int(f[1][0:2]),int(f[1][2:4]),int(f[1][4:]))
            pt = Point(f2[3],f2[4],f2[5],t,'XYZ',f2[6],f2[7],f2[8])
            pt.FLHset(f2[12],f2[13],f2[14])
            pt.ENUset(f2[16],f2[15],f2[17],f2[19],f2[18],f2[20])
            pt.anex['Rxy'] = f2[9]
            pt.anex['Rxz'] = f2[10]
            pt.anex['Ryz'] = f2[11]
            pt.anex['Rne'] = f2[-4]
            pt.anex['Rnu'] = f2[-3]
            pt.anex['Reu'] = f2[-2]
            tsout.add_point(pt)
        if line[0] == '*':
            header = False
    tsout.boolENU = True
    tsout.meta_set(filein)
    tsout.stat = os.path.basename(filein)[0:4]
    return tsout


def read_sonardyne_attitude(filein):
    reader = pd.read_csv(open(filein),skip_footer=1)

    T = [ dateutil.parser.parse(e) + dt.timedelta(seconds=dUTCGPS) for e in list(reader['PTPTime'])]

    pitch = np.array(reader['Pitch'])
    roll  = np.array(reader['Roll'])
    head  = np.array(reader['Heading'])

    initype = 'FLH'

    AHRS = list(reader['AHRS Id'])
    nbunit = geok.guess_seq_len(AHRS)
    print("nbre de devices ID : " ,nbunit)
    AHRS = reader['AHRS Id']

    tsout_list = []

    for n in range(nbunit):
        tsout = TimeSerieObs()
        tsout.typeobs = 'RPY'
        tsout_list.append(tsout)

    i = 0
    for p,r,h,t,a in zip(roll,pitch,head,T,AHRS):
        att = Attitude(p,r,h,t,devID=a)
        n = np.mod(i,nbunit)
        tsout_list[n].add_obs(att)
        i=i+1

    for ts in tsout_list:
        ts.meta_set(filein,devID=list(set([e.devID for e in ts.obs])))
        ts.interp_set()

    return tsout_list


def interp_sndy_SYS_UTC(time_conv_file_in):
    timeConv = genefun.read_mat_file(timepath)
    timeSYS  = timeConv[0,:]
    timeUTC  = timeConv[1,:]
    IntSYSUTC = scipy.interpolate.interp1d(timeSYS,timeUTC,kind='slinear',
                                           bounds_error=0)
    TSout = TimeSerieObs(time_conv_file_in)

    return IntSYSUTC


def read_sndy_mat_att(filein,IntSYSUTCin=None):
    if IntSYSUTCin == None:
        print("WARN : read_sndy_mat_att : No Interpolator")
    attmat  = genefun.read_mat_file(filein)
    Tsys_att = attmat[0,:]
    Tposix_att = IntSYSUTCin(Tsys_att)
#    Tdt_att = np.array(geok.posix2dt(Tposix_att))
    roll  = attmat[1,:]
    pitch = attmat[2,:]
    head  = attmat[3,:]

    TSout = TimeSerieObs('RPY',filein)

    for r,p,h,t in zip(roll,pitch,head,Tposix_att):
        att = Attitude(r,p,h,t)
        TSout.add_obs(att)
    return TSout


def read_sndy_mat_nav(filein,IntSYSUTCin=None):
    if IntSYSUTCin == None:
        print("WARN : read_sndy_mat_att : No Interpolator")
    navmat  = genefun.read_mat_file(filein)
    Tsys_nav = navmat[0,:]
    Tposix_nav = IntSYSUTCin(Tsys_nav)
#    Tdt_att = np.array(geok.posix2dt(Tposix_att))
    lat = navmat[2,:]
    lon = navmat[3,:]
    h   = navmat[4,:]

    TSout = TimeSeriePoint()

    for f,l,h,t in zip(lat,lon,h,Tposix_nav):
        pt = Point(f,l,h,t,initype='FLH')
        TSout.add_point(pt)
    return TSout

def read_hector_neu(filein):
    print("WARN : XYZ/FLH conversion not implemented")
    M = np.loadtxt(filein)
    stat = genefun.grep(filein,'Site :',only_first_occur=True).split()[3]
    tsout = ts_from_list(M[:,2],M[:,1],M[:,3],geok.year_decimal2dt(M[:,0]),
                         'ENU',M[:,4],M[:,5],M[:,6],stat=stat,name=stat)

    return tsout


def read_bull_B(path):

    if not genefun.is_iterable(path):
        path = [path]

    path = sorted(path)

    DFstk = []

    for path_solo in path:
        S = genefun.extract_text_between_elements(path_solo,"1 - DAILY FINAL VALUES" ,
                                         "2 - DAILY FINAL VALUES" )

        L = S.replace('\t','\n').split("\n")

        L2 = []
        for e in L:
            if len(e) > 0:
                if e[0] !=  " ":
                    L2.append(e)
        L3 = []
        for e in L2:
            L4 = []
            for ee in e.split():
                L4.append(float(ee))
            L3.append(L4)

        DF = pd.DataFrame(np.vstack(L3))
        DFstk.append(DF)

    DFout = pd.concat(DFstk)
    DFout.columns = ["year","month","day","MJD","x","y","UT1-UTC","dX","dY",
                  "x err","y err","UT1 err","X err","Y err"]

    return DFout
##
def read_erp(path,return_array=False):
    """
    This function is discontinued, use read_erp1 instead
    """
    F = open(path)
    L = []
    for l in F:
        if "MJD" in l:
            keys = l.split()
        try:
            float(l[0])
        except:
            continue
        l=genefun.str_2_float_line(l.strip())
        L.append(l)

    M = np.vstack(L)

    if return_array:
        return M
    else:
        return pd.DataFrame(M)

def read_erp_multi(path_list , return_array=False,
                   smart_mode=True):
    """
    This function is discontinued, use read_erp1 instead

    Input :
        path_list : a list of ERP files
        smart_mode : keep only the latest value (True is recommended)
    """
    path_list = sorted(path_list)
    Lstk = []
    for path in path_list:
        L = read_erp(path,return_array)
        Lstk.append(L)

    M = np.vstack(Lstk)

    if smart_mode:
        Msmart = np.array([])
        for ilm , lm in enumerate(np.flipud(M)):
            if ilm == 0:
                Msmart = np.array([lm])
            elif lm[0] in Msmart[:,0]:
                continue
            else:
                Msmart = np.vstack((Msmart,lm))

        M = np.flipud(Msmart)

    if return_array:
        return M
    else:
        return pd.DataFrame(M)

def read_sp3(file_path_in,returns_pandas = True, name = '',
             epoch_as_pd_index = False):
    """
    Read a SP3 file (GNSS Orbits standard file) and return X,Y,Z coordinates
    for each satellite and for each epoch

    Parameters
    ----------
    file_path_in : str
        path of the SP3 file

    returns_pandas : bool
        if True, return a Pandas DataFrame.
        if False, return separated lists.

    name : str
        a manual name for the file

    epoch_as_pd_index : bool
        if True, the index of the output dataframe contains
        if False, it contains generic integer indexs


    Returns
    -------
    df : Pandas DataFrame
        if returns_pandas == True

    epoch_stk ,  Xstk , Ystk , Zstk , Clkstk : lists
        if returns_pandas == False

    """

    AC_name =  os.path.basename(file_path_in)[:3]

    fil = open(file_path_in)

    header = True

    epoch_stk = []
    Xstk , Ystk , Zstk , Clkstk = [],[],[],[]

    data_stk  = []

    for l in fil:
        if l[0] == '*':
            header = False

        if header:
            continue
        if 'EOF' in l:
            continue

        if l[0] == '*':
            epoc   = geok.tup_or_lis2dt(l[1:].strip().split())
        else:
            sat_nat = l[1:2].strip()
            sat_sv  = int(l[2:4].strip())
            sat_sat = l[1:4].strip()

            X   = float(l[4:18])
            Y   = float(l[18:32])
            Z   = float(l[32:46])
            Clk = float(l[46:60])

            if returns_pandas:
                line_data = [epoc,sat_sat,sat_nat,sat_sv,X,Y,Z,Clk,AC_name]
                data_stk.append(line_data)
            else:
                epoch_stk.append(epoc)
                Xstk.append(X)
                Ystk.append(Y)
                Zstk.append(Z)
                Clkstk.append(Clk)


    AC_name_stk = [AC_name] * len(Xstk)

    if returns_pandas:
        df = pd.DataFrame(data_stk, columns=['epoch','sat', 'const', 'sv',
                                             'x','y','z','clk','AC'])
        if epoch_as_pd_index:
            df.set_index('epoch',inplace=True)
        df.filename = os.path.basename(file_path_in)
        df.path = file_path_in

        if name != '':
            df.name = name
        else:
            df.name = os.path.basename(file_path_in)

        return df
    else:
        print("INFO : return list, very beta : no Sat. Vehicule Number info ...")
        return  epoch_stk ,  Xstk , Ystk , Zstk , Clkstk , AC_name_stk


def read_sp3_header(sp3_path):
    """
    Read a SP3 file header and return a Pandas DataFrame
    with sat. PRNs and sigmas contained in the header

    Parameters
    ----------
    sp3_path : str
        path of the SP3 file


    Returns
    -------
    df : Pandas DataFrame
        2 columns "sat", "sigma"

    Note
    -------
    More infos about the sigma
    http://acc.igs.org/orbacc.txt
    """


    F = open(sp3_path)
    ac_name = os.path.basename(sp3_path)[:3]


    Lines = F.readlines()

    Sat_prn_list = []
    Sat_sig_list = []

    for il , l in enumerate(Lines):
        if il == 1:
            date = geok.MJD2dt(int(l.split()[4]))
        if l[:2] == "+ ":
            Sat_prn_list.append(l)
        if l[:2] == "++":
            Sat_sig_list.append(l)
        if l[0] == "*":
            break

    ### PRN part
    Sat_prn_list_clean = []
    for prn_line in Sat_prn_list:
        prn_line_splited = prn_line.split()
        prn_line_splited = [e for e in prn_line_splited if not "+" in e]
        prn_line_splited = [e for e in prn_line_splited if not  e == "0"]
        Sat_prn_list_clean = Sat_prn_list_clean + prn_line_splited

    sat_nbr = int(Sat_prn_list_clean[0])

    Sat_prn_list_clean = Sat_prn_list_clean[1:]

    Sat_prn_string = "".join(Sat_prn_list_clean)

    Sat_prn_list_final = []
    for i in range(sat_nbr):
        Sat_prn_list_final.append(Sat_prn_string[i*3:i*3+3])

    ### Sigma part
    Sat_sig_list_clean = []
    for sig_line in Sat_sig_list:
        sig_line_splited = sig_line.split()
        sig_line_splited = [e for e in sig_line_splited if not "+" in e]
        Sat_sig_list_clean = Sat_sig_list_clean + sig_line_splited

    Sat_sig_list_final = [int(e) for e in Sat_sig_list_clean[:sat_nbr]]


    ### Export part
    AC_name_list = [ac_name] * sat_nbr
    Date_list    = [date] * sat_nbr

    Header_DF = pd.DataFrame(list(zip(AC_name_list,Sat_prn_list_final,
                                      Sat_sig_list_final,Date_list)),
                             columns=["AC","sat","sigma","epoch"])

    return Header_DF



def sp3_decimate(file_in,file_out,step=15):

    Fin = open(file_in)

    good_line = True
    outline = []

    for l in Fin:
        if l[0] == "*":
            epoc   = geok.tup_or_lis2dt(l[1:].strip().split())
            if np.mod(epoc.minute , step) == 0:
                good_line = True
            else:
                good_line = False

        if good_line:
            outline.append(l)

    line1     = outline[1]
    step_orig = outline[1].split()[3]

    step_ok = "{:14.8f}".format(step * 60).strip()

    line1b = line1.replace(step_orig,step_ok)
    outline[1] = line1b


    with open(file_out,"w+") as Fout:
        for l in outline:
            Fout.write(l)

    return file_out




def write_sp3(SP3_DF_in , outpath):
    """
    Write DOCSTRING
    """
    ################## MAIN DATA
    LinesStk = []


    SP3_DF_in.sort_values(["epoch","sat"],inplace=True)

    EpochList  = SP3_DF_in["epoch"].unique()
    SatList    = sorted(SP3_DF_in["sat"].unique())
    SatListSet = set(SatList)

    for epoc in EpochList:
        SP3epoc   = pd.DataFrame(SP3_DF_in[SP3_DF_in["epoch"] == epoc])
        ## Missing Sat
        MissingSats = SatListSet.difference(set(SP3epoc["sat"]))
        
        for miss_sat in MissingSats:
            miss_line = SP3epoc.iloc[0].copy()
            miss_line["sat"]   = miss_sat
            miss_line["const"] = miss_sat[0]
            miss_line["x"]     = 0.000000
            miss_line["y"]     = 0.000000
            miss_line["z"]     = 0.000000
            miss_line["clk"]   = 999999.999999
            
            SP3epoc = SP3epoc.append(miss_line)

        SP3epoc.sort_values("sat",inplace=True)
        timestamp = geok.dt2sp3_timestamp(geok.numpy_dt2dt(epoc)) + "\n"

        LinesStk.append(timestamp)

        linefmt = "P{:}{:14.6f}{:14.6f}{:14.6f}{:14.6f}\n"


        for ilin , lin in SP3epoc.iterrows():
            line_out = linefmt.format(lin["sat"],lin["x"],lin["y"],lin["z"],lin["clk"])

            LinesStk.append(line_out)



    ################## HEADER
    ######### SATELLITE LIST

    Satline_stk   = []
    Sigmaline_stk = []

    for i in range(5):
        SatLine = SatList[17*i:17*(i+1)]
        if len(SatLine) < 17:
            complem = " 00" * (17 - len(SatLine))
        else:
            complem = ""

        if i == 0:
            nbsat4line = len(SatList)
        else:
            nbsat4line = ''

        satline = "+  {:3}   ".format(nbsat4line) + "".join(SatLine) + complem + "\n"
        sigmaline = "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"

        Satline_stk.append(satline)
        Sigmaline_stk.append(sigmaline)


    ######### 2 First LINES
    start_dt = geok.numpy_dt2dt(EpochList.min())

    header_line1 = "#cP" + geok.dt2sp3_timestamp(start_dt,False) + "     {:3}".format(len(EpochList)) + "   u+U IGSXX FIT  XXX\n"

    delta_epoch = int(genefun.most_common(np.diff(EpochList) * 10**-9))
    MJD  = geok.dt2MJD(start_dt)
    MJD_int = int(np.floor(MJD))
    MJD_dec = MJD - MJD_int
    gps_wwww , gps_sec = geok.dt2gpstime(start_dt,False,"gps")

    header_line2 = "## {:4} {:15.8f} {:14.8f} {:5} {:15.13f}\n".format(gps_wwww,gps_sec,delta_epoch,MJD_int,MJD_dec)


    ######### HEADER BOTTOM
    header_bottom = """%c G  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
%f  1.2500000  1.025000000  0.00000000000  0.000000000000000
%f  0.0000000  0.000000000  0.00000000000  0.000000000000000
%i    0    0    0    0      0      0      0      0         0
%i    0    0    0    0      0      0      0      0         0
/* PCV:IGSXX_XXXX OL/AL:FESXXXX  NONE     YN CLK:CoN ORB:CoN
/*     GeodeZYX Toolbox Output
/*
/*
"""


    ################## FINAL STACK

    FinalLinesStk = []

    FinalLinesStk.append(header_line1)
    FinalLinesStk.append(header_line2)
    FinalLinesStk = FinalLinesStk + Satline_stk + Sigmaline_stk
    FinalLinesStk.append(header_bottom)
    FinalLinesStk = FinalLinesStk + LinesStk + ["EOF"]

    FinalStr = "".join(FinalLinesStk)

    F = open(outpath,"w+")
    F.write(FinalStr)

def clk_decimate(file_in,file_out,step=5):

    Fin = open(file_in)

    good_line = True
    outline = []

    step = 5

    for l in Fin:
        good_line = True
        if l[0:2] in ("AR","AS"):
            epoc   = geok.tup_or_lis2dt(l[8:34].strip().split())
            if np.mod(epoc.minute , step) == 0:
                good_line = True
            else:
                good_line = False

        if good_line:
            outline.append(l)

    with open(file_out,"w+") as Fout:
        for l in outline:
            Fout.write(l)

    return file_out


def AC_equiv_vals(AC1,AC2):
    ### 1) Merge the 2 DF to find common lines
    ACmerged = pd.merge(AC1 , AC2 , how='inner', on=['epoch', 'sat'])

    ### 2) Extract merged epoch & sv
    common_epoch = ACmerged["epoch"]
    common_sat   = ACmerged["sat"]
    ### 3) Create a boolean line based on common epoch / sv
    common_sat_epoch_AC1 = (AC1["epoch"].isin(common_epoch)) & (AC1["sat"].isin(common_sat))
    common_sat_epoch_AC2 = (AC2["epoch"].isin(common_epoch)) & (AC2["sat"].isin(common_sat))
    ### 4) Get epoch and sv in the combined sol which correspond to the SP3
    AC1new = AC1[common_sat_epoch_AC1].copy()
    AC2new = AC2[common_sat_epoch_AC2].copy()
    ### 5) A sort to compare the same things
    AC1new.sort_values(by=['sat','sv'],inplace=True)
    AC2new.sort_values(by=['sat','sv'],inplace=True)

    ### Check for > 99999 vals
    AC1_bad_bool_9  = (AC1new["x"] > 9999) & (AC1new["y"] > 9999) & (AC1new["z"] > 9999)
    AC1_bad_bool_9  = np.logical_not(np.array(AC1_bad_bool_9))

    AC2_bad_bool_9  = (AC2new["x"] > 9999) & (AC2new["y"] > 9999) & (AC2new["z"] > 9999)
    AC2_bad_bool_9  = np.logical_not(np.array(AC2_bad_bool_9))

    AC12_bad_bool_9 = np.array(np.logical_and(AC1_bad_bool_9 , AC2_bad_bool_9))

    ### Check for NaN vals
    AC1_bad_bool_nan  = np.isnan(AC1new["x"]) & np.isnan(AC1new["y"]) & np.isnan(AC1new["z"])
    AC1_bad_bool_nan  = np.logical_not(np.array(AC1_bad_bool_nan))

    AC2_bad_bool_nan  = np.isnan(AC2new["x"]) & np.isnan(AC2new["y"]) & np.isnan(AC2new["z"])
    AC2_bad_bool_nan  = np.logical_not(np.array(AC2_bad_bool_nan))

    AC12_bad_bool_nan = np.array(np.logical_and(AC1_bad_bool_nan , AC2_bad_bool_nan))

    AC12_bad_bool = np.logical_and(AC12_bad_bool_9, AC12_bad_bool_nan)

    AC1_ok = AC1new[AC12_bad_bool]
    AC2_ok = AC2new[AC12_bad_bool]

    return AC1_ok , AC2_ok , ACmerged


############################# READ CLK 30
def read_clk1(file_path_in, returns_pandas = True, interval=None):
    """
    General description

    Parameters
    ----------
    file_path_in :  str
        Path of the file in the local machine.

    returns_pandas :  bool
        Define if pandas function will be used to put the data on tables

    interval :  int
        Defines which interval should be used in the data tables. The interval is always in minutes unit.

    Returns
    -------
    out1 : float or int oandas table
        Returns a panda table format with the data extracted from the file.

    out2 :  list
       In case of the pandas function are not in use, the return will be a list with the data extract from the file.


    """

    file = open(file_path_in,'r')
    fil = file.readlines()
    file.close()


    types_m, name_m, year_m, month_m, day_m, hr_m, minute_m, seconds_m, epoch_m, offset_m, rms_m = [],[],[],[],[],[],[],[],[],[],[]

    Clk_read  = []

        #####IGNORES HEADER

    le = 0
    i = 0
    count = 0
    for le in range(len(fil)):
     linhaatual = linecache.getline(file_path_in, le)
     header_end = (linhaatual[60:73])
     count +=1
     if header_end =="END OF HEADER":
        i = count
#############################################################
    while i <= len(fil):
          linhaatual = linecache.getline(file_path_in, i)
          types = (linhaatual[0:2]) # Column refers to AS or AR
          name = (linhaatual[3:7])# Name of the station or satellite
          year = int((linhaatual[8:12]))
          month = int((linhaatual[13:15]))
          day = int((linhaatual[16:18]))
          hr= int((linhaatual[19:22]))# hour
          minute = int((linhaatual[22:25]))
          seconds = float((linhaatual[25:33]))
          offset = float(((linhaatual[40:59])))# Clock offset
          ##### check if there is a value for thr rms
          if (linhaatual[65:70]) == '     ' or len(linhaatual[65:70])==0:
              rms = 0
          else:
              rms = float(linhaatual[61:79])
          ############
          epoch = dt.datetime(year,month,day,hr,minute,int(seconds),int(((int(seconds)-seconds)*10**-6))) # epoch as date.time function

        ############ Data in a defined interval
          if interval:
              if (float(minute)%interval) == 0 and float(seconds) == 0:

                  if returns_pandas:
                            Clk_dados = [types,name,year,month,day,hr,minute,seconds,epoch,offset,rms]
                            Clk_read.append(Clk_dados)
                  else:
                            types_m.append(types)
                            name_m.append(name)
                            year_m.append(year)
                            month_m.append(month)
                            day_m.append(day)
                            hr_m.append(hr)
                            minute_m.append(minute)
                            seconds_m.append(seconds)
                            epoch_m.append(epoch)
                            offset_m.append(offset)
                            rms_m.append(rms)
         #########################################################
        ########################################################### standart file epochs
          else:
              if returns_pandas:
                            Clk_dados = [types,name,year,month,day,hr,minute,seconds,epoch,offset,rms]
                            Clk_read.append(Clk_dados)
              else:
                            types_m.append(types)
                            name_m.append(name)
                            year_m.append(year)
                            month_m.append(month)
                            day_m.append(day)
                            hr_m.append(hr)
                            minute_m.append(minute)
                            seconds_m.append(seconds)
                            epoch_m.append(epoch)
                            offset_m.append(offset)
                            rms_m.append(rms)


          i +=1

  ################################################################################################
   ############################################ put on pandas table format
    if returns_pandas:
     Clk_readed = pd.DataFrame(Clk_read, columns=['type','name','year','month','day','hr','minutes','seconds','epoch','offset','rms'])
     Clk_readed.path = file_path_in

     return Clk_readed

          ###############################
    else:
           print(" ...")
           return  types_m, name_m, year_m, month_m, day_m, hr_m, minute_m, seconds_m, offset_m, rms_m




############################# END FUNCTION READ CLK 30




################################################READ ERP GUS
def read_erp1(caminho_arq,ac):

    le = open(caminho_arq, 'r')
    letudo = le.readlines()
    le.close()
    tamanho = len(letudo) #usado para saber quantas linhas tem o arquivo



    para = tamanho #onde o arquivo para de ser lido

    numeros = ['0','1','2','3','4','5','6','7','8','9']
    le = 0
    numlin = 0 #numero de linhas da matriz de epocas
    numcol = 16 #numero de colunas que a matriz final deve ter


    while le <= para:
     linhaatual = linecache.getline(caminho_arq, le)
     if linhaatual[0:1] in numeros:
         numlin +=1
     le +=1

    dt = np.dtype(object)
    ERP = np.empty((numlin,numcol), dtype=dt)
    n = 0
    j = 0
    l = 0
    g = 0
#####################################
    if ac == 'wum':
        while l <=para:
         linhaatual = linecache.getline(caminho_arq, l)
         if linhaatual[0:1] in numeros:
             ERP[j,0] = float(linhaatual[0:8])
             ERP[j,1] = round((float(linhaatual[11:17])*(10**-6)),6)
             ERP[j,2] = float(linhaatual[18:25])*(10**-6)
             ERP[j,4] = float(linhaatual[26:34])*(10**-7)
             ERP[j,3] = float(linhaatual[35:41])*(10**-7)
             ERP[j,5] = float(linhaatual[44:47])*(10**-6)
             ERP[j,6] = float(linhaatual[50:53])*(10**-6)
             ERP[j,7] = float(linhaatual[57:61])*(10**-7)
             ERP[j,8] = float(linhaatual[63:69])*(10**-7)
             ERP[j,9] = float(linhaatual[69:73])
             ERP[j,10] = float(linhaatual[74:76])
             ERP[j,11] = float(linhaatual[76:80])
             ERP[j,12] = float(linhaatual[81:86])*(10**-6)
             ERP[j,13] = float(linhaatual[88:93])*(10**-6)
             ERP[j,14] = float(linhaatual[96:101])*(10**-6)
             ERP[j,15] = float(linhaatual[102:107])*(10**-6)

             j +=1
         l +=1

        Erp_end = pd.DataFrame(ERP, columns=['MJD','X-P', 'Y-P', 'UT1UTC','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
                                             'X-RT','Y-RT','S-XR','S-YR'])
    #    Erp_end.set_index('MJD',inplace=True)
        return Erp_end
#######################################
    if ac == 'COD':
        while n <=para:
         linhaatual = linecache.getline(caminho_arq, n)
         if linhaatual[0:1] in numeros:
             ERP[j,0] = float(linhaatual[0:8])
             ERP[j,1] = round((float(linhaatual[12:17])*(10**-6)),6)
             ERP[j,2] = float(linhaatual[20:27])*(10**-6)
             ERP[j,4] = float(linhaatual[38:41])*(10**-7)
             ERP[j,3] = float(linhaatual[29:35])*(10**-7)
             ERP[j,5] = float(linhaatual[47:49])*(10**-6)
             ERP[j,6] = float(linhaatual[53:55])*(10**-6)
             ERP[j,7] = float(linhaatual[59:61])*(10**-7)
             ERP[j,8] = float(linhaatual[65:66])*(10**-7)
             ERP[j,9] = float(linhaatual[67:70])
             ERP[j,10] = float(linhaatual[71:74])
             ERP[j,11] = float(linhaatual[75:77])
             ERP[j,12] = float(linhaatual[81:85])*(10**-6)
             ERP[j,13] = float(linhaatual[88:92])*(10**-6)
             ERP[j,14] = float(linhaatual[94:98])*(10**-6)
             ERP[j,15] = float(linhaatual[100:104])*(10**-6)

             j +=1
         n +=1

        Erp_end = pd.DataFrame(ERP, columns=['MJD','X-P', 'Y-P', 'UT1UTC','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
                                             'X-RT','Y-RT','S-XR','S-YR'])
    #    Erp_end.set_index('MJD',inplace=True)
        return Erp_end
################################################
    if ac == 'gbm':
        while g <=para:
         linhaatual = linecache.getline(caminho_arq, g)
         if linhaatual[0:1] in numeros:
             ERP[j,0] = float(linhaatual[0:9])
             ERP[j,1] = round((float(linhaatual[12:18])*(10**-6)),6)
             ERP[j,2] = float(linhaatual[19:25])*(10**-6)
             ERP[j,4] = float(linhaatual[26:37])*(10**-7)  ###UT1 - TAI!!!!!!!!!!!!
             ERP[j,3] = float(linhaatual[39:44])*(10**-7)
             ERP[j,5] = float(linhaatual[46:49])*(10**-6)
             ERP[j,6] = float(linhaatual[50:54])*(10**-6)
             ERP[j,7] = float(linhaatual[56:59])*(10**-7)
             ERP[j,8] = float(linhaatual[60:64])*(10**-7)
             ERP[j,9] = float(linhaatual[65:70])
             ERP[j,10] = float(linhaatual[71:72])
             ERP[j,11] = float(linhaatual[73:76])
             ERP[j,12] = float(linhaatual[77:82])*(10**-6)
             ERP[j,13] = float(linhaatual[83:88])*(10**-6)
             ERP[j,14] = float(linhaatual[92:96])*(10**-6)
             ERP[j,15] = float(linhaatual[100:104])*(10**-6)

             j +=1
         g +=1

        Erp_end = pd.DataFrame(ERP, columns=['MJD','X-P', 'Y-P', 'UT1UTC','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
                                             'X-RT','Y-RT','S-XR','S-YR'])
    #    Erp_end.set_index('MJD',inplace=True)
        return Erp_end
        #out = '/home/mansur/Documents/MGEX/Saida/outERP_'+caminho_arq[28:-4]+'.txt'
        #np.savetxt(out, ERP, fmt="%02s %04s %05s %05s %05s %05s %05s %10s %10s %10s %10s %10s %10s %10s %10s %10s")
################################################### END READ ERP GUS
##########################################################################################################################
def read_erp2(caminho_arq,ac):
    """
    General description
    Units: ('MJD','X-P (arcsec)', 'Y-P (arcsec)', 'UT1UTC (E-7S)','LOD (E-7S/D)','S-X (E-6" arcsec)','S-Y (E-6" arcsec)',
    'S-UT (E-7S)','S-LD (E-7S/D)','NR (E-6" arcsec)', 'NF (E-6" arcsec)', 'NT (E-6" arcsec)',
    'X-RT (arcsec/D)','Y-RT (arcsec/D)','S-XR (E-6" arcsec/D)','S-YR (E-6" arcsec/D)', 'C-XY', 'C-XT',
    'C-YT', 'DPSI', 'DEPS','S-DP','S-DE')

    Parameters
    ----------
    file_path_in :  str
        Path of the file in the local machine.

    which AC :  str
        The analisys center that will be used


    Returns
    -------
    out1 :  pandas table
        Returns a panda table format with the data extracted from the file.


    """

    #### FIND DELIVERY DATE
    name = os.path.basename(caminho_arq)

    if len(name) == 12:
        dt_delivery = geok.sp3name2dt(caminho_arq)
    elif len(name) == 38:
        dt_delivery = geok.sp3name_v3_2dt(caminho_arq)
    else:
        dt_delivery = geok.posix2dt(0)


    le = open(caminho_arq, 'r')
    letudo = le.readlines()
    le.close()
    tamanho = len(letudo) #usado para saber quantas linhas tem o arquivo

    #para = tamanho #onde o arquivo para de ser lido

    numeros = ['0','1','2','3','4','5','6','7','8','9']
    #le = 0
    #numlin = 0 #numero de linhas da matriz de epocas
    #numcol = 16 #numero de colunas que a matriz final deve ter

    ERP=[]


    if caminho_arq[-3:] in ('snx','ssc'):
        file = open(caminho_arq)
        Lines =  file.readlines()
        XPO_stk  = []
        XPO_std_stk = []
        YPO_stk  = []
        YPO_std_stk = []
        LOD_stk  = []
        LOD_std_stk = []
        MJD_stk = []
        marker = False

        for i in range(len(Lines)):

            if len(Lines[i].strip()) == 0:
                continue
            else:

                if Lines[i].split()[0] == '+SOLUTION/ESTIMATE':
                    marker = True

                if Lines[i].split()[0] == '-SOLUTION/ESTIMATE':
                    marker = False

                if cmg.contains_word(Lines[i],'XPO') and marker:
                    Doy = (Lines[i][30:33])
                    Year = (Lines[i][27:29])
                    Pref_year = '20'
                    Year = int(Pref_year+Year)
                    Date = geok.doy2dt(Year,Doy)
                    XPO = float(Lines[i][47:68])*(10**-3)
                    XPO_std = float(Lines[i][69:80])*(10**-3)
                    XPO_stk.append(XPO)
                    XPO_std_stk.append(XPO_std)
                    MJD_stk.append(cmg.jd_to_mjd(cmg.date_to_jd(Date.year,Date.month,Date.day)))

                if cmg.contains_word(Lines[i],'YPO') and marker:
                    Doy = (Lines[i][30:33])
                    Year = str(Lines[i][27:29])
                    Pref_year = '20'
                    Year = int(Pref_year+Year)
                    Date = geok.doy2dt(Year,Doy)
                    YPO = float(Lines[i][47:68])*(10**-3)
                    YPO_std = float(Lines[i][69:80])*(10**-3)
                    YPO_stk.append(YPO)
                    YPO_std_stk.append(YPO_std)
                    MJD_stk.append(cmg.jd_to_mjd(cmg.date_to_jd(Date.year,Date.month,Date.day)))

                if cmg.contains_word(Lines[i],'LOD') and marker:
                    Doy = (Lines[i][30:33])
                    Year = str(Lines[i][27:29])
                    Pref_year = '20'
                    Year = int(Pref_year+Year)
                    Date = geok.doy2dt(Year,Doy)
                    LOD = float(Lines[i][47:68])*(10**+4)
                    LOD_std = float(Lines[i][69:80])*(10**+4)
                    LOD_stk.append(LOD)
                    LOD_std_stk.append(LOD_std)
                    MJD_stk.append(cmg.jd_to_mjd(cmg.date_to_jd(Date.year,Date.month,Date.day)))

        MJD = list(set(MJD_stk))
        if len(LOD_stk) == 0:
                LOD_stk = ['0']*len(MJD)
                LOD_std_stk = ['0']*len(MJD)

        for i in range(len(MJD)):

            ERP_data = [ac, MJD[i], XPO_stk[i], YPO_stk[i], 0, LOD_stk[i], XPO_std_stk[i], YPO_std_stk[i],
                     0, LOD_std_stk[i], 0, 0, 0, 0, 0, 0, 0, dt_delivery]
            ERP.append(ERP_data)



    if ac in ('COD','cod','com', 'cof', 'grg', 'mit', 'sio'):
        for i in range(tamanho+1):
            linhaatual = linecache.getline(caminho_arq, i)
            if linhaatual[0:1] in numeros:
                ERP_data = linhaatual.split()
                for j in range(len(ERP_data)):
                    ERP_data[j] = float(ERP_data[j])
                ERP_data.insert(0,ac)
                ERP_data[2] =  ERP_data[2]*(10**-6)
                ERP_data[3] =  ERP_data[3]*(10**-6)
                ERP_data[13] = ERP_data[13]*(10**-6)
                ERP_data[14] = ERP_data[14]*(10**-6)
                del ERP_data[17:]
                ERP_data.append(dt_delivery)

                ERP.append(ERP_data)

#        Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
#                                                 'X-RT','Y-RT','S-XR','S-YR',"Delivery_date"])
#        return Erp_end
#

    if ac in ('wum','grg','esa', 'mit', 'ngs', 'sio'):
        for i in range(tamanho+1):
            linhaatual = linecache.getline(caminho_arq, i)
            if linhaatual[0:1] in numeros:
                ERP_data = linhaatual.split()
                for j in range(len(ERP_data)):
                    ERP_data[j] = float(ERP_data[j])
                ERP_data.insert(0,ac)
                ERP_data[2] =  ERP_data[2]*(10**-6)
                ERP_data[3] =  ERP_data[3]*(10**-6)
                ERP_data[13] = ERP_data[13]*(10**-6)
                ERP_data[14] = ERP_data[14]*(10**-6)
                ERP_data.append(dt_delivery)

                ERP.append(ERP_data)
#        Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
#                                                 'X-RT','Y-RT','S-XR','S-YR'])
#        return Erp_end
#
    if ac in ('gbm', 'gfz'):
        for i in range(tamanho+1):
            linhaatual = linecache.getline(caminho_arq, i)
            if linhaatual[0:1] in numeros:
                ERP_data = linhaatual.split()
                for j in range(len(ERP_data)):
                    ERP_data[j] = float(ERP_data[j])
                ERP_data.insert(0,ac)
                ERP_data[2] =  ERP_data[2]*(10**-6)
                ERP_data[3] =  ERP_data[3]*(10**-6)
                ERP_data[13] = ERP_data[13]*(10**-6)
                ERP_data[14] = ERP_data[14]*(10**-6)
                ERP_data.append(dt_delivery)

                ERP.append(ERP_data)
#        Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
#                                                 'X-RT','Y-RT','S-XR','S-YR'])  ##EH TBM O RATE XY POR DIA??????
#        return Erp_end

    header = []
    if ac in ('emr'):
        for i in range(tamanho+1):
            linhaatual = linecache.getline(caminho_arq, i)
            if linhaatual == 'EOP  SOLUTION':
                del ERP_data[:]
                header = ['EOP  SOLUTION']
            if linhaatual[0:1] in numeros and 'EOP  SOLUTION' in header:
                ERP_data = linhaatual.split()
                for j in range(len(ERP_data)):
                    ERP_data[j] = float(ERP_data[j])
                ERP_data.insert(0,ac)
                ERP_data[2] =  ERP_data[2]*(10**-6)
                ERP_data[3] =  ERP_data[3]*(10**-6)
                ERP_data[13] = ERP_data[13]*(10**-6)
                ERP_data[14] = ERP_data[14]*(10**-6)
                del ERP_data[17:]
                ERP_data.append(dt_delivery)

                ERP.append(ERP_data)

    Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD',
                                         'NR', 'NF', 'NT',
                                         'X-RT','Y-RT','S-XR','S-YR',
                                         'Delivered_date'])
    return Erp_end

########################################################################################################################################
############################################### READ IERS GUS

def read_IERS(file_path_in):
#    file_path_in =  '/dsk/mansur/MGEX_C/ERP_IERS'

    file = open(file_path_in,'r')
    fil = file.readlines()
    file.close()
    return fil
def read_IERS_info(fil, mjd):
    #    mjd = int(mjd)
    X, Y, ut1utc, dx, dy, xerr, yerr, ut1utcerr, dxerr, dyerr, LOD, sig_lod, OMEGA, sig_ome = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
    PoleXY = []

#    mjd = 58230

    for i in range(len(fil)):
    #    linhaatual = linecache.getline(file_path_in, i)
        linhaatual = fil[i]
        if (linhaatual[15:20]) == str(mjd) and (linhaatual[104:105]) != '':
          X = float(linhaatual[23:30])*(10**-3)
          Y = float(linhaatual[31:39])*(10**-3)
          ut1utc = float(linhaatual[40:49])*(10**-3)
          dx = float(linhaatual[51:58])*(10**-3)
          dy = float(linhaatual[59:65])*(10**-3)
          xerr = float(linhaatual[68:74])*(10**-3)
          yerr = float(linhaatual[77:82])*(10**-3)
          ut1utcerr = float(linhaatual[86:93])*(10**-3)
          dxerr = float(linhaatual[94:100])*(10**-3)
          dyerr = float(linhaatual[101:106])*(10**-3)
        elif (linhaatual[15:20]) == str(mjd) and (linhaatual[66:68]) != '':
          LOD = float(linhaatual[22:29])*(10**-3)
          sig_lod = float(linhaatual[30:37])*(10**-3)
          OMEGA = float(linhaatual[38:53])*(10**-3)
          sig_ome = float(linhaatual[56:70])*(10**-3)

    line_data = [mjd,X,Y,ut1utc,dx,dy,xerr,yerr,ut1utcerr,dxerr,dyerr,LOD,sig_lod,OMEGA,sig_ome]
    PoleXY.append(line_data)



    PoleXY_data = pd.DataFrame(PoleXY, columns=['mjd','X','Y','ut1-utc','dx','dy','xerr','yerr','ut1utcerr','dxerr','dyerr','LOD','sig_lod','OMEGA','sig_ome'])
    PoleXY_data.set_index('mjd',inplace=True)
    return PoleXY_data




################################################## END READ IERS GUS




##################################### MAKE A LIST OF FILES IN A FOLDER GUS

def list_files(dire,file = None):
    """
    written for MGEX combi by GM
    should be discontinued
    """
    List_arq = os.listdir(dire)

    path = []

    for l in range(len(List_arq)):
        atual = List_arq[l]
        if file == 'sp3':
            if atual[-3:] == 'SP3' or atual[-3:] == 'sp3':
                path.append(os.path.join(dire,List_arq[l]))

        if file == 'clk':
            if atual[-3:] == 'CLK' or atual[-3:] == 'clk':
                path.append(os.path.join(dire,List_arq[l]))

        if file == 'erp':
            if atual[-3:] == 'ERP' or atual[-3:] == 'erp':
                path.append(os.path.join(dire,List_arq[l]))

        if file == 'eph':
            if atual[-3:] == 'EPH' or atual[-3:] == 'eph':
                path.append(os.path.join(dire,List_arq[l]))


        if file == 'snx':
            if atual[-3:] == 'SNX' or atual[-3:] == 'snx':
                path.append(os.path.join(dire,List_arq[l]))


    return path

################################################ END MAKE A LIST OF FILES IN A FOLDER GUS




def diff_pandas(DF,col_name):
    """
    Differentiate a Pandas DataFrame, if index is time

    Parameters
    ----------
    DF : Pandas DataFrame
         input DataFrame

    col_name : str
        the column of the DataFrame you want to differentiate

    Returns
    -------
    DSout : Pandas DataFrame
        Differenciated column of the input DataFrame

    """
    DSout = DF[col_name].diff() / DF[col_name].index.to_series().diff().dt.total_seconds()
    return DSout


def compar_orbit(Data_inp_1,Data_inp_2,step_data = 900,
                 sats_used_list = ['G'],
                 name1='',name2='',use_name_1_2_for_table_name = False,
                 RTNoutput = True,convert_ECEF_ECI=True,
                 clean_null_values = True):
    """
    Compares 2 GNSS orbits files (SP3), and gives a summary plot and a
    statistics table

    Parameters
    ----------
    Data_inp_1 & Data_inp_2 : str or Pandas DataFrame
        contains the orbits or path (string) to the sp3

    step_data : int
        per default data sampling

    sats_used_list : list of str
        used constellation or satellite : G E R C ... E01 , G02 ...
        Individuals satellites are prioritary on whole constellations
        e.g. ['G',"E04"]


    RTNoutput : bool
        select output, Radial Transverse Normal or XYZ

    convert_ECEF_ECI : bool
        convert sp3 ECEF => ECI, must be True in operational !

    name1 & name2 : str (optionals)
        optional custom names for the 2 orbits

    use_name_1_2_for_table_name : bool
        False : use name 1 and 2 for table name, use datafile instead

    clean_null_values : bool or str
        if True or "all" remove sat position in all X,Y,Z values
        are null (0.000000)
        if "any", remove sat position if X or Y or Z is null
        if False, keep everything

    Returns
    -------
    Diff_sat_all : Pandas DataFrame
    contains differences b/w Data_inp_1 & Data_inp_2
    in Radial Transverse Normal OR XYZ frame

        Attributes of Diff_sat_all :
            Diff_sat_all.name : title of the table

    Note
    ----
    clean_null_values if useful (and necessary) only if
    convert_ECEF_ECI = False
    if convert_ECEF_ECI = True, the cleaning will be done by
    a side effect trick : the convertion ECEF => ECI will generate NaN
    for a zero-valued position
    But, nevertheless, activating  clean_null_values = True is better
    This Note is in fact usefull if you want to see bad positions on a plot
    => Then convert_ECEF_ECI = False and clean_null_values = False

    Source
    ------
    "Coordinate Systems", ASEN 3200 1/24/06 George H. Born

    """

    # selection of both used Constellations AND satellites
    const_used_list = []
    sv_used_list    = []
    for sat in sats_used_list:
        if len(sat) == 1:
            const_used_list.append(sat)
        elif len(sat) == 3:
            sv_used_list.append(sat)
            if not sat[0] in const_used_list:
                const_used_list.append(sat[0])

    # Read the files or DataFrames
    # metadata attributes are not copied
    # Thus, manual copy ...
    # (Dirty way, should be impoved without so many lines ...)
    if type(Data_inp_1) is str:
        D1orig = read_sp3(Data_inp_1,epoch_as_pd_index=True)
    else:
        D1orig = Data_inp_1.copy(True)
        try:
            D1orig.name = Data_inp_1.name
        except:
            D1orig.name = "no_name"
        try:
            D1orig.path = Data_inp_1.path
        except:
            D1orig.path = "no_path"
        try:
            D1orig.filename = Data_inp_1.filename
        except:
            D1orig.filename = "no_filename"

    if type(Data_inp_2) is str:
        D2orig = read_sp3(Data_inp_2,epoch_as_pd_index=True)
    else:
        D2orig = Data_inp_2.copy(True)
        try:
            D2orig.name = Data_inp_2.name
        except:
            D2orig.name = "no_name"
        try:
            D2orig.path = Data_inp_2.path
        except:
            D2orig.path = "no_path"
        try:
            D2orig.filename = Data_inp_2.filename
        except:
            D2orig.filename = "no_filename"

    #### NB : It has been decided with GM that the index of a SP3 dataframe
    ####      will be integers, not epoch datetime anymore
    ####      BUT here, for legacy reasons, the index has to be datetime

    if isinstance(D1orig.index[0], (int, np.integer)):
        D1orig.set_index("epoch",inplace=True)

    if isinstance(D2orig.index[0], (int, np.integer)):
        D2orig.set_index("epoch",inplace=True)

    Diff_sat_stk = []

    # This block is for removing null values
    if clean_null_values:
        if clean_null_values == "all":
            all_or_any = np.all
        elif clean_null_values == "any":
            all_or_any = np.any
        else:
            all_or_any = np.all

        xyz_lst = ['x','y','z']

        D1_null_bool = all_or_any(np.isclose(D1orig[xyz_lst],0.),axis=1)
        D2_null_bool = all_or_any(np.isclose(D2orig[xyz_lst],0.),axis=1)

        D1 = D1orig[np.logical_not(D1_null_bool)]
        D2 = D2orig[np.logical_not(D2_null_bool)]

        if np.any(D1_null_bool) or np.any(D2_null_bool):
            print("WARN : Null values contained in SP3 files : ")
            print("f1:" , np.sum(D1_null_bool) , genefun.join_improved(" " ,
                  *list(set(D1orig[D1_null_bool]["sat"]))))
            print("f2:" , np.sum(D2_null_bool) , genefun.join_improved(" " ,
                  *list(set(D2orig[D2_null_bool]["sat"]))))

    else:
        D1 = D1orig.copy()
        D2 = D2orig.copy()

    for constuse in const_used_list:
        D1const = D1[D1['const'] == constuse]
        D2const = D2[D2['const'] == constuse]

        # checking if the data correspond to the step
        bool_step1 = np.mod((D1const.index - np.min(D1.index)).seconds,step_data) == 0
        bool_step2 = np.mod((D2const.index - np.min(D2.index)).seconds,step_data) == 0

        D1window = D1const[bool_step1]
        D2window = D2const[bool_step2]

        # find common sats and common epochs
        sv_set   = sorted(list(set(D1window['sv']).intersection(set(D2window['sv']))))
        epoc_set = sorted(list(set(D1window.index).intersection(set(D2window.index))))

        # if special selection of sats, then apply it
        # (it is late and this selection is incredibely complicated ...)
        if np.any([True  if constuse in e else False for e in sv_used_list]):
            # first find the selected sats for the good constellation
            sv_used_select_list = [int(e[1:]) for e in sv_used_list if constuse in e]
            #and apply it
            sv_set = sorted(list(set(sv_set).intersection(set(sv_used_select_list))))

        for svv in sv_set:
            # First research : find corresponding epoch for the SV
            # this one is sufficent if there is no gaps (e.g. with 0.00000) i.e.
            # same nb of obs in the 2 files
            # NB : .reindex() is smart, it fills the DataFrame
            # with NaN
            try:
                D1sv_orig = D1window[D1window['sv'] == svv].reindex(epoc_set)
                D2sv_orig = D2window[D2window['sv'] == svv].reindex(epoc_set)
            except Exception as exce:
                print("ERR : Unable to re-index with an unique epoch")
                print("      are you sure there is no multiple-defined epochs for the same sat ?")
                print("      it happens e.g. when multiple ACs are in the same DataFrame ")
                print("TIP : Filter the input Dataframe before calling this fct with")
                print("      DF = DF[DF['AC'] == 'gbm']")
                raise exce

            # Second research, it is a security in case of gap
            # This step is useless, because .reindex() will fill the DataFrame
            # with NaN
            if len(D1sv_orig) != len(D2sv_orig):
                print("INFO : different epochs nbr for SV",svv,len(D1sv_orig),len(D2sv_orig))
                epoc_sv_set = sorted(list(set(D1sv_orig.index).intersection(set(D2sv_orig.index))))
                D1sv = D1sv_orig.loc[epoc_sv_set]
                D2sv = D2sv_orig.loc[epoc_sv_set]
            else:
                D1sv = D1sv_orig
                D2sv = D2sv_orig

            P1     = D1sv[['x','y','z']]
            P2     = D2sv[['x','y','z']]

            # Start ECEF => ECI
            if convert_ECEF_ECI:
                # Backup because the columns xyz will be reaffected
                D1sv_bkp = D1sv.copy()
                D2sv_bkp = D2sv.copy()

                P1b = geok.ECEF2ECI(np.array(P1),geok.dt_gpstime2dt_utc(P1.index.to_pydatetime(),out_array=True))
                P2b = geok.ECEF2ECI(np.array(P2),geok.dt_gpstime2dt_utc(P2.index.to_pydatetime(),out_array=True))

                D1sv[['x','y','z']] = P1b
                D2sv[['x','y','z']] = P2b

                P1  = D1sv[['x','y','z']]
                P2  = D2sv[['x','y','z']]
            # End ECEF => ECI

            if not RTNoutput:
                # Compatible with the documentation +
                # empirically tested with OV software
                # it is  P1 - P2 (and not P2 - P1)
                Delta_P = P1 - P2


                Diff_sat = Delta_P.copy()
                Diff_sat.columns = ['dx','dy','dz']

            else:
                rnorm = np.linalg.norm(P1,axis=1)

                Vx = diff_pandas(D1sv,'x')
                Vy = diff_pandas(D1sv,'y')
                Vz = diff_pandas(D1sv,'z')

                V =  pd.concat((Vx , Vy , Vz),axis=1)
                V.columns = ['vx','vy','vz']

                R = P1.divide(rnorm,axis=0)
                R.columns = ['xnorm','ynorm','znorm']

                H      = pd.DataFrame(np.cross(R,V),columns=['hx','hy','hz'])
                hnorm  = np.linalg.norm(H,axis=1)

                C         = H.divide(hnorm,axis=0)
                C.columns = ['hxnorm','hynorm','hznorm']

                I         = pd.DataFrame(np.cross(C,R),columns=['ix','iy','iz'])

                R_ar = np.array(R)
                I_ar = np.array(I)
                C_ar = np.array(C)

                R_ar[1]
                Beta = np.stack((R_ar,I_ar,C_ar),axis=1)

                # Compatible with the documentation +
                # empirically tested with OV software
                # it is  P1 - P2 (and not P2 - P1)
                Delta_P = P1 - P2

                # Final determination
                Astk = []

                for i in range(len(Delta_P)):
                    A = np.dot(Beta[i,:,:],np.array(Delta_P)[i])
                    Astk.append(A)

                Diff_sat = pd.DataFrame(np.vstack(Astk),
                                   index = P1.index,columns=['dr','dt','dn'])

            Diff_sat = Diff_sat * 1000 # metrer conversion

            Diff_sat['const'] = [constuse] * len(Diff_sat.index)
            Diff_sat['sv']    = [svv]      * len(Diff_sat.index)
            Diff_sat['sat']   = [constuse + str(svv).zfill(2)] * len(Diff_sat.index)

            Diff_sat_stk.append(Diff_sat)


    Diff_sat_all = pd.concat(Diff_sat_stk)
    Date = Diff_sat.index[0]

    # Attribute definition
    if RTNoutput:
        Diff_sat_all.frame_type = 'RTN'

        # Pandas donesn't manage well iterable as attribute
        # So, it is separated
        Diff_sat_all.frame_col_name1 = 'dr'
        Diff_sat_all.frame_col_name2 = 'dt'
        Diff_sat_all.frame_col_name3 = 'dn'

    else:
        # Pandas donesn't manage well iterable as attribute
        # So, it is separated
        Diff_sat_all.frame_col_name1 = 'dx'
        Diff_sat_all.frame_col_name2 = 'dy'
        Diff_sat_all.frame_col_name3 = 'dz'

        if convert_ECEF_ECI:
            Diff_sat_all.frame_type = 'ECI'
        else:
            Diff_sat_all.frame_type = 'ECEF'


    # Name definitions
    if name1:
        Diff_sat_all.name1 = name1
    else:
        Diff_sat_all.name1 = D1orig.name

    if name2:
        Diff_sat_all.name2 = name2
    else:
        Diff_sat_all.name2 = D2orig.name

    Diff_sat_all.filename1 = D1orig.filename
    Diff_sat_all.filename2 = D2orig.filename

    Diff_sat_all.path1 = D1orig.path
    Diff_sat_all.path2 = D2orig.path

    Diff_sat_all.name = ' '.join(('Orbits comparison ('+Diff_sat_all.frame_type +') b/w',
                                  Diff_sat_all.name1 ,'(ref.) and',
                                  Diff_sat_all.name2 ,',',Date.strftime("%Y-%m-%d"),
                                  ', doy', str(geok.dt2doy(Date))))


    return Diff_sat_all


def compar_orbit_plot(Diff_sat_all_df_in,
                      save_plot=False,
                      save_plot_dir="",
                      save_plot_name="auto",
                      save_plot_ext=(".pdf",".png",".svg")):
    """
    General description

    Parameters
    ----------
    Diff_sat_all_df_in : DataFrame
        a DataFrame produced by compar_orbit

    Returns
    -------
    None if no save is asked
    export path (str) if save is asked
    but plot a plot anyway
    """

    import matplotlib.dates as mdates
    fig,[axr,axt,axn] = plt.subplots(3,1,sharex='all')

    satdispo = natsorted(list(set(Diff_sat_all_df_in['sat'])))

    SymbStk = []

    cm = plt.get_cmap('viridis')
    NUM_COLORS = len(satdispo)
    Colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]

    # Pandas donesn't manage well iterable as attribute
    # So, it is separated
    col_name0 = Diff_sat_all_df_in.frame_col_name1
    col_name1 = Diff_sat_all_df_in.frame_col_name2
    col_name2 = Diff_sat_all_df_in.frame_col_name3

    for satuse,color in zip(satdispo,Colors):
        Diffuse = Diff_sat_all_df_in[Diff_sat_all_df_in['sat'] == satuse]

        Time = Diffuse.index
        R    = Diffuse[col_name0]
        T    = Diffuse[col_name1]
        N    = Diffuse[col_name2]

        #fig.fmt_xdata = mdates.DateFormatter('%Y-%m-%d')

        Symb = axr.plot(Time,R,label=satuse,c=color)
        axt.plot(Time,T,label=satuse,c=color)
        axn.plot(Time,N,label=satuse,c=color)

        SymbStk.append(Symb[0])

        fig.autofmt_xdate()

    if Diff_sat_all_df_in.frame_type == 'RTN':
        axr.set_ylabel('Radial diff. (m)')
        axt.set_ylabel('Transverse diff. (m)')
        axn.set_ylabel('Normal diff. (m)')

    else:
        axr.set_ylabel(Diff_sat_all_df_in.frame_type + ' X diff. (m)')
        axt.set_ylabel(Diff_sat_all_df_in.frame_type + ' Y diff. (m)')
        axn.set_ylabel(Diff_sat_all_df_in.frame_type + ' Z diff. (m)')


    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    axr.yaxis.set_major_formatter(y_formatter)
    axt.yaxis.set_major_formatter(y_formatter)
    axn.yaxis.set_major_formatter(y_formatter)

    import matplotlib.dates as mdates
    fig.fmt_xdata = mdates.DateFormatter('%Y-%m-%d')

    lgd = fig.legend(tuple(SymbStk), satdispo , loc='lower center',ncol=8,
                     columnspacing=1)

    fig.set_size_inches(8.27,11.69)
    plt.suptitle(Diff_sat_all_df_in.name)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(bottom=0.15)

    if save_plot:
        if save_plot_name == "auto":
            save_plot_name = "_".join((Diff_sat_all_df_in.name1,
                                      Diff_sat_all_df_in.name2,
                                      Diff_sat_all_df_in.index.min().strftime("%Y-%m-%d")))

        for ext in save_plot_ext:
            save_plot_path = os.path.join(save_plot_dir,save_plot_name)
            plt.savefig(save_plot_path + ext)
    else:
        save_plot_path = None

    return save_plot_path

def compar_orbit_table(Diff_sat_all_df_in,GRGS_style = True,
                       light_tab  = False):
    """
    Generate a table with statistical indicators for an orbit comparison
    (RMS mean, standard dev, ...)
    Parameters
    ----------
    Diff_sat_all_df_in : Pandas DataFrame
        a DataFrame produced by compar_orbit

    GRGS_style : bool
        RMS calc based on the GRGS definition of the RMS (OV help)

    light_tab : bool
        produce a table with only RMS, with min/max/arithmetic instead

    Returns
    -------
    Compar_tab_out : DataFrame
        Statistical results of the comparison

    Note
    ----
    you can pretty print the output DataFrame using tabular module
    here is a template:

    >>> from tabulate import tabulate
    >>> print(tabulate(ComparTable,headers="keys",floatfmt=".4f"))
    """

    sat_list = genefun.uniq_and_sort(Diff_sat_all_df_in['sat'])

    # Pandas donesn't manage well iterable as attribute
    # So, it is separated
    col_name0 = Diff_sat_all_df_in.frame_col_name1
    col_name1 = Diff_sat_all_df_in.frame_col_name2
    col_name2 = Diff_sat_all_df_in.frame_col_name3

    rms_stk = []

    for sat in sat_list:
        Diffwork = genefun.df_sel_val_in_col(Diff_sat_all_df_in,'sat',sat)

        if not GRGS_style:
            rms_A = geok.rms_mean(Diffwork[col_name0])
            rms_B = geok.rms_mean(Diffwork[col_name1])
            rms_C = geok.rms_mean(Diffwork[col_name2])
        else:
            rms_A = geok.rms_mean(Diffwork[col_name0] - Diffwork[col_name0].mean())
            rms_B = geok.rms_mean(Diffwork[col_name1] - Diffwork[col_name1].mean())
            rms_C = geok.rms_mean(Diffwork[col_name2] - Diffwork[col_name2].mean())

        RMS3D = np.sqrt(rms_A**2 + rms_B**2 + rms_C**2)

        min_A = Diffwork[col_name0].min()
        min_B = Diffwork[col_name1].min()
        min_C = Diffwork[col_name2].min()

        max_A = Diffwork[col_name0].max()
        max_B = Diffwork[col_name1].max()
        max_C = Diffwork[col_name2].max()

        mean_A = Diffwork[col_name0].mean()
        mean_B = Diffwork[col_name1].mean()
        mean_C = Diffwork[col_name2].mean()

        if light_tab:
            rms_stk.append([rms_A,rms_B,rms_C])
        else:
            rms_stk.append([rms_A,rms_B,rms_C,RMS3D,
                            min_A,max_A,mean_A,
                            min_B,max_B,mean_B,
                            min_C,max_C,mean_C])


    #################################
             # ALL SATS

    rms_A = geok.rms_mean(Diff_sat_all_df_in[col_name0] - Diff_sat_all_df_in[col_name0].mean())
    rms_B = geok.rms_mean(Diff_sat_all_df_in[col_name1] - Diff_sat_all_df_in[col_name1].mean())
    rms_C = geok.rms_mean(Diff_sat_all_df_in[col_name2] - Diff_sat_all_df_in[col_name2].mean())

    RMS3D = np.sqrt(rms_A**2 + rms_B**2 + rms_C**2)

    min_A = Diff_sat_all_df_in[col_name0].min()
    min_B = Diff_sat_all_df_in[col_name1].min()
    min_C = Diff_sat_all_df_in[col_name2].min()

    max_A = Diff_sat_all_df_in[col_name0].max()
    max_B = Diff_sat_all_df_in[col_name1].max()
    max_C = Diff_sat_all_df_in[col_name2].max()

    mean_A = Diff_sat_all_df_in[col_name0].mean()
    mean_B = Diff_sat_all_df_in[col_name1].mean()
    mean_C = Diff_sat_all_df_in[col_name2].mean()

    if light_tab:
        rms_stk.append([rms_A,rms_B,rms_C])
    else:
        rms_stk.append([rms_A,rms_B,rms_C,RMS3D,
                        min_A,max_A,mean_A,
                        min_B,max_B,mean_B,
                        min_C,max_C,mean_C])

             # ALL SATS
    #################################

    if  Diff_sat_all_df_in.frame_type == 'RTN':
        if light_tab:
            cols_nam = ["rmsR","rmsT","rmsN","rms3D"]
        else:
            cols_nam = ["rmsR","rmsT","rmsN","rms3D",
                        "minR","maxR","meanR",
                        "minT","maxT","meanT",
                        "minN","maxN","meanN"]

    else:
        if light_tab:
            cols_nam = ["rmsX","rmsY","rmsZ","rms3D"]
        else:
            cols_nam = ["rmsX","rmsY","rmsZ","rms3D",
                        "minX","maxX","meanX",
                        "minY","maxY","meanY",
                        "minZ","maxZ","meanZ"]

    Compar_tab_out     = pd.DataFrame(rms_stk,index=sat_list + ['ALL'],
                                      columns=cols_nam)

    return Compar_tab_out


def compar_orbit_frontend(DataDF1,DataDF2,ac1,ac2):
    K = compar_orbit(DataDF1[DataDF1["AC"] == ac1],DataDF2[DataDF2["AC"] == ac2])
    compar_orbit_plot(K)
    return K


def compar_sinex(snx1 , snx2 , stat_select = None, invert_select=False,
                 out_means_summary=True,out_meta=True,out_dataframe = True,
                 manu_wwwwd=None):

    if type(snx1) is str:
        week1 = genefun.split_improved(os.path.basename(snx1),"_",".")[:]
        week2 = genefun.split_improved(os.path.basename(snx2),"_",".")[:]
        if week1 != week2:
            print("WARN : Dates of 2 input files are differents !!! It might be very bad !!!",week1,week2)
        else:
            wwwwd = week1
        D1 = gfc.read_sinex(snx1,True)
        D2 = gfc.read_sinex(snx2,True)
    else:
        print("WARN : you are giving the SINEX input as a DataFrame, wwwwd has to be given manually using manu_wwwwd")
        D1 = snx1
        D2 = snx2


    if manu_wwwwd:
        wwwwd = manu_wwwwd


    STATCommon  = set(D1["STAT"]).intersection(set(D2["STAT"]))

    if stat_select:

        STATCommon_init = list(STATCommon)

        if invert_select:
            select_fct = lambda x : not x
        else:
            select_fct = lambda x : x

        if type(stat_select) is str:
            STATCommon = [sta for sta in STATCommon_init if select_fct(re.search(stat_select, sta)) ]
        elif genefun.is_iterable(stat_select):
            STATCommon = [sta for sta in STATCommon_init if select_fct(sta in stat_select) ]
        else:
            print("WARN : check type of stat_select")

    D1Common = D1[D1["STAT"].isin(STATCommon)].sort_values("STAT").reset_index(drop=True)
    D2Common = D2[D2["STAT"].isin(STATCommon)].sort_values("STAT").reset_index(drop=True)


    Ddiff = pd.DataFrame()
    Ddiff = Ddiff.assign(STAT=D1Common["STAT"])

    #### XYZ Part
    for xyz in ("x","y","z"):

        dif = pd.to_numeric((D2Common[xyz] - D1Common[xyz]))

        Ddiff = Ddiff.assign(xyz=dif)
        Ddiff = Ddiff.rename(columns={"xyz": xyz})

    D3D = np.sqrt((Ddiff["x"]**2 + Ddiff["y"]**2 + Ddiff["z"]**2 ).astype('float64'))

    Ddiff = Ddiff.assign(d3D_xyz=D3D)

    ### ENU Part
    E , N ,U = [] , [] , []
    enu_stk = []

    for (_,l1) , (_,l2) in zip( D1Common.iterrows() , D2Common.iterrows() ):
        enu   = geok.XYZ2ENU_2(l1["x"],l1["y"],l1["z"],l2["x"],l2["y"],l2["z"])
        enu_stk.append(np.array(enu))


    if len(enu_stk) == 0:
        E,N,U = np.array([]) , np.array([]) , np.array([])
    else:
        ENU = np.hstack(enu_stk)
        E,N,U = ENU[0,:] , ENU[1,:] , ENU[2,:]


    D2D = np.sqrt((E**2 + N**2).astype('float64'))
    D3D = np.sqrt((E**2 + N**2 + U**2 ).astype('float64'))

    Ddiff = Ddiff.assign(e=E)
    Ddiff = Ddiff.assign(n=N)
    Ddiff = Ddiff.assign(u=U)
    Ddiff = Ddiff.assign(d2D_enu=D2D)
    Ddiff = Ddiff.assign(d3D_enu=D3D)

    #    E,N,U    = geok.XYZ2ENU_2((X,Y,Z,x0,y0,z0))
    #    E,N,U    = geok.XYZ2ENU_2((X,Y,Z,x0,y0,z0))

    if out_dataframe:
        out_meta = True


    if not out_means_summary:
        print("INFO : this is not used operationally and it can be improved")
        return Ddiff
    else:
        output = []

        col_names = ("x","y","z","d3D_xyz",
                     "e","n","u","d2D_enu","d3D_enu")

        for xyz in col_names:
            output.append(geok.rms_mean(Ddiff[xyz]))
        for xyz in col_names:
            output.append(np.nanmean(Ddiff[xyz]))
        for xyz in col_names:
            output.append(np.nanstd(Ddiff[xyz]))

        if out_meta:
            print(wwwwd)
            nstat = len(STATCommon)
            week   = int(wwwwd[:4])
            day    = int(wwwwd[4:])
            output = [week , day ,nstat] + output


        if not out_dataframe:
            return tuple(output)
        else:

            output_DF = pd.DataFrame(output).transpose()

            output_DF.columns = ["week","dow","nbstat",
             "x_rms","y_rms","z_rms","d3D_xyz_rms",
             "e_rms","n_rms","u_rms","d2D_enu_rms","d3D_enu_rms",
             "x_ari","y_ari","z_ari","d3D_xyz_ari",
             "e_ari","n_ari","u_ari","d2D_enu_ari","d3D_enu_ari",
             "x_ari","y_std","z_std","d3D_xyz_std",
             "e_ari","n_std","u_std","d2D_enu_std","d3D_enu_std"]

            return output_DF



def read_clk(file_path_in):
    """
    This function is discontinued, use read_clk1 instead
    """

    i  = 0
    with open(file_path_in) as oF:
        for l in oF:
            i += 1
            if 'END OF HEADER'  in l:
                break

    colnam = ['type' , 'name' , 'year' ,'month' ,'day' , 'h' ,'m' ,'s','nb_val' ,
              'clk_bias','clk_bias_std','clk_rate','clk_rate_std','clk_accel','clk_accel_std']
    DF = pd.read_table(file_path_in,skiprows=i,header = -1,names = colnam,delim_whitespace=True)

    DF['epoc'] =  pd.to_datetime(DF[[ 'year' ,'month' ,'day','h','m','s']])

    # Removing NaN columns
    for colnam in ['clk_bias','clk_bias_std','clk_rate',
                   'clk_rate_std','clk_accel','clk_accel_std']:
        if np.all(np.isnan(DF[colnam])):
            DF = DF.drop(colnam,axis=1)

    return DF

def clk_diff(file_1,file_2):
    if type(file_1) is str:
        DF1 = read_clk(file_1)
    else:
        DF1 = file_1

    if type(file_2) is str:
        DF2 = read_clk(file_2)
    else:
        DF2 = file_2

    name_list = sorted(set(DF1['name']).intersection(set(DF2['name'])))

    epoc_common_stk = []
    name_list_stk   = []
    type_list_stk   = []
    dClk_stk        = []

    for nam in name_list:
        DF1nam_bool = (DF1['name'] == nam)
        DF2nam_bool = (DF2['name'] == nam)

        DF1nam = DF1[DF1nam_bool]
        DF2nam = DF2[DF2nam_bool]

        epoc_common = set(DF1nam['epoc']).intersection(set(DF2nam['epoc']))

        DF1a = DF1nam[ DF1nam['epoc'].isin(epoc_common) ]
        DF2a = DF2nam[ DF2nam['epoc'].isin(epoc_common) ]

        dClk = np.array(DF1a['clk_bias']) - np.array(DF2a['clk_bias'])

        type_list_stk   += list(DF1a['type'])
        name_list_stk   += [nam] * len(dClk)
        epoc_common_stk += list(epoc_common)
        dClk_stk        += list(dClk)

    DFdiff = pd.DataFrame((list(zip(name_list_stk,type_list_stk,
                                   epoc_common_stk,dClk_stk))),
                              columns = ['name','type','date','diff'])

    return DFdiff


### compar_plot
def stations_in_EPOS_sta_coords_file_mono(coords_file_path):
    """
    Gives stations in a EPOS coords. file
    like : YYYY_DDD_sta_coordinates

    Args :
        coords_file_path : path of the EPOS coords. file
    Returns :
        epoch : the main mean epoch in the EPOS coords. file
        stats_list : list of 4 char station list
    """

    SITE_line_list = genefun.grep(coords_file_path , " SITE            m")

    stats_list = []
    mean_mjd_list = []
    for l in SITE_line_list:
        stat = l.split()[8].lower()
        stats_list.append(stat)
        mean_mjd = np.mean([float(l.split()[6]) , float(l.split()[7])])
        mean_mjd_list.append(mean_mjd)

    mjd_final = genefun.most_common(mean_mjd_list)
    epoch = geok.MJD2dt(mjd_final)

    return epoch , stats_list


def stations_in_sinex_mono(sinex_path):
    """
    Gives stations list in a SINEX
    Args :
        sinex_path : path of the SINEX file
    Returns :
        epoch : the main mean epoch in the SINEX
        stats_list : list of 4 char station list

    """
    extract = genefun.extract_text_between_elements(sinex_path,'+SITE/ID','-SITE/ID')
    extract = extract.split('\n')
    extract2 = []
    for e in extract:
        if e != '' and e[0] == ' ' and e != '\n':
            extract2.append(e)

    stats_list = [e.split()[0].lower() for e in extract2]

    extract = genefun.extract_text_between_elements(sinex_path,'+SOLUTION/EPOCHS','-SOLUTION/EPOCHS')
    extract = extract.split('\n')
    extract2 = []
    for e in extract:
        if e != '' and e[0] == ' ' and e != '\n':
            extract2.append(e)

    epoch = geok.datestr_sinex_2_dt(genefun.most_common([e.split()[-1] for e in extract2]))

    return epoch , stats_list


####################### Compare troposphere delay  ####################################################
def compare_trop_ties(input_file,STA1,STA2,coord_file="",grid_met="",apply_ties=False,mode="DF",mode_coor="sinex",coord_t="static"):
    """
    Calculate differences of tropospheric delay and gradients between selected stations (Atmospheric ties)
    Args  :
          sinex_file : troposphere solutions from sinex
          STA1 : Reference station
          STA2 : Rover station
          coord_file : Station coordinates file path
          grid_met : Grid file for meteological information from Global Pressture Temperature (GPT)
          apply_ties : Apply height coorections from standard ties (Default No)
          mode : data in DataFrame (DF) or SINEX (SINEX) (Default DataFrame)
          mode_coor : Coordinates file format in SINEX or EPOS or DataFrame format (Argruments : sinex , epos , df)
          coord_t : static (default) or kinematic
    Return :
         trop_diff : difference of tropospheric delay and gradients between selected stations (Atmospheric ties)
                     and uncertainty of atmospheric ties and gradients ties
    """

    if mode == "SINEX":
        trop_pd = gfc.read_snx_trop(input_file)
    elif mode == "DF":
        trop_pd = input_file
    else:
        import sys
        print("No this option for troposphere solution")
        sys.exit()

    trop_ref = trop_pd[trop_pd.STAT == STA1]
    trop_rov = trop_pd[trop_pd.STAT == STA2]

    if trop_ref.empty:
        print("No solution for "+STA1+" Reference station")
        return None,None

    if trop_rov.empty:
        print("No solution for "+STA2+" Rover station")
        return None,None

    diff_pd = pd.merge(trop_ref,trop_rov,how='outer',on='epoc')
    diff_pd = diff_pd.dropna()
    # Tropospheric ties
    diff_pd['Trop_ties'] = diff_pd['tro_x'] - diff_pd['tro_y']
    diff_pd['STrop_ties'] = np.nan # add blank column before input values
    diff_pd['STrop_ties'] = diff_pd.apply(lambda x: np.round(np.sqrt(x.stro_x**2 + x.stro_y**2),2),axis=1)

    # North gradients ties
    diff_pd['Ngra_ties'] = diff_pd['tgn_x'] - diff_pd['tgn_y']
    diff_pd['SNgra_ties'] = np.nan # add blank column before input values
    diff_pd['SNgra_ties'] = diff_pd.apply(lambda x: np.round(np.sqrt(x.stgn_x**2 + x.stgn_y**2),2),axis=1)

    # East gradients ties
    diff_pd['Egra_ties'] = diff_pd['tge_x'] - diff_pd['tge_y']
    diff_pd['SEgra_ties'] = np.nan # add blank column before input values
    diff_pd['SEgra_ties'] = diff_pd.apply(lambda x: np.round(np.sqrt(x.stge_x**2 + x.stge_y**2),2),axis=1)

     # Apply height corrections from Standard ties
    if apply_ties == True:

        if isinstance(grid_met,str):
            import sys
            print("Please read grid file before use this function")
            sys.exit()

        if mode_coor == "epos" and coord_t == "static":
            coord = read_epos_sta_coords_mono(coord_file)
            # Extract coordinates Ref station in Lat. Lon. height
            coord_ref = coord[coord.site==STA1]
            lat_ref , lon_ref , h_ref = geok.XYZ2GEO(coord_ref.x,coord_ref.y,coord_ref.z,False)
            # Extract coordinates Rov station in Lat. Lon. height
            coord_rov = coord[coord.site==STA2]
            lat_rov , lon_rov , h_rov = geok.XYZ2GEO(coord_rov.x,coord_rov.y,coord_rov.z,False)

            #Merge coordinates results
            coord_res = pd.merge(coord_ref,coord_rov,how='outer',on='MJD_ref')
            coord_res['lat_ref'] , coord_res['lon_ref'] , coord_res['h_ref'] = float(lat_ref) , float(lon_ref) , float(h_ref)
            coord_res['lat_rov'] , coord_res['lon_rov'] , coord_res['h_rov'] = float(lat_rov) , float(lon_rov) , float(h_rov)

            #drop unnecessary column
            coord_res = coord_res.drop([ 'site_num_x', 'tecto_plate_x', 'MJD_start_x',
                                        'MJD_end_x', 'Vx_x',
                                        'Vy_x', 'Vz_x', 'sVx_x', 'sVy_x', 'sVz_x', 'site_num_y',
                                        'tecto_plate_y', 'MJD_start_y', 'MJD_end_y', 'Vx_y', 'Vy_y', 'Vz_y', 'sVx_y', 'sVy_y',
                                        'sVz_y'],axis=1)

        elif mode_coor == "sinex" and coord_t == "static":
            coord = gfc.read_sinex(coord_file,True)
             # Extract coordinates Ref station in Lat. Lon. height
            coord_ref = coord[coord.STAT==STA1]
            lat_ref , lon_ref , h_ref = geok.XYZ2GEO(coord_ref.x,coord_ref.y,coord_ref.z,False)
            # Extract coordinates Rov station in Lat. Lon. height
            coord_rov = coord[coord.STAT==STA2]
            lat_rov , lon_rov , h_rov = geok.XYZ2GEO(coord_rov.x,coord_rov.y,coord_rov.z,False)

            #Merge coordinates results
            coord_res = pd.merge(coord_ref,coord_rov,how='outer',on='epoc')
            coord_res['lat_ref'] , coord_res['lon_ref'] , coord_res['h_ref'] = float(lat_ref) , float(lon_ref) , float(h_ref)
            coord_res['lat_rov'] , coord_res['lon_rov'] , coord_res['h_rov'] = float(lat_rov) , float(lon_rov) , float(h_rov)

             #drop unnecessary columns
            coord_res = coord_res.drop(['AC_x', 'soln_x', 'vx_x', 'vy_x', 'vz_x', 'svx_x', 'svy_x', 'svz_x', 'start_x',
                                       'end_x', 'AC_y', 'soln_y', 'vx_y', 'vy_y', 'vz_y', 'svx_y', 'svy_y', 'svz_y',
                                       'start_y', 'end_y'],axis=1)

        elif mode_coor == "epos" and coord_t == "kinematic":
            coord = read_epos_sta_kinematics(coord_file)
            # Extract coordinates Ref station in Lat. Lon. height
            coord_ref = coord[coord.site==STA1.lower()]
            f_x_ref , f_y_ref , f_z_ref = interpolator_with_extrapolated(coord_ref.MJD_epo,coord_ref.x,coord_ref.y,coord_ref.z)
            x_ref_new = f_x_ref(geok.dt2MJD(diff_pd.epoc))
            y_ref_new = f_y_ref(geok.dt2MJD(diff_pd.epoc))
            z_ref_new = f_z_ref(geok.dt2MJD(diff_pd.epoc))
            lat_ref , lon_ref , h_ref = geok.XYZ2GEO(x_ref_new,y_ref_new,z_ref_new,False)

            # Extract coordinates Rov station in Lat. Lon. height
            coord_rov = coord[coord.site==STA2.lower()]
            f_x_rov , f_y_rov , f_z_rov = interpolator_with_extrapolated(coord_rov.MJD_epo,coord_rov.x,coord_rov.y,coord_rov.z)
            x_rov_new = f_x_rov(geok.dt2MJD(diff_pd.epoc))
            y_rov_new = f_y_rov(geok.dt2MJD(diff_pd.epoc))
            z_rov_new = f_z_rov(geok.dt2MJD(diff_pd.epoc))
            lat_rov , lon_rov , h_rov = geok.XYZ2GEO(x_rov_new,y_rov_new,z_rov_new,False)

            diff_pd['lat_ref'] , diff_pd['lon_ref'] , diff_pd['h_ref'] = lat_ref , lon_ref , h_ref
            diff_pd['lat_rov'] , diff_pd['lon_rov'] , diff_pd['h_rov'] = lat_rov , lon_rov , h_rov

            #Merge coordinates results
            coord_res = diff_pd[['STAT_x','STAT_y','epoc','lat_ref','lon_ref','h_ref','lat_rov','lon_rov','h_rov']].copy()
            coord_res = coord_res.rename(index=str,columns={"STAT_x":"STAT_ref","STAT_y":"STAT_rov"})
        else:
            import sys
            print("No this option for coordinates")
            sys.exit()
        #Extract standard ties

        grid = grid_met
        if coord_t == "kinematic":
            diff_pd['stand_ties'] = diff_pd.apply(lambda x: gtro.calc_stand_ties(x['epoc'], x.lat_ref , x.lon_ref , x.h_ref,x.lat_rov , x.lon_rov , x.h_rov,grid),axis=1)
        elif coord_t == "static":
            diff_pd['stand_ties'] = diff_pd.apply(lambda x: gtro.calc_stand_ties(x['epoc'], float(lat_ref) , float(lon_ref) , float(h_ref) , float(lat_rov) , float(lon_rov) , float(h_rov),grid),axis=1)


        diff_pd['Trop_ties_corr'] = diff_pd.apply(lambda x: x['tro_x'] - (x['tro_y']+x['stand_ties']),axis=1)

    #drop unnecessary column
    if coord_t == "kinematic":
        diff_pd = diff_pd.drop(['tro_x', 'stro_x', 'tgn_x', 'stgn_x', 'tge_x', 'stge_x','tro_y', 'stro_y', 'tgn_y', 'stgn_y', 'tge_y','stge_y','lat_ref','lon_ref','h_ref','lat_rov','lon_rov','h_rov'],axis=1)
    else:
        diff_pd = diff_pd.drop(['tro_x', 'stro_x', 'tgn_x', 'stgn_x', 'tge_x', 'stge_x','tro_y', 'stro_y', 'tgn_y', 'stgn_y', 'tge_y','stge_y'],axis=1)

    #Change column name of station
    diff_pd = diff_pd.rename(index=str,columns={"STAT_x":"STAT_ref","STAT_y":"STAT_rov"})

    return diff_pd,coord_res

def stat_summary_trop_ties(df):
    wmean_no_ties = np.round(np.average(df.Trop_ties,weights = 1/df.STrop_ties),3)
    wmean_wt_ties = np.round(np.average(df.Trop_ties_corr,weights = 1/df.STrop_ties),3)
    rms_mean_no_ties = np.round(geok.rms_mean(df.Trop_ties),3)
    rms_mean_wt_ties = np.round(geok.rms_mean(df.Trop_ties_corr),3)

    return [wmean_no_ties,wmean_wt_ties,rms_mean_no_ties,rms_mean_wt_ties]

def plot_trop_ties(df,ref_sta,rov_sta,analy_coor=False,analy_num_obs=False,df_coord="",savePlot=False,filePath="",fileName=""):
    """
    Plot tropospheric ties function
    Input:
        df : DataFrame from "compare_trop_ties" function
        ref_sta : Refernce station
        rov_sta : Rover station
        analy_coor : Add height difference information from "compare_trop_ties"
        analy_num_obs : Add plot number of observation from "compare_trop_ties"
        df_coord : list of height difference
        savePlot : save figure
        filePath : Directory to save
        fileName : Filename of figure
    """
    epo_plt = geok.dt2year_decimal(df.epoc)
    h_diff = df_coord.h_ref - df_coord.h_rov
    if analy_coor:

        fig, ax = plt.subplots(3,1,sharex=True)
        axA = ax[0]
        axB = ax[1]
        axC = ax[2]
        axA.plot(epo_plt,df.Trop_ties , marker="P",linestyle="--",label="Trop.ties")
        axA.plot(np.unique(epo_plt), np.poly1d(np.polyfit(epo_plt, df.Trop_ties, 1))(np.unique(epo_plt)),label="Trop. ties fitline")
        if 'Trop_ties_corr' in df.columns:
            axA.plot(epo_plt,df.Trop_ties_corr,marker="*",linestyle="-.",label="Trop. ties apply height corr.")
            axA.plot(np.unique(epo_plt), np.poly1d(np.polyfit(epo_plt, df.Trop_ties_corr, 1))(np.unique(epo_plt)),label="Trop. ties apply height fitline")
        axA.grid()
        axA.legend()
        axA.set_ylabel("Trop. ties (mm)")
        axA.set_xlabel("Time")
        axA.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))

        axB.plot(epo_plt,h_diff,marker="P",linestyle="--",label="Height difference")
        axB.plot(np.unique(epo_plt), np.poly1d(np.polyfit(epo_plt, h_diff, 1))(np.unique(epo_plt)),label="Height diff. fitline")
        axB.grid()
        axB.legend()
        axB.set_ylabel("Height difference (m)")
        axB.set_xlabel("Time")
        axB.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        if analy_num_obs:
            if len(epo_plt) != len(df.num_obs_ref):
                print("No plot number of observations")
                fig.delaxes(axC)
            else:
                axC.plot(epo_plt,df.num_obs_ref,marker="P",linestyle="--",label="Num obs. Ref sta")
                axC.plot(epo_plt,df.num_obs_rov,marker="*",linestyle="-.",label="Num obs. Rov sta")
                axC.grid()
                axC.legend()
                axC.set_ylabel("Num Obs.")
                axC.set_xlabel("Time")
                axC.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    else:
        fig , ax = plt.subplots(2,1,sharex=True)
        axA = ax[0]
        axC = ax[1]
        axA.plot(epo_plt,df.Trop_ties , marker="P",linestyle="--",label="Trop.ties")
        axA.plot(np.unique(epo_plt), np.poly1d(np.polyfit(epo_plt, df.Trop_ties, 1))(np.unique(epo_plt)),label="Trop. ties fitline")
        if 'Trop_ties_corr' in df.columns:
            axA.plot(epo_plt,df.Trop_ties_corr,marker="*",linestyle="-.",label="Trop. ties apply height corr.")
            axA.plot(np.unique(epo_plt), np.poly1d(np.polyfit(epo_plt, df.Trop_ties_corr, 1))(np.unique(epo_plt)),label="Trop. ties apply height fitline")
        axA.grid()
        axA.legend()
        axA.set_ylabel("Trop. ties (mm)")
        axA.set_xlabel("Time")
        axA.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))

        if analy_num_obs:
            if len(epo_plt) != len(df.num_obs_ref):
                print("No plot number of observations")
                fig.delaxes(axC)
            else:
                axC.plot(epo_plt,df.num_obs_ref,marker="P",linestyle="--",label="Num obs. Ref sta")
                axC.plot(epo_plt,df.num_obs_rov,marker="*",linestyle="-.",label="Num obs. Rov sta")
                axC.grid()
                axC.legend()
                axC.set_ylabel("Num Obs.")
                axC.set_xlabel("Time")
                axC.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))

    plt.tight_layout()
    plt.suptitle("Total delay ties of " + ref_sta + "-" + rov_sta)
    if savePlot:
        export_ts_figure_pdf(fig,filePath,fileName,True)

    plt.show()

    return None

##########################################################################################

#def stations_in_sinex_multi(sinex_path_list):
#    """
#    Gives stations list in a SINEX
#    Args :
#        sinex_path_list : path of the SINEX files list (e.g. made with glob)
#    Returns :
#        datadico : a dico with stat 4 char. code as key and epoch list as values (for timeline_plotter)
#    """
#    epoch_list_stk , stats_list_stk = [] , []
#    for sinex_path in sinex_path_list:
#        try:
#            epoch , stats_list = stations_in_sinex_mono(sinex_path)
#            epoch_list = [epoch] * len(stats_list)
#            epoch_list_stk = epoch_list_stk + epoch_list
#            stats_list_stk = stats_list_stk + stats_list
#        except:
#            continue
#
#    datadico = dict()
#    for e,s in zip(epoch_list_stk,stats_list_stk):
#        if not s in datadico.keys():
#            datadico[s] = []
#        else:
#            datadico[s].append(e)
#
#    for v in datadico.values():
#        v.sort()
#
#    return datadico

def stations_in_coords_file_multi(files_path_list,files_type = "sinex"):
    """
    Gives stations list in a SINEX or EPOS coords. file

    Args :
        files_path_list : path of the SINEX/coords files list (e.g. made with glob)
        files_type       : string, "sinex" or "EPOS_sta_coords"

    Returns :
        datadico : a dico with stat 4 char. code as key and epoch list as values (for timeline_plotter)
    """

    if   files_type == "sinex":
        extract_fct = stations_in_sinex_mono
    elif files_type == "EPOS_sta_coords":
        extract_fct = stations_in_EPOS_sta_coords_file_mono
    else:
        print("ERR : check station file type !!")
        return None

    epoch_list_stk , stats_list_stk = [] , []

    for file_path in files_path_list:
        try:
            epoch , stats_list = extract_fct(file_path)
            epoch_list = [epoch] * len(stats_list)
            epoch_list_stk = epoch_list_stk + epoch_list
            stats_list_stk = stats_list_stk + stats_list
        except:
            print('WARN : something wrong with' , file_path)
            print('TIPS : good files type ? ("sinex" or "EPOS_sta_coords") : ', files_type)

            continue

    datadico = dict()
    for e,s in zip(epoch_list_stk,stats_list_stk):
        if not s in datadico.keys():
            datadico[s] = []
        else:
            datadico[s].append(e)

    for v in datadico.values():
        v.sort()

    return datadico



def stations_in_sinex_multi(sinex_path_list):
    """
    Gives stations list in a SINEX

    Just a wrapper of stations_in_sinex_or_EPOS_sta_coords_file_multi !!!
    This other function should be used !!!

    Args :
        sinex_path_list : path of the SINEX files list (e.g. made with glob)
    Returns :
        datadico : a dico with stat 4 char. code as key and epoch list as values (for timeline_plotter)
    """

    return stations_in_coords_file_multi(sinex_path_list,"sinex")



def sinex_bench_antenna_DF_2_disconts(DFantenna_in,stat,return_full=False):
    DFantenna_work = DFantenna_in[DFantenna_in["Code"] == stat]
    Start_List     = geok.datestr_sinex_2_dt(DFantenna_work["_Data_Start"])
    End_list       = geok.datestr_sinex_2_dt(DFantenna_work["_Data_End__"])
    if return_full:
        return Start_List,End_list
    else:
        Clean_list = sorted(list(set(Start_List + End_list)))
        Clean_list = [e for e in Clean_list if e != dt.datetime(1970, 1, 1, 0, 0)]
        return Clean_list


def rinex_check_epochs_availability(rinex_path_list):
    """
    Args :
        A list of rinex paths
    Returns :
        T : a table with results
    """

    results_stk = []

    for rinex_path in rinex_path_list:

        rinex_name = os.path.basename(rinex_path)

        QC = softs_runner.teqc_qc(rinex_path)

        if not QC:
            continue

        epoc_all  = int(genefun.egrep_big_string("Poss. # of obs epochs" ,QC,only_first_occur=True).split()[-1])
        epoc_disp = int(genefun.egrep_big_string("Epochs w/ observations",QC,only_first_occur=True).split()[-1])

        dt_rnx = geok.rinexname2dt(rinex_name)

        date_str = geok.dt2str(dt_rnx,"%F")

        percentage = (float(epoc_disp) / float(epoc_all)) * 100.

        results = [rinex_name,date_str,epoc_disp,epoc_all,percentage]

        results_stk.append(results)

    header = ['RINEX','date','Avbl.', 'Poss.', '%']
    T = tabulate.tabulate(results_stk,headers=header)

    return T

def write_sndy_light_dat(ts_in,outdir,outprefix):
    """pas fini"""
    fil = open(os.path.join(outdir,outprefix),'w+')
    if isinstance(ts_in,TimeSeriePoint):
        if ts_in.initype() == 'FLH':
            for pt in ts_in.pts:
                lin = ' '.join([str(e) for e in [pt.F , pt.L , pt.H , pt.T , pt.sF , pt.sL , pt.sH ]])
                fil.write(lin + '\n')
    elif isinstance(ts_in,TimeSerieObs):
        if ts_in.typeobs == 'RPY':
            for att in ts_in.obs:
                lin = ' '.join([str(e) for e in [att.R , att.P , att.Y , att.T , att.Q.w , att.Q.x , att.Q.y , att.Q.z ]])
                fil.write(lin + '\n')
    fil.close()


def mean_posi_multi(tstup):
    mp = tstup[0].mean_posi()
    for ts in tstup:
        ts.ENUcalc(mp)
    return None

def refENU_for_tslist(tslist_in, tsref_marker = 0):
    """
    tsref_marker : indice of the reference time serie OR the 'all' keyword
    in this case all the time series mean position will be averaged
    """
    if tsref_marker == 'all':
        refENU_stk =  []
        for ts in tslist_in:
            refENUtmp = ts.mean_posi()
            refENU_stk.append(refENUtmp)

        X = np.mean([renu.X for renu in refENU_stk])
        Y = np.mean([renu.Y for renu in refENU_stk])
        Z = np.mean([renu.Z for renu in refENU_stk])

        refENU = Point(X,Y,Z,0,initype='XYZ')

    else:
        refENU = tslist_in[tsref_marker].mean_posi()


    for ts in tslist_in:
        ts.ENUcalc(refENU)

    return refENU



def time_win_multi(inplis):
    strt = inplis[0].startdate()
    end  = inplis[0].enddate()
    tsoutlis = [inplis[0]]
    for ts in list(inplis)[1:]:
        tsout = time_win(ts,[[strt,end]])
        tsoutlis.append(tsout)
    return tsoutlis

def ts_from_list(A,B,C,T,initype,sA=[],sB=[],sC=[],stat='STAT',name='NoName'):
    tsout = TimeSeriePoint()
    if len(sA) == 0:
        for a,b,c,t in zip(A,B,C,T):
            pt = Point(a,b,c,t,initype)
            tsout.add_point(pt)
    else:
        for a,b,c,t,sa,sb,sc in zip(A,B,C,T,sA,sB,sC):
            pt = Point(a,b,c,t,initype,sa,sb,sc)
            tsout.add_point(pt)
    tsout.stat = stat
    tsout.name = name
    if initype == 'ENU':
        tsout.boolENU = True
    else:
        tsout.boolENU = False
    return tsout





# ================================== FONCTION DE TRAITEMENT ==================================

def print4compar(dA,dB,dC,dD,coortype):

    if coortype == 'ENU':
        Astr,Bstr,Cstr = 'E','N','U'
    elif coortype == 'XYZ':
        Astr,Bstr,Cstr = 'X','Y','Z'
    else:
        Astr,Bstr,Cstr = 'A','B','C'

    print("moyenne aritm & RMS et std composante " + Astr)
    print(str(np.nanmean(dA)))
    print(str(geok.RMSmean(dA)))
    print(str(np.nanstd(dA)))
    print('')

    print("moyenne aritm & RMS et std composante " + Bstr)
    print(str(np.nanmean(dB)))
    print(str(geok.RMSmean(dB)))
    print(str(np.nanstd(dB)))
    print('')

    print("moyenne aritm & RMS et std composante " + Cstr)
    print(str(np.nanmean(dC)))
    print(str(geok.RMSmean(dC)))
    print(str(np.nanstd(dC)))
    print('')

    print("moyenne aritm & RMS et std composante D")
    print(str(np.nanmean(dD)))
    print(str(geok.RMSmean(dD)))
    print(str(np.nanstd(dD)))
    print('')

    print("RMS3D : sqrt((RMS_{}**2 + RMS_{}**2 + RMS_{}**2)/3 ) ".format(Astr,Bstr,Cstr))
    print(geok.RMSmean([geok.RMSmean(dA),geok.RMSmean(dB),geok.RMSmean(dC)]))
    print("RMS2D : uniquement sur les 2 composantes plani")
    print(geok.RMSmean([geok.RMSmean(dA),geok.RMSmean(dB)]))
    print('')

def print4compar_tabular(dicolist,split=0,print_2D3D_if_any=True):

    dAlist,dBlist,dClist,dDlist,dDblist,statlist = [],[],[],[],[],[]
    for dic in dicolist:
        if "dD2D" in list(dic.keys()):
            D2Dn3D = True
        else:
            D2Dn3D = False

        if not print_2D3D_if_any:
            D2Dn3D = False

        coortype = dic['coortype']
        Dtype = dic['Dtype']

        dAlist.append(dic['dA'])
        dBlist.append(dic['dB'])
        dClist.append(dic['dC'])

        if not D2Dn3D:
            dDlist.append(dic['dD'])
        else:
            dDlist.append(dic['dD2D'])
            dDblist.append(dic['dD3D'])

        statlist.append(dic['name'])

    if coortype == 'ENU':
        Astr,Bstr,Cstr = 'E','N','U'
    elif coortype == 'XYZ':
        Astr,Bstr,Cstr = 'X','Y','Z'
    else:
        Astr,Bstr,Cstr = 'A','B','C'

    datacol = ("Arithm. mean", "RMS mean", "std. dev.")
    datacol = ("Arit.", "RMS", geok.greek_alphabet()[17])

    headerline = ['Experiment']

    if not D2Dn3D:
        headervartup = (Astr,Bstr,Cstr,'D ' + Dtype)
    else:
        headervartup = (Astr,Bstr,Cstr,'D 2D', 'D 3D')

    for abc in headervartup:
        for dc in datacol:
            headerline.append(' '.join((dc,abc)))

#    headerline.append("RMS3D : sqrt((RMS_{}**2 + RMS_{}**2 + RMS_{}**2)/3 )".format(Astr,Bstr,Cstr))
#    headerline.append("RMS2D : uniquement sur les 2 composantes plani")

    headerline.append("RMS3D")
    headerline.append("RMS2D")

    LINES_STK = [headerline]

    if not D2Dn3D:
        dDblist = [None] * len(dDlist)

    for (dA,dB,dC,dD,dDb,stat) in zip(dAlist,dBlist,dClist,dDlist,dDblist,statlist):

        if np.isclose(np.sum(dA) + np.sum(dB) + np.sum(dC),0):
            continue

        line = []
        line.append(stat)

        line.append(np.nanmean(dA))
        line.append(geok.RMSmean(dA))
        line.append(np.nanstd(dA))

        line.append(np.nanmean(dB))
        line.append(geok.RMSmean(dB))
        line.append(np.nanstd(dB))

        line.append(np.nanmean(dC))
        line.append(geok.RMSmean(dC))
        line.append(np.nanstd(dC))

        line.append(np.nanmean(dD))
        line.append(geok.RMSmean(dD))
        line.append(np.nanstd(dD))

        if D2Dn3D:
            line.append(np.nanmean(dDb))
            line.append(geok.RMSmean(dDb))
            line.append(np.nanstd(dDb))

        line.append(geok.RMSmean([geok.RMSmean(dA),geok.RMSmean(dB),geok.RMSmean(dC)]))
        line.append(geok.RMSmean([geok.RMSmean(dA),geok.RMSmean(dB)]))

        LINES_STK.append(line)

    LINES_STK_orig = list(LINES_STK)

    if split != 0:
        LINES_arr = np.vstack(LINES_STK)
        LINES1 = LINES_arr[:,:split]
        LINES2 = np.column_stack((LINES_arr[:,0],LINES_arr[:,split:]))
        return LINES1,LINES2
    else:
        return LINES_STK

def compar_plot(dico_list_in, namest = 0, namend    = 10  ,
                alpha = 0.8 , diapt = 1.5 , new_style = True,
                colormap = 'gnuplot'):


    if not new_style:
        fig , ax_arr = plt.subplots(2,2)
        fig.set_size_inches(16.53,11.69)
    else:
        fig , ax_arr = plt.subplots(4,1)
        fig.set_size_inches(11.69,16.53)

    cm = plt.get_cmap(colormap)
    NUM_COLORS = len(dico_list_in)
    color_list = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]


    for color , D in zip(color_list,dico_list_in):

        nametmp  = D['name']
        nametmp  = nametmp.ljust(namend - namest)
        coortype = D['coortype']
        Dtype    = D['Dtype']

        if coortype == 'ENU':
            Ati = '\u0394' + 'East component'
            Bti = '\u0394' + 'North component'
            Cti = '\u0394' + 'Up component'
            Dti = '\u0394' + 'Distance (' + Dtype + ') component'
            ylbl = 'meters'
        elif coortype == 'XYZ':
            Ati = 'delta X'
            Bti = 'delta Y'
            Cti = 'delta Z'
            Dti = 'delta D ' + Dtype
            ylbl = 'meters'
        elif coortype == 'FLH':
            Ati = 'delta Latitude'
            Bti = 'delta Longitude'
            Cti = 'delta Height'
            Dti = 'delta D ' + Dtype + ' (no sens)'
            ylbl = 'degrees'
        else:
            Ati = 'delta A'
            Bti = 'delta B'
            Cti = 'delta C'
            Dti = 'delta D ' + Dtype
            ylbl = '???'

        if not new_style:
            axA = ax_arr[0][0]
            axB = ax_arr[0][1]
            axC = ax_arr[1][0]
            axD = ax_arr[1][1]
        else:
            axA = ax_arr[0]
            axB = ax_arr[1]
            axC = ax_arr[2]
            axD = ax_arr[3]

        axA.plot(D['TA'],D['dA'], '.-' ,label=nametmp[namest:namend] ,
                 markersize=diapt, alpha=alpha,color=color)
        axA.set_ylabel(ylbl)
        axA.legend()
        axA.set_title(Ati)

        axB.plot(D['TB'] ,D['dB'] ,'.-',label=nametmp[namest:namend] ,
                 markersize=diapt, alpha=alpha,color=color)
        axB.set_ylabel(ylbl)
        axB.legend()
        axB.set_title(Bti)

        axC.plot(D['TC'] ,D['dC'],'.-',label=nametmp[namest:namend] ,
                 markersize=diapt, alpha=alpha,color=color)
        axC.set_ylabel(ylbl)
        axC.legend()
        axC.set_title(Cti)

        axD.plot(D['TD'] ,D['dD'],'.-',label=nametmp[namest:namend] ,
                 markersize=diapt, alpha=alpha,color=color)
        axD.set_ylabel(ylbl)
        axD.legend()
        axD.set_title(Dti)

    return fig

def compar(tstup , coortype='ENU' , seuil=3. , win=[] , mode='keep' ,
           Dtype='2D3D', namest=0,namend=10,alpha = 5 , diapt = 5 , verbose=True,
           print_report=True,plot=True):
    """
        si seuil == 0, pas de nettoyage

        On preconise compar en mode keep :
        Ainsi le Tref est strictement cantonné aux bonnes valeurs des Ti
        En mode del, les Ti sont extrapolé aux valeur de Tref => nan => pbs pour la mad


        if 2D3D /3D2D : selectionne 2D ou 3D (selon) en prioritaire pour dD mais
        détermine 2D et 3D et les introduit dans des clés adhoc dD2D dD3D
        '''
    """

    if len(Dtype) == 4:
        D2n3 = True
        Dtype = Dtype[0:2]
    else:
        D2n3 = False

    diapt=10
    alpha=0.8

    if len(tstup) <= 1:
        print("ERR : compar : len(tstup) <= 1 , on ne compar pas un élément seul ;) !!!!")

    if verbose:
        print("")
        print("BILAN DE LA COMPARAISON")
        print("reference : " +  tstup[0].name)
        print("coordonnées : " +  coortype)
        print("")


    dicolist = []

    for ivar,tsvar in enumerate(tstup):


        #dico de sortie
        dicovar = dict()

        if verbose:
            print("===========================================================")
            print(tsvar.name)
            print("------------------------------")


#        tsout = TimeSeriePoint()
#        tsout = copy.copy(tsvar)
#        tsout.del_data()

        if tsvar.bool_interp_uptodate == False:
            tsvar.interp_set()

        A, B, C, T, sA , sB , sC = tsvar.to_list(coortype=coortype)

        tsref = tstup[0].interp_get(T,coortype=coortype)
        Aref, Bref, Cref , Tref , sA , sB , sC = tsref.to_list(coortype=coortype)


        # WTF IS THOSE LINES
        # autowin = time_gap(tstup[0],marge=1.1,mode='keep')
        # win = autowin
        #
        # outbool = []


        # (winlis == [] or (ivar in winlis)) and
        if win != []:
            if verbose:
                print("Application du fenetrage")
                print("------------------------------")
            bl = time_win_T(T,win,mode=mode)
            if verbose:
                print("Nb de pts : " , len(bl))
                print("Pts validés :" , np.sum(bl))
                print('')

        else:
            bl = np.asarray([True] * len(Tref))

        dA = A[bl] - Aref[bl]
        dB = B[bl] - Bref[bl]
        dC = C[bl] - Cref[bl]
        Tref = Tref[bl]

        #print "N of diff close to 0 => Bug ..."
        #print 'dA' , np.sum(np.isclose(dA,0))
        #print 'dB' , np.sum(np.isclose(dB,0))
        #print 'dC' , np.sum(np.isclose(dC,0))


        dD2D = np.sqrt(dA ** 2 + dB ** 2)
        dD3D = np.sqrt(dA ** 2 + dB ** 2 + dC ** 2)

        if Dtype == '2D':
            dD = np.sqrt(dA ** 2 + dB ** 2)
        elif Dtype == '3D':
            dD = np.sqrt(dA ** 2 + dB ** 2 + dC ** 2)

        else:
            print("ERR : compar : wrong Dtype")


#        print "Stats sans Nettoyage"
#        print "------------------------------"
#        print4compar(dA,dB,dC,dD,coortype)

        # nettoyage par la MAD
        if seuil != 0:
            if verbose:
                print("Nettoyage par la MAD")
                print("------------------------------")
            dAin = dA
            dA,bb = geok.outiler_mad(dA,seuil=seuil)
            dAout = dA
            TA = geok.posix2dt(np.array(Tref[bb]))

            dBin = dB
            dB,bb = geok.outiler_mad(dB,seuil=seuil)
            dBout = dB
            TB = geok.posix2dt(np.array(Tref[bb]))

            dCin = dC
            dC,bb = geok.outiler_mad(dC,seuil=seuil)
            dCout = dC
            TC = geok.posix2dt(np.array(Tref[bb]))

            dDin = dD
            dD,bb = geok.outiler_mad(dD,seuil=seuil)
            TD = geok.posix2dt(np.array(Tref[bb]))
            dDout = dD
            if D2n3:
                dD2Din = dD2D
                dD2D,bb = geok.outiler_mad(dD2D,seuil=seuil)
                TD2D = geok.posix2dt(np.array(Tref[bb]))
                dD2Dout = dD2D
                dD3Din = dD3D
                dD3D,bb = geok.outiler_mad(dD3D,seuil=seuil)
                TD3D = geok.posix2dt(np.array(Tref[bb]))
                dD3Dout = dD3D

            if verbose:
                print('')
                print("Stats aprÃ¨s Nettoyage")
                print("------------------------------")
            if print_report:
                print4compar(dA,dB,dC,dD,coortype)

        else:
            Tdt = geok.posix2dt(Tref)
            TA,TB,TC,TD = Tdt,Tdt,Tdt,Tdt
            dAin,dBin,dCin,dDin = dA,dB,dC,dD

        # pour la legende : les indices de debut et fin du nom

        # fabrication du Dico de sortie
        dicovar['name']     = tsvar.name
        dicovar['coortype'] = coortype
        dicovar['Dtype']    = Dtype

        dicovar['TA'] = np.array(TA)
        dicovar['dA'] = np.array(dA)
        dicovar['dAbrut'] = np.array(dAin)
        dicovar['TB'] = np.array(TB)
        dicovar['dB'] = np.array(dB)
        dicovar['dBbrut'] = np.array(dBin)
        dicovar['TC'] = np.array(TC)
        dicovar['dC'] = np.array(dC)
        dicovar['dCbrut'] = np.array(dCin)
        dicovar['TD'] = np.array(TD)
        dicovar['dD'] = np.array(dD)
        dicovar['dDbrut'] = np.array(dDin)
        dicovar['dDtype'] = np.array(Dtype)
        if D2n3:
            dicovar['TD2D'] = np.array(TD2D)
            dicovar['dD2D'] = np.array(dD2D)
            dicovar['dD2Dbrut'] = np.array(dD2Din)
            dicovar['dD2Dtype'] = np.array('2D')

            dicovar['TD3D'] = np.array(TD3D)
            dicovar['dD3D'] = np.array(dD2D)
            dicovar['dD3Dbrut'] = np.array(dD3Din)
            dicovar['dD3Dtype'] = np.array('3D')

        dicolist.append(dicovar)

    if plot:
        compar_plot(dicolist,namest,namend,alpha,diapt)
        try:
            suptit = ' '.join((tstup[-1].stat , str(tstup[-1].startdate()) ,
                               str(tstup[-1].enddate())))
            plt.suptitle(suptit)
        except:
            pass

    return dicolist

def compar2(tstup,coortype='ENU',seuil=3.,win=[],mode='keep',
           Dtype='3D',namest=0,namend=10,alpha = 5 , diapt = 5 ,verbose=True,
           print_report=True,plot=True):
    """
    160903 cette fonction semble discontinuée
    """

    tsref = tstup[0]
    dicolist = []

    for tsvar in tstup:
        Aref,Bref,Cref,Tref,_,_,_  = tsref.to_list()
        Avar,Bvar,Cvar,Tvar,_,_,_  = tsvar.to_list()

        Tref = np.round(Tref,1)
        Tvar = np.round(Tvar,1)

        Tcommon    = np.intersect1d(Tref,Tvar)
        Tcommon_dt = geok.posix2dt(Tcommon)

        ind_ref = np.nonzero(np.in1d(Tref , Tcommon))[0]
        ind_var = np.nonzero(np.in1d(Tvar , Tcommon))[0]

        dA = Avar[ind_var] - Aref[ind_ref]
        dB = Bvar[ind_var] - Bref[ind_ref]
        dC = Cvar[ind_var] - Cref[ind_ref]

        if Dtype == '2D':
            dD = np.sqrt(dA ** 2 + dB ** 2)
        elif Dtype == '3D':
            dD = np.sqrt(dA ** 2 + dB ** 2 + dC ** 2)

        dicovar = dict()

        dicovar['name']     = tsvar.name
        dicovar['coortype'] = coortype
        dicovar['Dtype']    = Dtype

        dicovar['TA'] = Tcommon_dt
        dicovar['dA'] = np.array(dA)
        dicovar['dAbrut'] = np.array(dA)
        dicovar['TB'] = Tcommon_dt
        dicovar['dB'] = np.array(dB)
        dicovar['dBbrut'] = np.array(dB)
        dicovar['TC'] = Tcommon_dt
        dicovar['dC'] = np.array(dC)
        dicovar['dCbrut'] = np.array(dC)
        dicovar['TD'] = Tcommon_dt
        dicovar['dD'] = np.array(dD)
        dicovar['dDbrut'] = np.array(dD)
        dicovar['dDtype'] = np.array(Dtype)

        dicolist.append(dicovar)

    if plot:
        compar_plot(dicolist,namest,namend,alpha,diapt)
        try:
            suptit = ' '.join((tstup[-1].stat , str(tstup[-1].startdate()) ,
                               str(tstup[-1].enddate())))
            plt.suptitle(suptit)
        except:
            pass

    return dicolist


def mad_cleaner(tsin,seuil=3.5,method='dist',coortype='ABC',
                detrend_first=False,output_detrended=False, verbose=False):
    '''
    method : methode d'élimination :
            dist : on élimine les point qu sont trop loin en distance de la posi de ref
            indep : on traite les point independaments

            dist est a privilégier
    output_detrended ne marche que si detrend_first est activé
    '''

    if coortype == 'ABC':
        coortype = tsin.initype()
    print(coortype)

    if detrend_first:
        tswork = detrend_ts(tsin,coortype)
    else:
        tswork = tsin

    A,B,C,T,sA,sB,sC = tswork.to_list(coortype)

    D = np.sqrt(A ** 2 + B ** 2 + C ** 2)

    if method == 'dist':
        Dout , bb = geok.outiler_mad(D,seuil)
    if method == 'indep':
        Aout , bbA = geok.outiler_mad(A,seuil)
        Bout , bbB = geok.outiler_mad(B,seuil)
        Cout , bbC = geok.outiler_mad(C,seuil)

        bb = bbA * bbB * bbC

    if output_detrended:
        tsout = bool_cleaner(tswork,bb)
    else:
        tsout = bool_cleaner(tsin,bb)

    if verbose:
        killratio = np.round(float(tsout.nbpts) / float(tsin.nbpts),4)

        print("INFO : mad_cleaner :" , tsin.stat , "clean ratio" , tsin.nbpts , tsout.nbpts , killratio)

    return tsout

def sigma_cleaner(tsin,seuil=3,coortype='ABC',cleantype='any', verbose=False):

    if coortype == 'ABC':
        coortype = tsin.initype()

    A,B,C,T,sA,sB,sC = tsin.to_list(coortype)

    sAout , bbsA = geok.outlier_sigma(sA,seuil)
    sBout , bbsB = geok.outlier_sigma(sB,seuil)
    sCout , bbsC = geok.outlier_sigma(sC,seuil)

    if cleantype == 'any':
        boolbad = bbsA * bbsB * bbsC
    elif cleantype == 'all':
        boolbad = bbsA + bbsB + bbsC
    tsout = bool_cleaner(tsin,boolbad)

    if verbose:
        killratio = np.round(float(tsout.nbpts) / float(tsin.nbpts),4)

        print("INFO : sigma_cleaner :" , tsin.stat , "clean ratio" , tsin.nbpts , tsout.nbpts , killratio)

    return tsout

def std_dev_cleaner(tsin,stddev_threshold,coortype="ABC",cleantype="any",
                    verbose=False):
    """
    A rebooted (1807) version of sigma_cleaner
    just remove values in a timeserie with a high sigma/std deviation
    """

    if coortype == 'ABC':
        coortype = tsin.initype()

    A,B,C,T,sA,sB,sC = tsin.to_list(coortype)

    bbsA = sA <= stddev_threshold
    bbsB = sB <= stddev_threshold
    bbsC = sC <= stddev_threshold

    if cleantype == 'any':
        boolbad = bbsA * bbsB * bbsC
    elif cleantype == 'all':
        boolbad = bbsA + bbsB + bbsC
    tsout = bool_cleaner(tsin,boolbad)

    if verbose:
        killratio = np.round(float(tsout.nbpts) / float(tsin.nbpts),4)

        print("INFO : std_dev_cleaner :" , tsin.stat , "clean ratio" , tsin.nbpts , tsout.nbpts , killratio)

    return tsout

def bool_cleaner(tsin,boollist, verbose=False):
    ''' A partir d'une liste de bool de meme longeur que le nbre de points
        on ne conserve que les points True '''
    tsout = copy.copy(tsin)
    ptslistin = tsin.pts
    if len(ptslistin) != len(boollist):
        print("ERR : bool_cleaner : len bool != len pts")
        return 0
    ptslistout = []
    for pt,b in zip(ptslistin,boollist):
         if b:
             ptslistout.append(pt)
    tsout.pts = ptslistout

    if verbose:
        killratio = np.round(float(tsout.nbpts) / float(tsin.nbpts),4)

        print("INFO : bool_cleaner :" , tsin.stat , "clean ratio" , tsin.nbpts , tsout.nbpts , killratio)

    return tsout

def linear_regress_find_coeff(tsin,coortype='ENU'):
    A, B, C , T , sA , sB , sC = tsin.to_list(coortype)
    outW = []
    for i,composante in enumerate((A,B,C)):
        M = np.array([ T , np.ones(len(T)) ])
        w  = scipy.linalg.lstsq(M.T,composante)[0]
        outW.append((w[0],w[1]))
    return outW

def detrend_ts(tsin,coortype='ENU'):
    A, B, C , T , sA , sB , sC = tsin.to_list(coortype)
    try:
        Xref , Yref , Zref = tsin.mean_posi()
    except:
        Xref , Yref , Zref = 0,0,0

    W = linear_regress_find_coeff(tsin,coortype)
    UnTr = [] # UnTr
    for i,composante in enumerate((A,B,C)):
        w = W[i]
        comp_line = w[0]*T+w[1]
        comp_untrend = composante - comp_line
        UnTr.append(comp_untrend)

    tsout = ts_from_list(UnTr[0],UnTr[1],UnTr[2],T,coortype,sA,sB,sC,tsin.stat,tsin.name)




    tsout.anex['Xref'] = Xref
    tsout.anex['Yref'] = Yref
    tsout.anex['Zref'] = Zref

    return tsout

def linear_regress_ts(tsin,coortype='ENU',titledetails = ''):
    """ doit être cablé ASAP linear_regress_find_coeff"""
    try:
        A, B, C , T , sA , sB , sC = tsin.to_list(coortype)
        fig , axes = plt.subplots(2,2)
        fig.suptitle(tsin.stat + ' ' + titledetails)
        axes = np.reshape(axes,-1)

        W = linear_regress_find_coeff(tsin,coortype)

        for i,composante in enumerate((A,B,C)):
            M = np.array([ T , np.ones(len(T)) ])
#            w  = scipy.linalg.lstsq(M.T,composante)[0]
            w = W[i]
            xi = np.arange(min(T),max(T),10**5)
            line = w[0]*xi+w[1]
            axes[i].plot(geok.posix2dt(xi),line,'r-')
            axes[i].plot(geok.posix2dt(T),composante,'+')
            axes[i].set_title(coortype[i])
            axes[i].set_xlabel('Date')
            axes[i].set_ylabel('Displacement (m)')
            v = w[0] * 31556926 / 10**-3
            v2 = round(v,2)
            print('composante ' + coortype[i] + ' = ' + str(v) + 'mm/year')
            axes[i].text(0.9, 0.1, str(v2) + 'mm/year' , horizontalalignment='center',verticalalignment='center',transform=axes[i].transAxes)
    except:
        pass

def linear_regress_ts_discont(tsin,coortype = 'ENU') :
    discont_improved = tsin.discont + [dt.datetime(2099,1,1)]
    for j in range(len(discont_improved)-1):
        print('période ',j,':',discont_improved[j],discont_improved[j+1])
        tsbis = time_win(tsin,[[discont_improved[j],discont_improved[j+1]]])
        datetitle = ' : ' + str(discont_improved[j]) + "\u2192" + str(discont_improved[j+1])
        linear_regress_ts(tsbis,coortype,titledetails=datetitle)

#def export_ts_as_pbo_pos(tsin):
#    # IN PROGRESS
#    tswork = copy.deepcopy(tsin)
#    outtfile.write('{:.5f}   {:+.6f}    {:+.6f}    {:+.6f} {:+.6f} {:+.6f} {:+.6f}\n'.format(t,n-n0,e-e0,u-u0,se,sn,su))
#

def export_ts_figure_pdf(fig,export_path,filename,close=False):
    """ fig can accept a int (id of a Figure)
         OR the figure Object itself """

    if type(fig) is int:
        f = plt.figure(fig)
    elif type(fig) is Figure:
        f = fig

    f.set_size_inches(16.53,11.69)
    # tralala pour avoir des dates propres
    # parce que dans des cas + simples f.autofmt_xdate() suffit
    for a in f.axes[1:]:
        labels = a.get_xticklabels()
        for l in labels:
            l.set_rotation(40)
            l.set_horizontalalignment('right')
    out_path = os.path.join(export_path,filename+'.pdf')
    f.tight_layout()
    f.subplots_adjust(top=0.94)
    f.savefig(out_path,papertype='a4',format='pdf')

    if close:
        plt.close(f)

    return None

def export_ts_plot(tsin,export_path,coortype='ENU',export_type=("pdf","png"),
                   plot_B = False,close_fig_after_export=True):
    """ Very beta ...
        to be implemented : merge w/ the export_figure_pdf fct """
    # plot A avec les barres de sigma
    plt.clf()
    tsin.plot(coortype,fig=1)
    tsin.plot_discont()
    f = plt.gcf()
    f.set_size_inches(16.53,11.69)
    # tralala pour avoir des dates propres
    # parce que dans des cas + simples f.autofmt_xdate() suffit
    for a in f.axes[1:]:
        labels = a.get_xticklabels()
        for l in labels:
            l.set_rotation(40)
            l.set_horizontalalignment('right')
    f.tight_layout()
    f.subplots_adjust(top=0.92)
    for typ in export_type:
        export_file = os.path.join(export_path,tsin.stat+'a.' + typ)
        f.savefig(export_file,
                  papertype='a4',format=typ)
        print("INFO : plot exported in " , export_file)
    if close_fig_after_export:
        plt.close(f)

    # plot B avec une droite de regression
    if plot_B:
        linear_regress_ts(tsin)
        f = plt.gcf()
        #    f.set_dpi(300)
        f.set_size_inches(16.53,11.69)
        for a in f.axes:
            labels = a.get_xticklabels()
            for l in labels:
                l.set_rotation(40)
                l.set_horizontalalignment('right')
        f.tight_layout()
        f.subplots_adjust(top=0.92)
        for typ in export_type:
            export_file = os.path.join(export_path,tsin.stat+'b.' + typ)
            f.savefig(export_file,
                      papertype='a4',format=typ)
        print("INFO : plot exported in " , export_file)
        if close_fig_after_export:
            plt.close(f)
    return None


def export_ts(ts,outdir,coordtype = 'ENU',outprefix='',write_header=False):
    """
    export the timeserie

    write_header not well implemented !!!
    """

    proto_str = '{:23} {:23} {:23} {:23} {:23} {:23} {:23} {} {:23} {:23} {}'

    A,B,C,T,sA,sB,sC = ts.to_list(coordtype)
    if outprefix != '':
        outprefix = outprefix + '_'

    outfilenam = outprefix + ts.stat + '.' + coordtype + '.ts.dat'
    outpath = os.path.join(outdir,outfilenam)

    filobj = open(outpath,'w+')

    if write_header:
        "#" + proto_str.format(("year","month","day","hour","minute","seconds","year_decimal","posix_time"))
        filobj.write()

    for a,b,c,t,sa,sb,sc in zip(A,B,C,T,sA,sB,sC):
        tt = geok.posix2dt(t)
        paramstr = [str(e) for e in [a,b,c,t,sa,sb,sc]]
        #paramstr = [e.ljust(18, '0') for e in paramstr]

        yr_dec_str = str(geok.dt2year_decimal(tt))
        posix_str = str(t)

        paramstr2 = paramstr + [tt.strftime("%Y %m %d %H %M %S"),yr_dec_str,posix_str ,'\n']


        outlin =  proto_str.format(*paramstr2)
        filobj.write(outlin)
    filobj.close()

    print('INFO : timeserie exported in ' + outpath)
    return None


def export_ts_as_neu(tsin,outdir,outprefix,coordtype = 'ENU'):
    """
    export to a HECTOR .neu compatible format

    outfile will be writed in
    /outdir/outprefixSTAT.neu
    """
    if not hasattr(tsin[0],'X'):
        print('WARN : export_ts_as_neu : no XYZ in ts')
        noXYZ = True
    else:
        noXYZ = False

    tswork = copy.deepcopy(tsin)
    if coordtype == 'XYZ':
        mp = tswork.mean_posi()
        tswork.ENUcalc(mp)
    outpath = outdir +'/' + outprefix + tswork.stat + '.neu'
    outfile = open(outpath,'w+')
    E,N,U,T,sE,sN,sU = tswork.to_list('ENU')
    first_pt = tswork.pts[0]
    if noXYZ:
        first_pt.X = 0.
        first_pt.Y = 0.
        first_pt.Z = 0.
        first_pt.L = 0.
        first_pt.H = 0.
        first_pt.F = 0.


    e0,n0,u0,t0 = list(zip(E,N,U,T))[0]
    # write the header
    outfile.write('# Site : {} \n'.format(tswork.stat))
    if 'calc_center' in list(tswork.anex.keys()):
        outfile.write('# Analysis Centre: {} \n'.format(tswork.anex['calc_center']))
    else:
        outfile.write('# Analysis Centre: N/A \n')

    outfile.write('# Solution code: GINS_PS \n')
    outfile.write('# Datum: ITRF2008\n')
    outfile.write('#\n')
    outfile.write('# Reference epoch: {}\n'.format(geok.toYearFraction(first_pt.Tdt)))
    outfile.write('# X : {}\n'.format(first_pt.X))
    outfile.write('# Y : {}\n'.format(first_pt.Y))
    outfile.write('# Z : {}\n'.format(first_pt.Z))
    outfile.write('#\n')
    outfile.write('# Longitude : {}\n'.format(first_pt.L))
    outfile.write('# Latitude  : {}\n'.format(first_pt.F))
    outfile.write('# Height    : {}\n'.format(first_pt.H))
    outfile.write('#\n')
    if tswork.bool_discont:
        outfile.write('# type_of_offset : from discontinuties got from a station.info\n')
        outfile.write('#\n')
        for disc in sorted(tswork.discont):
            outfile.write('# offset {} 7\n'.format(geok.toYearFraction(disc)))
        outfile.write('#\n')
    # write the data
    outfile.write('#  Year         DN           DE           DH        SDN       SDE       SDH\n')
    for e,n,u,t,se,sn,su in zip(E,N,U,T,sE,sN,sU):
        t = geok.toYearFraction(geok.posix2dt(t))
        outfile.write('{:.5f}   {:+.6f}    {:+.6f}    {:+.6f} {:+.6f} {:+.6f} {:+.6f}\n'.format(t,n-n0,e-e0,u-u0,se,sn,su))

    print('INFO : timeserie exported in ' + outpath)
    return None


def export_ts_as_hector_enu(tsin,outdir,outprefix,coordtype = 'ENU'):
    """
    export to a HECTOR .enu (and not .neu !) compatible format
    This format is simpler : just gives MJD E N U

    This format is necessary to force a sampling period.

    outfile will be writed in
    /outdir/outprefixSTAT.enu
    """

    print("NOT IMPLEMENTED YET !")


    return None



def export_ts_as_midas_tenu(tsin,outdir,outprefix,coordtype = 'ENU',
                            export_step=True):
    """
    export to a MIDAS .tneu compatible format

    outfile will be writed in
    /outdir/outprefixSTAT.tneu

    if export_step == True:
    export a step file as
    /outdir/outprefixSTAT.step
    """
    if not hasattr(tsin[0],'X'):
        print('WARN : export_ts_as_midas_tneu : no XYZ in ts')
        noXYZ = True
    else:
        noXYZ = False

    tswork = copy.deepcopy(tsin)
    stat = tswork.stat
    if coordtype == 'XYZ':
        mp = tswork.mean_posi()
        tswork.ENUcalc(mp)
    outpath = outdir +'/' + outprefix + tswork.stat + '.tenu'
    outfile = open(outpath,'w+')
    E,N,U,T,sE,sN,sU = tswork.to_list('ENU')
    first_pt = tswork.pts[0]
    if noXYZ:
        first_pt.X = 0.
        first_pt.Y = 0.
        first_pt.Z = 0.
        first_pt.L = 0.
        first_pt.H = 0.
        first_pt.F = 0.

    e0,n0,u0,t0 = list(zip(E,N,U,T))[0]

    for e,n,u,t in zip(E,N,U,T):
        t = geok.toYearFraction(geok.posix2dt(t))
        #outfile.write('{} {:.5f} {:+.6f} {:+.6f} {:+.6f} \n'.format(stat,t,n-n0,e-e0,u-u0))
        outfile.write('{} {:.5f} {:+.6f} {:+.6f} {:+.6f} \n'.format(stat,t,e-e0,n-n0,u-u0))

    print('INFO : timeserie exported in ' + outpath)

    if export_step and tswork.bool_discont:
        outpath_step = outdir +'/' + outprefix + tswork.stat + '.step'
        outfile_step = open(outpath_step,'w+')
        for d in tswork.discont:
            d = geok.toYearFraction(d)
            line = tswork.stat + " " + str(d) + "\n"
            outfile_step.write(line)

            print('INFO : timeserie discont. (steps) exported in ' + outpath_step)

    return None




def decimate_cleaner(tsin,minval,in_place = False):
    """
    in_place DOES'NT WORK !!!
    """
    if not in_place:
        tsout = copy.deepcopy(tsin)
    else:
        tsout = tsin
    if minval == 0:
        return ts_out

    tsout.del_data()
    for pt in tsin.pts:
        #print pt.T ,minval , np.mod(pt.T,minval)
        if np.mod(np.round(pt.T,2),minval) < 10**-3:
            tsout.add_point(pt)
    tsout.interval_nominal()
    return tsout

def decimate_cleaner_2(tsin,N,in_place = False):
    """ keep a value every N vals """
    if in_place:
        tsout = copy.deepcopy(tsin)
    else:
        tsout = tsin
    tsout.del_data()
    tsout.pts = tsin.pts[:-1:N]
    tsout.interval_nominal()
    return tsout


def mean_list_of_pts(ptslisin):
    """ useful for merge fct
    ONLY IMPLEMENTED FOR ENU coords for the moment """
    E = np.nanmean([p.E for p in ptslisin])
    N = np.nanmean([p.N for p in ptslisin])
    U = np.nanmean([p.U for p in ptslisin])
    sE = np.nanmean([p.sE for p in ptslisin])
    sN = np.nanmean([p.sN for p in ptslisin])
    sU = np.nanmean([p.sU for p in ptslisin])

    Tlis = [p.T for p in ptslisin]
    T = ((np.max(Tlis) - np.min(Tlis)) / 2.) + np.min(Tlis)
    pt = Point(E,N,U,T,'ENU',sE,sN,sU,ptslisin[0].name)

    return pt


def merge(tsin,N):
    """ merge N points in one """
    tsin.sort()
    tsout = copy.deepcopy(tsin)
    tsout.del_data()
    slic = genefun.sliceIt(tsin.pts,N)
    for sl in slic:
        pt = mean_list_of_pts(sl)
        tsout.add_point(pt)
    tsout.interval_nominal()
    return tsout

def merge_ts(ts_list_in):
    pts_list_merged = []
    for ts in ts_list_in:
        pts_list_merged = pts_list_merged + ts.pts

    ts_out = TimeSeriePoint()
    ts_out.pts = pts_list_merged
    ts_out.anex = ts_list_in[0].anex
    ts_out.interval_nominal()
    ts_out.sort()
    return ts_out


def time_win(tsin,windows,mode='keep',outbool=False):
    if type(windows[0][0]) == dt.datetime:
        Ttype = 'Tdt'
    else:
        Ttype = 'T'

    tsout = copy.copy(tsin)
    tsout.del_data()

    if isinstance(tsout,TimeSeriePoint):
        data = 'pts'
    elif isinstance(tsout,TimeSerieObs):
        data = 'obs'
    else:
        print('BUG : time_win : ca chie')

    outboollist = []

    print(tsin)
    for i,e in enumerate(windows):
        print('windows', i)
        print(e[0] , e[1])


    for pt in getattr(tsin,data):
        boollist = []
        for win in windows:
            if win[1] < win[0]:
                print("WARN : time_win : le début de periode est + grand que la fin")
                print(win[0] , win[1])
                print(mode)

            if win[0] < getattr(pt,Ttype) and win[1] > getattr(pt,Ttype):
                boollist.append(True)
            else:
                boollist.append(False)

        if mode == 'keep' :
            finalbool = np.any(boollist)
        elif mode == 'del':
            finalbool = not np.any(boollist)

        if finalbool:
            tsout.add_point(pt)

        outboollist.append(finalbool)

    tsout.bool_interp_uptodate = False

    #tsout.interp_set()
    if outbool:
        outboollist = np.array(outboollist)
        return tsout,outboollist
    else:
        return tsout

def time_win_T(Tin,win,mode='del'):

    boolfinalist = []

    for t in Tin:
        booltemp = []
        for w in win:
            startw = w[0]
            endw = w[1]

            # CONVENTION : VRAI SI DANS L'INTERVALLE
            # QU'IMPORTE SI KEEP OU DEL
            if startw < t and t < endw :
                booltemp.append(True)
            else:
                booltemp.append(False)

        if mode == 'keep' :
            finalbool = np.any(booltemp)
        elif mode == 'del':
            finalbool = not np.any(booltemp)

        boolfinalist.append(finalbool)

    return np.array(boolfinalist)

def time_gap(tsin,marge=2,mode='del'):
    ''' ENTREE une TimeSerie
        SORTIE une window (liste de listes)'''

    inomi = tsin.i_nomi

    T = tsin.to_list()[3]

    Tdiff = np.diff(np.sort(T))

    boollist = []
    # Dans boollist :
    # False si c'est gap
    # True sinon

    for td in Tdiff:
        if td > inomi * marge:
            boollist.append(False)
        else:
            boollist.append(True)

    # Difference suivant le mode
    if mode == 'keep' :
        modbool = True
    elif mode == 'del':
        modbool = False

    indbool = np.where(np.array(boollist) == modbool)[0]
    protowin = geok.group_consecutives(indbool)

    winout = []

    for e in protowin:
        if len(e) == 1:
            eout = [T[e[0]],T[e[0]+1]]
        elif len(e) == 2:
            eout = [T[e[0]],T[e[1]+1]]
        else:
            "BUG : time_gap"
        winout.append(eout)

    return winout

def compar_elts_in_ts(ts1,ts2):
    ''' ts2 must contains less elts than ts1 (ts2 = cleaned one)'''
    T1 = ts1.to_list()[3]
    T2 = ts1.to_list()[3]

    for t2 in T2:
        T1.remove(t2)

    return T1


def rotate_pt_cls_solo(tsattin,pointin,Rtype='R1', xyzreftuple = ([1, 0, 0], [0, 1, 0], [0, 0, 1]),angtype='deg'):
    ''' ENTREE : tsattin : une TS d'attitude  (N angles)
                 pointin : UN Point en entrée

        SORTIE : une TSpoint de N points '''

    attlisttmp = tsattin.to_list()
    R = attlisttmp[0]
    P = attlisttmp[1]
    Y = attlisttmp[2]
    T = attlisttmp[3]

    A,B,C,_,_ = pointin()

    ptrotlis = geok.rotate_points(R,P,Y,[np.array([A,B,C])],Rtype,xyzreftuple,angtype)

    ptrotlis = genefun.shrink_listoflist(ptrotlis)

    tsout = TimeSeriePoint()

    for p,t in zip(ptrotlis,T):
        #print p
        tsout.add_point(Point(p[0],p[1],p[2],t,initype=pointin.initype))

    return tsout


def rotate_points_class(tsattin,ptslin, Rtype='R1', xyzreftuple = ([1, 0, 0], [0, 1, 0], [0, 0, 1]),angtype='deg'):

    listsout = []

    for pt in ptslin:
        ts = rotate_pt_cls_solo(tsattin,pt,Rtype,xyzreftuple,angtype)
        listsout.append(ts)

    return listsout

def add_offset_point(ptin,dA,dB,dC,coortype='ENU'):
    """
    ONLY IMPLEMENTED FOR ENU FOR THE MOMENT
    150415 : remark still necessary ???

    coortype == 'UXYZ' :
        specific case where we correct an Up offset directly in the XYZ coords
        very usefull for an antenna offset in for a moving GPS
        (but works only for the up)
    """
    ptout = copy.copy(ptin)
    if coortype == 'ENU':
        ptout.ENUset(ptin.E + dA,ptin.N + dB,ptin.U + dC,
                     ptin.sE,ptin.sN,ptin.sU)
    elif coortype == 'XYZ':
        ptout.XYZset(ptin.X + dA,ptin.Y + dB,ptin.Z + dC,
                     ptin.sX,ptin.sY,ptin.sZ)
    elif coortype == 'FLH':
        ptout.FLHset(ptin.F + dA,ptin.L + dB,ptin.H + dC,
                     ptin.sF,ptin.sL,ptin.sH)

    elif coortype == 'UXYZ':
        N = geok.normal_vector(ptin.F , ptin.L , ptin.H)
        dA2 , dB2 , dC2 = N * dC
        ptout.XYZset(ptin.X + dA2,ptin.Y + dB2,ptin.Z + dC2,
                     ptin.sX,ptin.sY,ptin.sZ)

    return ptout

def add_offset_ts(tsin,dA,dB,dC,coortype='ENU'):
    """
    return a copy of the tsin
    (tsin won't be affected)

    coortype == 'UXYZ' :
        specific case where we correct an Up offset directly in the XYZ coords
        very usefull for an antenna offset in for a moving GPS
        (but works only for the up)
    """
    tsout = copy.deepcopy(tsin)
    tsout.del_data()
    newpts = [add_offset_point(pt,dA,dB,dC,coortype) for pt in tsin.pts]
    for pt in newpts:
        tsout.add_point(pt)
    return tsout

def add_offset_smart_for_GINS_kine(tsin,tslist_offset_3ple,list_windows,coortype='XYZ'):
    """
    tslist_offset_3ple : list of len N containing (dX,dY,dZ) offsets
    list_windows       : list of len N-1 containing dates of changes
    """

    if type(list_windows[0]) is int:
        list_windows = [geok.posix2dt(e) for e in list_windows]

    list_windows = [dt.datetime(1980,1,1)] + list_windows + [dt.datetime(2099,1,1)]
    tsout = copy.deepcopy(tsin)
    tsout.del_data()
    olddX ,olddY , olddZ = 12345 , 12345 , 12345
    for pt in tsin:
        for i in range(len(list_windows)-1):
            if list_windows[i] <= pt.Tdt < list_windows[i+1]:
                curdX , curdY , curdZ = tslist_offset_3ple[i][0],tslist_offset_3ple[i][1],tslist_offset_3ple[i][2]
                if curdX != olddX or curdY != olddY or curdZ != olddZ:
                    olddX ,olddY , olddZ = curdX , curdY , curdZ
                    print('INFO : change of time window / offset')
                ptout = add_offset_point(pt,curdX,curdY,curdZ,coortype='XYZ')
                tsout.add_point(ptout)
                continue
    return tsout


def find_pts_from_ts_with_time(tin,tstupin,tol=0.001):
    ptsout = []
    for ts in tstupin:
        pt, i = ts.find_point(tin,tol=tol)
        if np.isnan(i):
            print("WARN : find_pts_from_ts_with_time : no pt find")
        ptsout.append(pt)
    return ptsout

def dist_btwn_2pts(ptA,ptB,coortype='XYZ'):

        A = np.array([ptA.X,ptA.Y,ptA.Z])
        B = np.array([ptB.X,ptB.Y,ptB.Z])

        return geok.dist(A,B)

def dist_diff_btwn_2pts(ptA,ptB):
    #dérive la distance entre 2 points A et B
        A = np.array([ptA.X,ptA.Y,ptA.Z])
        B = np.array([ptB.X,ptB.Y,ptB.Z])

        dAB   = A-B
        dist  = np.linalg.norm(dAB)

        diffA =   dAB / dist
        diffB = - dAB / dist

        return diffA, diffB

def helmert_trans(tsin,params='itrf2008_2_etrf2000', invert=False):
    tsout = copy.deepcopy(tsin)
    [pt.helmert_trans(params,invert) for pt in tsout.pts]
    return tsout


def velocity_trans(tsin , vx, vy, vz, epoc_init = 'auto' , epoc_end = 'auto'):
    tsout = copy.deepcopy(tsin)
    [pt.velocity_trans(vx, vy, vz, epoc_init , epoc_end) for pt in tsout.pts]
    return tsout


def interpolator_light(T,X,Y,Z):
    from scipy.interpolate import interp1d

    IX = interp1d(T,X,bounds_error=False,fill_value=np.nan)
    IY = interp1d(T,Y,bounds_error=False,fill_value=np.nan)
    IZ = interp1d(T,Z,bounds_error=False,fill_value=np.nan)

    return IX , IY , IZ

def interpolator_with_extrapolated(T,X,Y,Z):
    from scipy.interpolate import interp1d

    IX = interp1d(T,X,bounds_error=False,fill_value="extrapolate")
    IY = interp1d(T,Y,bounds_error=False,fill_value="extrapolate")
    IZ = interp1d(T,Z,bounds_error=False,fill_value="extrapolate")

    return IX , IY , IZ




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


class point_n_click_plot():
    """
    USAGE :
        Data have to be ploted already in a figure

        Then :

            PnC = point_n_click_plot()

            multi , cid = PnC(fig=6)

            PnC.selectedX

        i.e. :

            Create an object point_n_click_plot (here it is PnC in the exemple below)

            Call the object like a function with the id of the plot figure

            Make your selection using SPACE key

            Get your results in a list called PnC.selectedX


    IMPORTANT : cursor objects (i.e. multi & cid)
                must be stored as global variables when you call the method
                like this :

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

        print("INFO : press SPACE to record a X-value, \n       press R to Remove the previously recorded one")


        def onclick_discont(event):

            if event.key == ' ':
                if Xdata_are_time:
                    ix, iy = matplotlib.dates.num2date(event.xdata).replace(tzinfo=None),event.ydata
                else:
                    ix, iy = event.xdata , event.ydata

                print("INFO : X value recorded : " , ix)

                for ax in figobj.axes:
                    out_bar_list = geok.plot_vertical_bar_ax([ix],ax,"b",
                                                             linewidth=1)

                self.ver_bar_stk.append(out_bar_list[0])

                self.selectedX.append(ix)

                plt.draw()

            elif event.key == 'r' and len(self.selectedX) > 0:
                last = self.selectedX[-1]

                self.selectedX.remove(last)
                #for ax in figobj.axes:
                #    geok.plot_vertical_bar_ax([last],ax,"r")
                last_bar = self.ver_bar_stk[-1]
                self.ver_bar_stk.remove(last_bar)
                last_bar.remove()
                print("INFO : value removed : " , last )

                plt.draw()

            return None

        multi = MultiCursor(figobj.canvas, figobj.axes , color='k', lw=1)
        cid   = figobj.canvas.mpl_connect('key_press_event', onclick_discont)

        return multi , cid
