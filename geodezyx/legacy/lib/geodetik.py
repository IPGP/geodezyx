# -*- coding: utf-8 -*-
"""
Created on Fri May 30 12:50:26 2014

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

import genefun
import softs_runner
import transformations as trans

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
import re
import os
import multiprocessing as mp
import itertools


#if platform.node() in ('calipso' , "diamant"):
    #import cgkitmod.cgtypes as cgt
    ##import cgtypes as cgt

import sys
sys.dont_write_bytecode = True


############################################################################
### Coordinates and angle conversion, rotation matrices, and projections
############################################################################
### Import conversion_coords into geodetik for legacy reasons
from conversion_coords import *
############################################################################

############################################################################
### Time conversion
############################################################################
### Import conversion_time into geodetik for legacy reasons
from conversion_time import *
############################################################################

############################################################################
### Statistics functions
############################################################################
### Import geok_statistics into geodetik for legacy reasons
from geok_statistics import *
############################################################################

############################################################################
### Least Squares functions
############################################################################
### Import geok_least_squares into geodetik for legacy reasons
from geok_least_squares import *
############################################################################

############################################################################
### Euler pole / tectonic plate velocity determination
############################################################################
### Import geok_euler_pole_calc into geodetik for legacy reasons
from geok_euler_pole_calc import *
############################################################################

############################################################################
### Quaternion management and interpolation
############################################################################
### Import geok_quaternions into geodetik for legacy reasons
from geok_quaternions import *
############################################################################





def BL_from_points(listpointin):
    ''' 
    A partir d'une liste de points,
    retourne les baselines entre ces points dans une matrice 
    '''
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
    BL = np.empty((N,N))

    for i,pt1 in enumerate(listpointin):
        for j,pt2 in enumerate(listpointin):

            if i == j:
                BL[i,j] = 0
            else:
                BL[i,j] = np.linalg.norm(pt1 - pt2)

    return BL



def mat_poids(Sinp,Ninp,fuvinp=1):
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
        print(Sinp[i])
        Ktemp.append(np.eye(Ninp[i]) * Sinp[i]**2)

    K = scipy.linalg.block_diag(*Ktemp)
    Q = (1/fuvinp) * K
    P = scipy.linalg.inv(Q)

    return K , Q , P

def rotmat2(theta,angtype='deg'):

    if angtype == 'deg':
        theta = np.deg2rad(theta)

    rotmat = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])

    return rotmat


def rotmat3(alpha,beta,gamma,xyzreftuple = ([1, 0, 0], [0, 1, 0], [0, 0, 1]),angtype='deg'):

    xaxis, yaxis, zaxis = xyzreftuple

    if angtype == 'deg':
        alpha = np.deg2rad(alpha)
        beta  = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)

    Rx = trans.rotation_matrix(alpha, xaxis)
    Ry = trans.rotation_matrix(beta, yaxis)
    Rz = trans.rotation_matrix(gamma, zaxis)
    R = trans.concatenate_matrices(Rz, Ry, Rx)[:3,:3]

    return R

def rotate_points(alphal,betal,gammal,pointlin,Rtype='R1',
                  xyzreftuple = ([1, 0, 0], [0, 1, 0], [0, 0, 1]),
                  angtype='deg',fullout = False):
    '''
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
        liste de listes de P * [ points ]

    Sortie :
        liste de listes [ [ xA ] [ xA ] ... xP [ xA ] ]  '''

    xaxis, yaxis, zaxis = xyzreftuple

    if not genefun.is_iterable(alphal):
        alphal = np.array([alphal])
        betal = np.array([betal])
        gammal = np.array([gammal])
        boolnotiterable = True
    else:
        boolnotiterable = False

    pointlout = []
    R_out = []


    for pt in pointlin:

        if not genefun.is_iterable(pt) or len(pt) != 3:
            print("ERR : rotate_points : pts != 3 coords")
            return 0

        pointltmp = []

        for a,b,g in zip(alphal,betal,gammal):

            R1 = rotmat3(a,b,g,angtype=angtype,xyzreftuple=xyzreftuple)
            R2 = C_rpy2enu(a,b,g,angtype=angtype)

            if Rtype == 'R1':
                R = R1
            elif Rtype == 'R2':
                R = R2
            R_out.append(R)

            pointltmp.append(np.dot(R,pt))

        pointlout.append(pointltmp)

        if boolnotiterable:
            pointlout = pointltmp

        pointlout = np.array(pointlout)

    if fullout:
        return pointlout , R_out
    else:
        return pointlout


def guess_seq_len(seq):
    #source
    #http://stackoverflow.com/questions/11385718/python-finding-repeating-sequence-in-list-of-integers
    guess = 1

    if len(set(seq)) == 1:
        return 1

    max_len = len(seq) / 2
    for x in range(2, max_len):
        if seq[0:x] == seq[x:2*x] :
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
    outlon = np.mod(lon , 2*np.pi)
    outlon[np.logical_and(outlon == 0 , positiv)] = 2 * np.pi
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
    q =  np.logical_and((outlon < -np.pi) , (np.pi < outlon))
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

    positiveInput = (lon > 0)
    lon = np.mod(lon, 360)
    lon[(lon == 0) & positiveInput] = 360
    return lon



# Pas convaincu de son utilité
def unwrap180(anglist,angtype='deg'):

    if angtype == 'deg':
        seuil = 360

    angout = []

    for a in anglist:
        if a > seuil / 2:
            a = a - seuil
        angout.append(a)

    return angout


def wrap360(anglist,angtype='deg'):

    angout = []

    if angtype == 'deg':
        seuil = 360
    elif angtype == 'rad':
        seuil = 2*np.pi

    for a in anglist:
        if a < 0:
            a = a + seuil

        angout.append(a)

    return angout

class interp1d_ang():

    def __init__(self,T,A,angtype='deg',kind='linear',bounds_error=False):

        if angtype == 'deg':
            A = np.deg2rad(A)

        self.A = A
        self.T = T
        self.C = np.cos(A)
        self.S = np.sin(A)

        self.CfT = interpolate.interp1d(T,self.C,kind=kind,bounds_error=bounds_error)
        self.SfT = interpolate.interp1d(T,self.S,kind=kind,bounds_error=bounds_error)


    def __call__(self,T,angtype='deg'):

        I = np.arctan2(self.SfT(T) ,self.CfT(T) )
        I = wrap360(I,angtype='rad')

        if angtype == 'deg':
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
                result2.append([r[0],r[-1]])
            else:
                result2.append(r)

    return result2


def rinex_lister(path,add_new_names=True):
    """
    find all rinex in a folder and his subfolders
    path can be a string or a tuple of string => manage multi files :)

    is very similar with softs_runner.multi_finder_rinex and
    gins_runner.get_rinex_list
    """

    if type(path) is str:
        path = [path]

    paths_walked = []
    for p in path:
        paths_walked = paths_walked + list(os.walk(p))

    wholefilelist = []
    for tup in paths_walked:
        wholefilelist = wholefilelist + tup[-1]

    wholefilelist = list(set(wholefilelist))
    
    rinexfilelist           = [fil for fil in wholefilelist if re.search( softs_runner.rinex_regex() , os.path.basename(fil))]
    if add_new_names:
        rinexfilelist_new_names = [fil for fil in wholefilelist if re.search( softs_runner.rinex_regex_new_name() , os.path.basename(fil))]
        rinexfilelist = rinexfilelist + rinexfilelist_new_names
    
    print('INFO : ' , len(rinexfilelist) , 'RINEXs found')

    return rinexfilelist



def rinex_timeline(inputlist_or_paths,start = dt.datetime(1980,1,1) ,
                   end = dt.datetime(2099,1,1),use_rinex_lister = True,
                   dots_plot=False,jul_date_plot=False,return_figure=False):
    """
    if use_rinex_lister = True :
    inputlist_or_paths is a path OR a list of path
    where RINEX can be found (use the subfunction rinex_lister for that)
    else if use_rinex_lister = False:
    it's a list of RINEXs
    """

#    if use_rinex_lister:
#        filelist = rinex_lister(inputlist_or_paths)
#    else:
#        filelist = inputlist_or_paths
#
#    rinexfilelist = [fil for fil in filelist if re.match( '.*' + softs_runner.rinex_regex() + '$', fil)]
#
#    if not use_rinex_lister:
#        rinexfilelist = [os.path.basename(e) for e in rinexfilelist]
#
#    print('INFO : ', len(rinexfilelist), 'RINEXs will be ploted on the timeline')
#
#    statname_lis = sorted(list(set([rin[0:4] for rin in rinexfilelist])))
#
#    print('INFO : ', len(statname_lis), 'stations will be ploted on the timeline')
#
#    datadico = dict()
#
#    for stat in statname_lis:
#        datadico[stat] = []
#
#    for rnx in rinexfilelist:
#        try:
#            datadico[rnx[0:4]].append((rnx,rinexname2dt(rnx)))
#        except:
#            print('error with : ', rnx)
#

    datadico = rinex_timeline_datadico(inputlist_or_paths,
                                       use_rinex_lister = use_rinex_lister)

    fig  = timeline_plotter(datadico,start = start ,
                     end=end ,dots_plot=dots_plot,
                     jul_date_plot=jul_date_plot)

#    ax.xaxis.grid(True)
#    plotstat_lis = []
#    for i,stat in enumerate(reversed(sorted(datadico.keys()))):
#        print(stat)
#        T = [ e[-1] for e in datadico[stat] ]
#        T = [ t for t in T if start <= t <= end ]
#        T = sorted(T)
#
#        plotstat_lis.append(stat)
#        if dots_plot:
#            #old old with dots
#            ax.plot(T,i*np.ones(len(T)), '.')
#        else:
#            TGrp = genefun.consecutive_groupIt(dt2MJD(T),True)
#            for tgrp in TGrp:
#                if not jul_date_plot:
#                    tgrp = MJD2dt(tgrp)
#                if tgrp[0] == tgrp[1]:
#                    ax.plot(tgrp[0],i, '.')
#                else:
#                    ax.plot(tgrp,[i]*2, '-')
#
##    ax.set_yticks(np.arange(0,len(plotstat_lis)-1),plotstat_lis)
#    plt.yticks(np.arange(0,len(plotstat_lis)),plotstat_lis)
#    fig.autofmt_xdate()
#    fig.set_size_inches(11.69,i * 0.28) #16.53
#    ax.set_ylim([-1 , len(plotstat_lis) + 1])
##    return fig
    if not return_figure:
        return datadico
    else:
        return datadico , fig

def rinex_timeline_datadico(inputlist_or_paths,use_rinex_lister = True,
                            optional_info=''):
    """
    convention for RINEX datadico :
        datadico[stat] = [(rinexname1,optional1,date1) ... (rinexnameN,optionalN,dateN)]
    """
    if use_rinex_lister:
        filelist = rinex_lister(inputlist_or_paths)
    else:
        filelist = inputlist_or_paths

    #rinexfilelist = [fil for fil in filelist if re.search( '.*' + softs_runner.rinex_regex() + '$', fil)]
    rinexfilelist_old = [fil for fil in filelist if re.search( softs_runner.rinex_regex() , fil)]
    rinexfilelist_new = [fil for fil in filelist if re.search( softs_runner.rinex_regex_new_name() , fil)]

    rinexfilelist = rinexfilelist_old + rinexfilelist_new

    if not use_rinex_lister:
        rinexfilelist = [os.path.basename(e) for e in rinexfilelist]

    print('INFO : ', len(rinexfilelist), 'RINEXs will be ploted on the timeline')

    statname_lis = sorted(list(set([rin[0:4] for rin in rinexfilelist])))

    print('INFO : ', len(statname_lis), 'stations will be ploted on the timeline')

    datadico = dict()

    for stat in statname_lis:
        datadico[stat] = []

    for rnx in rinexfilelist:
        try:
            datadico[rnx[0:4]].append((rnx,optional_info,rinexname2dt(rnx)))
        except:
            print('error with : ', rnx)

    return datadico

def rinex_timeline_datadico_merge_not_very_smart(datadico_list,priority_list):
    """
    Merge different RINEXs datadico, produced by rinex_timeline_datadico
    coming from different archives
    Args :
        rinex_timeline_datadico : list of RINEX datadico
        priority_list : priority list of 'optional_info' (archive ID)
                        it will erase optional_info of lower priority
    Returns :
        datadico_out : a merged datadico
    """

    datadico_out  = dict()

    datadico_merged = genefun.dicts_of_list_merge(*datadico_list)

    for k , dataval in datadico_merged.items():

        rnxname_list = [e[0]  for e in dataval]
        archive_list = [e[1]  for e in dataval]
        date_list    = [e[-1] for e in dataval]

        out_date_list , out_all_list = [] , []
        for r,a,d in zip(rnxname_list,archive_list,date_list):
            if d not in out_date_list:
                out_date_list.append(d)
                out_all_list.append((r,a,d))
            else:
                ind_existing   = out_date_list.index(d)
                archd_existing = out_all_list[ind_existing][1]
                if priority_list.index(a) < priority_list.index(archd_existing):
                    out_date_list.remove(d)
                    out_all_list.remove(out_all_list[ind_existing])
                    out_date_list.append(d)
                    out_all_list.append((r,a,d))

        datadico_out[k] = out_all_list

    return datadico_out


def rinex_timeline_datadico_merge(datadico_list,priority_list=None):
    """
    Merge different RINEXs datadico, produced by rinex_timeline_datadico
    coming from different archives
    Args :
        rinex_timeline_datadico : list of RINEX datadico
        priority_list : priority list of 'optional_info' (archive ID)
                        it will erase optional_info of lower priority
                        NB : it is not very useful, just sort
                             datadico_list in the right order ...
    Returns :
        datadico_out : a merged datadico
    """

    datadico_out  = dict()

    datadico_merged = genefun.dicts_of_list_merge(*datadico_list)

    for k , dataval in datadico_merged.items():

        rnxname_list = [e[0]  for e in dataval]
        archive_list = [e[1]  for e in dataval]
        date_list    = [e[-1] for e in dataval]

        out_date_list , out_all_list = [] , []
        for r,a,d in zip(rnxname_list,archive_list,date_list):
            if d not in out_date_list:
                out_date_list.append(d)
                out_all_list.append((r,a,d))
            elif priority_list:
                ind_existing   = out_date_list.index(d)
                archd_existing = out_all_list[ind_existing][1]
                if priority_list.index(a) < priority_list.index(archd_existing):
                    out_date_list.remove(d)
                    out_all_list.remove(out_all_list[ind_existing])
                    out_date_list.append(d)
                    out_all_list.append((r,a,d))
            else:
                pass

        datadico_out[k] = out_all_list

    return datadico_out

def timeline_plotter(datadico,start = dt.datetime(1980,1,1) ,
                     end = dt.datetime(2099,1,1),dots_plot=False,
                     jul_date_plot=False,datadico_anex_list = [],
                     use_only_stats_of_main_datadico=False,
                     colordico_for_main_datadico=None):
    """
    A simpler version has been commited to geodezyx toolbox for archive
    on 20180118 15:59A
    """

    fig , ax = plt.subplots()
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)

    if not use_only_stats_of_main_datadico:
        stats_concat_list = list(datadico.keys()) + sum([list(e.keys()) for e in datadico_anex_list], [])
        stats_concat_list = list(reversed(sorted(list(set(stats_concat_list)))))
    else:
        stats_concat_list = list(reversed(sorted(list(datadico.keys()))))

    # the plot has not the same behavior if it is the morning or not 
    # (rinexs timelines wont be ploted if it is the morning)
    if dt.datetime.now().hour < 12:
        morning_shift = dt.timedelta(days=1)
    else:
        morning_shift = dt.timedelta(days=0)

    legend_list = [] # must be here before the loop
    for i,stat in enumerate(stats_concat_list):
        # PART 1 : PLOT MAIN DATADICO
        if not stat in datadico.keys():
            continue

        #T = Time, O = Station name (Observation)
        Torig = [ e[-1] for e in datadico[stat] ]
        Oorig = [ e[1]  for e in datadico[stat] ]

        # Time windowing
        T , O = [] , []
        for t , o in zip(Torig , Oorig):
            if ( start - morning_shift ) <= t <= end: 
                T.append(t)
                O.append(o)

        T,O = genefun.sort_binom_list(T,O)

        TMJD=dt2MJD(T)
        TGrpAll = genefun.consecutive_groupIt(TMJD,True) # Tuples (start,end) of continue period

        for tgrp in TGrpAll:
            color1 = 'C0'
            color2 = ''
            extra_archive = False

            ###*** managing colors
            if colordico_for_main_datadico:
                igrpstart    = TMJD.index(tgrp[0])
                igrpend      = TMJD.index(tgrp[1])
                Ogrp         = O[igrpstart:igrpend+1]
                Tgrp         = TMJD[igrpstart:igrpend+1] # all dates in the current continue period

                opt_set      = list(set(Ogrp))

                if len(opt_set) == 1: # Regular case : only one archive
                    if opt_set[0] in colordico_for_main_datadico:
                        color1 = colordico_for_main_datadico[opt_set[0]]
                else: # several archives ... so the line has to be splited in segments
                    extra_archive = True
                    OSubgrp = genefun.identical_groupIt(Ogrp)
                    TSubgrp = genefun.sublistsIt(Tgrp , [len(e) for e in OSubgrp])

                    Tgrp_plt      = [ (e[0] , e[-1] + 1)    for e in TSubgrp ] # +1 because the end boundary day is not included
                    Ogrp_plt      = [ list(set(e))[0] for e in OSubgrp     ]
                    Color_grp_plt = [ colordico_for_main_datadico[e] if e in colordico_for_main_datadico else color2 for e in Ogrp_plt ]
            ###*** End of managing colors

            tgrp = (tgrp[0] , tgrp[1] + 1) # +1 because the end boundary day is not included
                                           # must stay there, in case of not colordico_for_main_datadico

            if not jul_date_plot:
                tgrp = MJD2dt(tgrp)
                if extra_archive:
                    Tgrp_plt = [ (MJD2dt(e[0]) , MJD2dt(e[1])) for e in Tgrp_plt ]

            #PLOT part
            if tgrp[0] == tgrp[1] + dt.timedelta(days=1): # CASE NO PERIOD, ONLY ONE DAY
                ax.plot(tgrp[0],i , '.' + color1)
            else:                  # CASE PERIOD
                if not extra_archive: # ONE ARCHIVE REGULAR CASE
                    ax.plot(tgrp,[i]*2, '-' + color1)
                else:              # SUBCASE "EXTRA" : SEVERAL ARCHIVE
                    for tt , cc in zip(Tgrp_plt , Color_grp_plt):
                        ax.plot(tt,[i]*2 , '-' + cc)

        # PART 2 : PLOT ANEX DATADICO
        for idatadico_anex,datadico_anex in enumerate(datadico_anex_list):
            if stat in list(datadico_anex.keys()):
                T = datadico_anex[stat]
                T = [ t for t in T if start <= t <= end  ]
                T = sorted(T)

                if jul_date_plot:
                    T = MJD2dt(T)

                pale_blue_dot , = ax.plot(T,i*np.ones(len(T)), 'o', color='skyblue',label="final SNX")
                legend_list     = [pale_blue_dot]

    #### LEGEND
    if colordico_for_main_datadico:
        import matplotlib.lines as mlines
        for arcnam , col in colordico_for_main_datadico.items():
            legend_list.append(mlines.Line2D([], [], color=col,
                                                     label=arcnam))
            plt.legend(handles=legend_list,loc='upper left',ncol=3,
                       columnspacing=1)


#    ax.set_yticks(np.arange(0,len(plotstat_lis)-1),plotstat_lis)
    plt.yticks(np.arange(0,len(stats_concat_list)),stats_concat_list)
    fig.autofmt_xdate()
    koef = np.sqrt(2) * 1
    fig.set_size_inches(koef*11.69,koef*len(stats_concat_list) * 0.28) #16.53
    ax.set_ylim([-1 , len(stats_concat_list) + 1])

    return fig


# def timeline_plotter_old_bkp(datadico,start = dt.datetime(1980,1,1) ,
#                      end = dt.datetime(2099,1,1),dots_plot=False,
#                      jul_date_plot=False,datadico_anex_list = [],
#                      use_only_stats_of_main_datadico=False,
#                      colordico_for_main_datadico=None):
#     """
#     A simpler version has been commited to geodezyx toolbox for archive
#     on 20180118 15:59A
#     """
#
#     fig , ax = plt.subplots()
#     ax.xaxis.grid(True)
#     ax.yaxis.grid(True)
#
#     if not use_only_stats_of_main_datadico:
#         stats_concat_list = list(datadico.keys()) + sum([list(e.keys()) for e in datadico_anex_list], [])
#         stats_concat_list = list(reversed(sorted(list(set(stats_concat_list)))))
#     else:
#         stats_concat_list = list(reversed(sorted(list(datadico.keys()))))
#
#     legend_list = [] # must be here before the loop
#     for i,stat in enumerate(stats_concat_list):
#         # PART 1 : PLOT MAIN DATADICO
#         if stat in datadico.keys():
#             T = [ e[-1] for e in datadico[stat] ]
#             T = [ t for t in T if start <= t <= end ]
#             T = sorted(T)
#
#             O = [ e[1] for e in datadico[stat] ]
#
#
#             if dots_plot:
#                 #old old with dots
#                 ax.plot(T,i*np.ones(len(T)), '.')
#             else:
#                 TMJD=dt2MJD(T)
#                 TGrpAll = genefun.consecutive_groupIt(TMJD,True)
#
#                 for tgrp in TGrpAll:
#
#                     color1 = ''
#                     color2 = ''
#                     extra_archive = False
#
#                     # managing colors
#                     if colordico_for_main_datadico:
#                         igrpstart    = TMJD.index(tgrp[0])
#                         igrpend      = TMJD.index(tgrp[1])
#                         Ogrp         = O[igrpstart:igrpend+1]
#                         Tgrp         = TMJD[igrpstart:igrpend+1]
#                         opt_set      = list(set(Ogrp))
#
#                         if len(opt_set) == 1: # Regular case : only one archive
#                             if opt_set[0] in colordico_for_main_datadico:
#                                 color1 = colordico_for_main_datadico[opt_set[0]]
#                         else: # several archives ... so the line has to be splited in segments
#                             extra_archive = True
#                             OSubgrp = genefun.identical_groupIt(Ogrp)
#                             TSubgrp = genefun.sublistsIt(Tgrp , [len(e) for e in OSubgrp])
#
#                             Tgrp_plt      = [ (e[0],e[-1] + 1)    for e in TSubgrp ] # +1 because the end boundary day is not included
#                             Ogrp_plt      = [ list(set(e))[0] for e in OSubgrp     ]
#                             Color_grp_plt = [ colordico_for_main_datadico[e] if e in colordico_for_main_datadico else color2 for e in Ogrp_plt ]
#
#                             print(Ogrp)
#
#                     tgrp = (tgrp[0] , tgrp[1] + 1) # +1 because the end boundary day is not included
#
#                     if not jul_date_plot:
#                         tgrp = MJD2dt(tgrp)
#                         if extra_archive:
#                             Tgrp_plt = [ (MJD2dt(e[0]) , MJD2dt(e[1])) for e in Tgrp_plt ]
#
#
#                     #PLOT part
#                     if tgrp[0] == tgrp[1]: # CASE NO PERIOD, ONLY ONE DAY
#                         ax.plot(tgrp[0],i , '.' + color1)
#                     else:                  # CASE PERIOD
#                         if not extra_archive: # ONE ARCHIVE REGULAR CASE
#                             ax.plot(tgrp,[i]*2, '-' + color1)
#                         else:              # SUBCASE "EXTRA" : SEVERAL ARCHIVE
#                             for tt , cc in zip(Tgrp_plt , Color_grp_plt):
#                                 print(i,stat,tt,cc)
#                                 ax.plot(tt,[i]*2 , '-' + cc)
#
#         # PART 2 : PLOT ANEX DATADICO
#         for idatadico_anex,datadico_anex in enumerate(datadico_anex_list):
#             if stat in list(datadico_anex.keys()):
#                 T = datadico_anex[stat]
#                 T = [ t for t in T if start <= t <= end  ]
#                 T = sorted(T)
#
#                 if jul_date_plot:
#                     T = MJD2dt(T)
#
#                 pale_blue_dot , = ax.plot(T,i*np.ones(len(T)), 'o', color='skyblue',label="final SNX")
#                 legend_list = [pale_blue_dot]
#
#     #### LEGEND
#     if colordico_for_main_datadico:
#         import matplotlib.lines as mlines
#         for arcnam , col in colordico_for_main_datadico.items():
#             legend_list.append(mlines.Line2D([], [], color=col,
#                                                      label=arcnam))
#             plt.legend(handles=legend_list,loc='upper left')
#
#
# #    ax.set_yticks(np.arange(0,len(plotstat_lis)-1),plotstat_lis)
#     plt.yticks(np.arange(0,len(stats_concat_list)),stats_concat_list)
#     fig.autofmt_xdate()
#     koef = np.sqrt(2) * 1
#     fig.set_size_inches(koef*11.69,koef*len(stats_concat_list) * 0.28) #16.53
#     ax.set_ylim([-1 , len(stats_concat_list) + 1])
#
#     return fig


# THIS FIRST SOLUTION PRINT DASHED LINES, NOT GOOD
#                        if len(opt_set) == 1: # Regular case : only one archive
#                            if opt_set[0] in colordico_for_main_datadico:
#                                color1 = colordico_for_main_datadico[opt_set[0]]
#                        else:
#                            dashed_line = True
#                            if opt_set[0] in colordico_for_main_datadico:
#                                color1 = colordico_for_main_datadico[opt_set[0]]
#                            if opt_set[1] in colordico_for_main_datadico:
#
# THIS SECOND SOLUTION PRINT X, NOT GOOD
#                            o_major = genefun.most_common(Ogrp)
#                            if o_major in colordico_for_main_datadico:
#                                color1 = colordico_for_main_datadico[o_major]
#
#                            Textra , Color_extra = [],[]
#                            Tgrp      = T[igrpstart:igrpend+1]
#
#                            for oo , tt in zip(Ogrp,Tgrp):
#                                if oo == o_major:
#                                    continue
#                                else:
#                                    if oo in colordico_for_main_datadico:
#                                        color2 = colordico_for_main_datadico[oo]
#                                    Textra.append(tt)
#                                    Color_extra.append(color2)


def listing_gins_timeline(path,stat_strt,stat_end,date_strt,date_end,suffix_regex=''):
    """ find all gins listings in a folder and his subfolders
        and plot timeline of the avaiable listings

        stat_strt,stat_end,date_strt,stat_end : where to find
        in the name the statname and the date"""

    fig , ax = plt.subplots()

    aaa = list(os.walk(path))
    wholefilelist = []
    for tup in aaa:
        wholefilelist = wholefilelist + tup[-1]

    wholefilelist = list(set(wholefilelist))
    lifilelist = [fil for fil in wholefilelist if re.search(suffix_regex +'.*\gins', fil)]

    statname_lis = sorted(list(set([li[stat_strt:stat_end] for li in lifilelist])))

    list(set(statname_lis))

    datadico = dict()

    for stat in statname_lis:
        datadico[stat] = []

    for li in lifilelist:
        try:
            tup = (li , jjulCNES2dt(li[date_strt:date_end] ))
            datadico[li[stat_strt:stat_end]].append(tup)
        except:
            print('error with : ', li)
    plotstat_lis = []
    for i,stat in enumerate(sorted(datadico.keys())):
        print(i , stat)
        T = [e[-1] for e in datadico[stat]]
        plotstat_lis.append(stat)
        ax.plot(T,i*np.ones(len(T)), '.')

    i_list = np.arange(0,len(plotstat_lis))

    print(i_list)
    plt.yticks(i_list,plotstat_lis)

    return fig





def randomwalk_normal(N=100, d=2 , moy = 0 , sigma = 1):
    """
    d = dimension
    """
    return np.cumsum(moy + np.random.randn(N,d) * sigma)

def randomwalk_uniform(N=100, d=2 , bound = 0.5):
    """
    d = dimension
    bound = contraint of the random walk
    """
    return np.cumsum(np.random.uniform(-bound,bound,(N,d)))


def circle_draw(xc,yc,R,N):
    theta = np.linspace(0,2 * np.pi,N)
    X = np.cos(theta) * R + xc
    Y = np.sin(theta) * R + yc
    return X,Y


def random_walk_in_a_circle(x0 , y0 , xc , yc ,
                            R , N , step_size ,  param = 1 ,
                            polar = True , uniform_or_normal = 'n',
                            rand_seed = -1):

    """
    random : normal ou uniform
    coords : polaire ou cartésien

    param est un paramètre très versatile pour controler le random :
    plage pour le uniform
    sigma pour le normal

    on recommande plutot le polaire normal : on a un pas constant & une derive réaliste sur le cap

    Returns :
        X,Y, Xcircle , Ycircle

    Exemple :

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

    Xcircle,Ycircle = circle_draw(xc,yc,R,500)

    for i in range(N-1):
        D = R+1
        iwhil = 0
        while D > R:
            iwhil += 1
            if iwhil > 500:
                print('WARN : infinite loop in random_walk_in_a_circle ...' , iwhil)
            if polar:
                if uniform_or_normal == 'u':
                    dalpha = RAND.uniform(-param,param) * 2 * np.pi
                else:
                    dalpha = RAND.normal(0,param)       * 2 * np.pi
                drho = step_size
                dx = drho * np.cos(dalpha)
                dy = drho * np.sin(dalpha)
                #print dx,dy
            else:
                if uniform_or_normal == 'u':
                    dx = np.random.uniform(-param,param)
                    dy = np.random.uniform(-param,param)
                else:
                    dx = np.random.normal(0,param)
                    dy = np.random.normal(0,param)
                print(dx , dy)
            xtemp = X[-1] + dx
            ytemp = Y[-1] + dy
            D = np.sqrt((xtemp - xc)**2 + (ytemp - yc)**2)

        X.append(xtemp)
        Y.append(ytemp)

    X = np.array(X)
    Y = np.array(Y)

    return X,Y, Xcircle , Ycircle



def randn_bool(N,true_ratio = 0.5,RandGene = None):
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


def points_circle_border(Npts,r,r_sigma,az_type_normal=True,
                         main_dir=3.14159,dir_range=3.14159,seed=None):
    if not seed:
        seed = np.random.randint(10000)

    S = np.random.RandomState(seed)

    if not az_type_normal:
        Az = S.rand(Npts)  * 2 * np.pi
    else:
        Az = S.randn(Npts) * dir_range + main_dir


    R = np.array(Npts * [r]) - np.abs(S.randn(Npts) * r_sigma)

    X,Y = polar2cartesian(R,Az,'rad')

    return X , Y


def estimated_autocorrelation(x):
    """
    http://stackoverflow.com/q/14297012/190597
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    """
    n = len(x)
    variance = x.var()
    x = x-x.mean()
    r = np.correlate(x, x, mode = 'full')[-n:]
    assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
    result = r/(variance*(np.arange(n, 0, -1)))
    return result

def savage_buford_formula(Vs,X,d):
    """
    X : distance à la faille , un iterable pour toutes le profil,
    un nombre pour la longeur max

    d : profondeur de la faille
    retourne X , et Vdeform(X)

    X et d doivent être dans la même unité, Vs pas forcément
    """

    if not genefun.is_iterable(X):
        X = np.arange(-X,X,1)
    return X , ( Vs / np.pi ) * np.arctan2(X,d)




def R2_calc(y_obs,y_fit,with_r2_bis=False):
    #https://en.wikipedia.org/wiki/Coefficient_of_determination
    ybar = np.mean(y_obs)
    SStot = np.sum((y_obs - ybar)**2)
    SSreg = np.sum((y_fit - ybar)**2)
    SSres = np.sum((y_obs - y_fit)**2)

    r2 = 1. - ( SSres / SStot)
    r2bis = ( SSreg / SStot)

    if not with_r2_bis:
        return r2
    else:
        return r2 , r2bis

def R2_from_a_line_regress(Xobs,Yobs,a,b):
    #https://en.wikipedia.org/wiki/Coefficient_of_determination
    Xfit , Yfit = geok.linear_reg_getvalue(Xobs,a,b)
    r2 = R2_calc(Yobs,Yfit)
    return r2




def greek_alphabet(num=None,maj=False):
    if not num:
        a = ['\u03B1',
        '\u03B2',
        '\u03B3',
        '\u03B4',
        '\u03B5',
        '\u03B6',
        '\u03B7',
        '\u03B8',
        '\u03B9',
        '\u03BA',
        '\u03BB',
        '\u03BC',
        '\u03BD',
        '\u03BE',
        '\u03BF',
        '\u03C0',
        '\u03C1',
        '\u03C3',
        '\u03C4',
        '\u03C5',
        '\u03C6',
        '\u03C7',
        '\u03C8',
        '\u03C9']

        A=['\u0391',
        '\u0392',
        '\u0393',
        '\u0394',
        '\u0395',
        '\u0396',
        '\u0397',
        '\u0398',
        '\u0399',
        '\u039A',
        '\u039B',
        '\u039C',
        '\u039D',
        '\u039E',
        '\u039F',
        '\u03A0',
        '\u03A1',
        '\u03A3',
        '\u03A4',
        '\u03A5',
        '\u03A6',
        '\u03A7',
        '\u03A8',
        '\u03A9']

        if maj:
            return A
        else:
            return a
    else:
        return greek_alphabet()[num]


def gebco_bathy_grid_extractor(dataset,latmin,latmax,lonmin,lonmax):
    """
    for safety reasons, lat and lon input MUST BE in the dataset,
    replaced by the closest elsewhere

    return latnew , lonnew , Znew
    """

    import genefun as gf

    lon = dataset['lon'][:]
    lat = dataset['lat'][:]
    Z   = dataset['elevation'][:]

    latmin_dec = latmin
    latmax_dec = latmax
    lonmin_dec = lonmin
    lonmax_dec = lonmax

    #latmin_dec = geok.dms2dec_num(*latmin)
    #latmax_dec = geok.dms2dec_num(*latmax)
    #lonmin_dec = geok.dms2dec_num(*lonmin)
    #lonmax_dec = geok.dms2dec_num(*lonmax)


    if not np.any(latmin_dec == lat):
        print("ERR : ... replacing to the nearest")
        latmin_dec = gf.find_nearest(lat,latmin_dec)[0]

    if not np.any(latmax_dec == lat):
        print("ERR : ... replacing to the nearest")
        latmax_dec = gf.find_nearest(lat,latmax_dec)[0]

    if not np.any(lonmin_dec == lon):
        print("ERR : ... replacing to the nearest")
        lonmin_dec = gf.find_nearest(lon,lonmin_dec)[0]

    if not np.any(lonmax_dec == lon):
        print("ERR : ... replacing to the nearest")
        lonmax_dec = gf.find_nearest(lon,lonmax_dec)[0]

    boollat = np.logical_and(latmin_dec <= lat, lat <= latmax_dec)
    boollon = np.logical_and(lonmin_dec <= lon, lon <= lonmax_dec)


    gridlat = np.tile( boollat[:,np.newaxis]  , (1,len(lon)))
    gridlon = np.tile( boollon , (len(lat),1))

    Znew   = Z[ gridlon * gridlat ]
    latnew = lat[boollat]
    lonnew = lon[boollon]
    Znew   = Znew.reshape(len(latnew),len(lonnew))

    return latnew , lonnew , Znew





