#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
from scipy.interpolate import interp1d
from scipy.stats import median_abs_deviation

#################################
# Variable d'environnement bash
#################################
WORK_DIR = os.getenv('WORK_DIR')


##########################
# Custom cursor for plot
##########################
class SnaptoCursor(object):
    def __init__(self, ax, x, y):
        self.ax = ax
        self.ly = ax.axvline(color='k', alpha=0.2)  # the vert line
        self.marker, = ax.plot([0],[0], marker="o", color="crimson", zorder=3) 
        self.x = x
        self.y = y
        self.txt = ax.text(0.7, 0.9, '')

    def mouse_move(self, event):
        if not event.inaxes: return
        x, y = event.xdata, event.ydata
        indx = np.searchsorted(self.x, [x])[0]
        x = self.x[indx]
        y = self.y[indx]
        self.ly.set_xdata(x)
        self.marker.set_data([x],[y])
        self.txt.set_text(' t=%.6f\n y=%.2f' % (x, y))
        self.txt.set_position((x,y))
        self.ax.figure.canvas.draw_idle()


######################
# Heavyside function
######################
def hvsd (t,t0=0):
    if t >= t0:
        res=1
    if t < t0:
        res=0
    return res
heavy = np.vectorize(hvsd)


################
# Fit function
################
def fit_func(t,t0,fit_name,SAUTS,*param):
    fit = np.zeros(len(t))
    num = 0
    if 'T' in fit_name: # ----------------------------- SLOPE
        fit += param[num] + param[num+1]*(t-t0)
        num += 2
    if 'C' in fit_name: # ----------------------------- ACCELERATION
        fit += param[num]*(t-t0)*(t-t0)
        num += 2
    if  ('a' in fit_name) | ('A' in fit_name): # ------ ANNUAL SIN/COS
        fit += param[num]*np.sin(2*np.pi*t) + param[num+1]*np.cos(2*np.pi*t)
        num += 2
    if 'A' in fit_name: # ----------------------------- SEMI-ANNUAL SIN/COS
        fit += param[num]*np.sin(4*np.pi*t) + param[num+1]*np.cos(4*np.pi*t)
        num += 2
    if ('d' in fit_name) | ('D' in fit_name): # ------- DRACONITIC FUNDAMENTAL SIN/COS
        drac = 1.04 #draconitic frequency
        fit += param[num]*np.sin(drac*2*np.pi*t) + param[num+1]*np.cos(drac*2*np.pi*t)
        num += 2
    if 'D' in fit_name: # ----------------------------- DRACONITIC HARMONICS SIN/COS
        drac = 1.04 #draconitic frequency
        for i in np.arange(7):
            fit += param[num]*np.sin((i+2)*drac*2*np.pi*t) + param[num+1]*np.cos((i+2)*drac*2*np.pi*t)
            num += 2
    if ('O' in fit_name) & (len(SAUTS) != 0): # ------- OFFSETS
        for i,s in enumerate(SAUTS):
            fit += param[num]*heavy(t,t0=s)
            num += 1
    return fit


######################
# Slice fit function
######################
def slice_fit(fit_name,fit_part,SAUTS):
    slice_ = []
    num=0
    if set(fit_name).intersection(fit_part) != set(fit_part):
        exit()

    if 'T' in fit_name: # ----------------------------- SLOPE
        if 'T' in fit_part:
            slice_.append(num+0)
            slice_.append(num+1)
        num += 2
    if 'C' in fit_name: # ----------------------------- ACCELERATION
        if 'C' in fit_part:
            slice_.append(num+0)
            slice_.append(num+1)
        num += 2
    if  ('a' in fit_name) | ('A' in fit_name): # ------ ANNUAL SIN/COS
        if ('a' in fit_part) | ('A' in fit_part) :
            slice_.append(num+0)
            slice_.append(num+1)
        num += 2
    if 'A' in fit_name: # ----------------------------- SEMI-ANNUAL SIN/COS
        if 'A' in fit_part :
            slice_.append(num+0)
            slice_.append(num+1)
        num += 2
    if ('d' in fit_name) | ('D' in fit_name): # ------- DRACONITIC FUNDAMENTAL SIN/COS
        if ('d' in fit_part) | ('D' in fit_part) :
            slice_.append(num+0)
            slice_.append(num+1)
        num += 2
    if 'D' in fit_name: # ----------------------------- DRACONITIC HARMONICS SIN/COS
        if 'D' in fit_part :
            slice_.append(num+0)
            slice_.append(num+1)
            slice_.append(num+2)
            slice_.append(num+3)
            slice_.append(num+4)
            slice_.append(num+5)
            slice_.append(num+6)
            slice_.append(num+7)
            slice_.append(num+8)
            slice_.append(num+9)
            slice_.append(num+10)
            slice_.append(num+11)
            slice_.append(num+12)
            slice_.append(num+13)
        num += 7*2
    if ('O' in fit_name) & (len(SAUTS) != 0): # ------- OFFSETS
        if 'O' in fit_part :
            for i,s in enumerate(SAUTS):
                slice_.append(num+0)
                num += 1
    return slice_



############################
# Residuals of time series
############################
def residuals_serie(time_series_file,fit):
    path="/".join(time_series_file.split("/")[:-2])
    name=time_series_file.split("/")[-1][:9]

    param_east=path+"/MODEL/"+name+"_PARAM_"+fit+".east"
    param_nort=path+"/MODEL/"+name+"_PARAM_"+fit+".nort"
    param_vert=path+"/MODEL/"+name+"_PARAM_"+fit+".vert"

    SOL = np.loadtxt(time_series_file)
    CATS_east = np.loadtxt(param_east,usecols=[3])
    CATS_nort = np.loadtxt(param_nort,usecols=[3])
    CATS_vert = np.loadtxt(param_vert,usecols=[3])
    CATS_name = np.loadtxt(param_vert,usecols=[0],dtype='str')
    t0 = CATS_east[0]

    CATS_east=np.delete(CATS_east,np.where((CATS_name=="(T_REF") | (CATS_name=="(PL") | (CATS_name=="(WH") | (CATS_name=="(INDEX")),0)
    CATS_nort=np.delete(CATS_nort,np.where((CATS_name=="(T_REF") | (CATS_name=="(PL") | (CATS_name=="(WH") | (CATS_name=="(INDEX")),0)
    CATS_vert=np.delete(CATS_vert,np.where((CATS_name=="(T_REF") | (CATS_name=="(PL") | (CATS_name=="(WH") | (CATS_name=="(INDEX")),0)

    with open(param_east) as file:
        lines = file.readlines()
        SAUTS = [float(line.rstrip().split()[-1]) for line in lines if "(OFFSET : )" in line]

    if ".cats" in time_series_file:
        time=0
        east_pos=2
        nort_pos=1
        vert_pos=3
    elif ".spotgins" in time_series_file:
        time=1
        east_pos=3
        nort_pos=4
        vert_pos=5
    else:
        print("ERROR : time series file format is not recognized.")
        print("Only use .cats & .spotgins files !")
        exit()

    ###########################################
    # Residuals between series and CATS model
    ###########################################
    res_east = SOL[:,east_pos]-1.e-3*fit_func(SOL[:,time], t0, fit, SAUTS, *CATS_east)
    res_nort = SOL[:,nort_pos]-1.e-3*fit_func(SOL[:,time], t0, fit, SAUTS, *CATS_nort)
    res_vert = SOL[:,vert_pos]-1.e-3*fit_func(SOL[:,time], t0, fit, SAUTS, *CATS_vert)

    return np.column_stack((SOL[:,time],res_east,res_nort,res_vert))



############################
# Residuals of some points
############################
def residuals_points(ac,name,fit,points):
    path=WORK_DIR+"/DATA/"+ac+"/"+name

    param_east=path+"/MODEL/"+name+"_PARAM_"+fit+".east"
    param_nort=path+"/MODEL/"+name+"_PARAM_"+fit+".nort"
    param_vert=path+"/MODEL/"+name+"_PARAM_"+fit+".vert"

    CATS_east = np.loadtxt(param_east,usecols=[3])
    CATS_nort = np.loadtxt(param_nort,usecols=[3])
    CATS_vert = np.loadtxt(param_vert,usecols=[3])
    CATS_name = np.loadtxt(param_vert,usecols=[0],dtype='str')
    t0 = CATS_east[0]

    CATS_east=np.delete(CATS_east,np.where((CATS_name=="(T_REF") | (CATS_name=="(PL") | (CATS_name=="(WH") | (CATS_name=="(INDEX")),0)
    CATS_nort=np.delete(CATS_nort,np.where((CATS_name=="(T_REF") | (CATS_name=="(PL") | (CATS_name=="(WH") | (CATS_name=="(INDEX")),0)
    CATS_vert=np.delete(CATS_vert,np.where((CATS_name=="(T_REF") | (CATS_name=="(PL") | (CATS_name=="(WH") | (CATS_name=="(INDEX")),0)

    with open(param_east) as file:
        lines = file.readlines()
        SAUTS = [float(line.rstrip().split()[-1]) for line in lines if "(OFFSET : )" in line]

    ###########################################
    # Residuals between series and CATS model
    ###########################################
    res_east = points[:,1]-1.e-3*fit_func(points[:,0], t0, fit, SAUTS, *CATS_east)
    res_nort = points[:,2]-1.e-3*fit_func(points[:,0], t0, fit, SAUTS, *CATS_nort)
    res_vert = points[:,3]-1.e-3*fit_func(points[:,0], t0, fit, SAUTS, *CATS_vert)

    return np.column_stack((points[:,0],res_east,res_nort,res_vert))



#######################################################
# Remove outliers from raw data (median week average)
#######################################################
def remove_outliers_raw(time_series_file,tol=3.):
    TS_raw = np.loadtxt(time_series_file)

    if ".cats" in time_series_file:
        time=0
        east_pos=2
        nort_pos=1
        vert_pos=3
        east_err=5
        nort_err=4
        vert_err=6
    elif ".spotgins" in time_series_file:
        time=1
        east_pos=3
        nort_pos=4
        vert_pos=5
        east_err=6
        nort_err=7
        vert_err=8
    else:
        print("ERROR : time series file format is not recognized.")
        print("Only use .cats & .spotgins files !")
        exit()

    ###################################################
    # Delete obvious outliers with huge formal errors
    ###################################################
    cond=(TS_raw[:,east_err] <= 0.01) & (TS_raw[:,nort_err] <= 0.01) & (TS_raw[:,vert_err] <= 0.01)
    TS = TS_raw[cond,:]

    #############################################
    # Make week-median series to detect ouliers
    #############################################
    TS_week=np.zeros((1,4))
    step=0.02
    weeks=np.arange(TS[0,time]+step,TS[-1,time],step)
    weeks[-1]=TS[-1,time]+0.0001
    BEG=TS[0,time]
    for END in weeks:
        cond=(TS[:,time] >= BEG) & (TS[:,time] <= END)
        if np.any(cond):
            med = np.median(TS[cond,:][:,[time,east_pos,nort_pos,vert_pos]],axis=0)
            TS_week = np.concatenate((TS_week,med[np.newaxis,:]))
        BEG=END
    TS_week=np.delete(TS_week,0,0)

    #############################################
    # Rebound the week series for interpolation
    #############################################
    TS_week[0,0]=TS[0,time]
    TS_week[-1,0]=TS[-1,time]

    ########################################
    # Interpolate the week series linearly
    ########################################
    east_interp = interp1d(TS_week[:,0],TS_week[:,1],kind='linear')
    nort_interp = interp1d(TS_week[:,0],TS_week[:,2],kind='linear')
    vert_interp = interp1d(TS_week[:,0],TS_week[:,3],kind='linear')

    #########################################################################
    # Calculate the residuals between day and interpolated_week time series
    #########################################################################
    res_east = TS[:,east_pos] - east_interp(TS[:,time])
    res_nort = TS[:,nort_pos] - nort_interp(TS[:,time])
    res_vert = TS[:,vert_pos] - vert_interp(TS[:,time])

    ###################################################################
    # Condition for points to not be considered as outliers [< tol*MAD]
    ###################################################################
    cond_east = (np.abs(res_east-np.median(res_east)) <= tol*1.4826*median_abs_deviation(res_east))
    cond_nort = (np.abs(res_nort-np.median(res_nort)) <= tol*1.4826*median_abs_deviation(res_nort))
    cond_vert = (np.abs(res_vert-np.median(res_vert)) <= tol*1.4826*median_abs_deviation(res_vert))
    cond = cond_east & cond_nort & cond_vert

    return TS[cond,:]


#######################################
# Remove outliers from fit CATS model
#######################################
def remove_outliers_fit(time_series_file,fit,tol=3.):
    res = residuals_serie(time_series_file,fit)
    res_east = res[:,1]
    res_nort = res[:,2]
    res_vert = res[:,3]
    cond_east = (np.abs(res_east-np.median(res_east)) <= tol*1.4826*median_abs_deviation(res_east))
    cond_nort = (np.abs(res_nort-np.median(res_nort)) <= tol*1.4826*median_abs_deviation(res_nort))
    cond_vert = (np.abs(res_vert-np.median(res_vert)) <= tol*1.4826*median_abs_deviation(res_vert))
    cond = cond_east & cond_nort & cond_vert

    SOL = np.loadtxt(time_series_file)
    return SOL[cond,:]


####################################################
# Compare newest solutions stats with previous fit
####################################################
def compare_std2sol(time_series_file,fit,points):
    ac=time_series_file.replace(WORK_DIR+"/DATA/","").split("/")[0]
    name=time_series_file.replace(WORK_DIR+"/DATA/","").split("/")[1]
    
    res = residuals_serie(time_series_file,fit)
    res_east = res[:,1]
    res_nort = res[:,2]
    res_vert = res[:,3]
    std_east = 1.4826*median_abs_deviation(res_east)
    std_nort = 1.4826*median_abs_deviation(res_nort)
    std_vert = 1.4826*median_abs_deviation(res_vert)

    res_points = residuals_points(ac,name,fit,points)

    for i in range(len(res_points[:,0])):
        print(abs(int(res_points[i,1]/std_east))+1, abs(int(res_points[i,2]/std_nort))+1, abs(int(res_points[i,3]/std_vert))+1)

