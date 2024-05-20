#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.reffram contains functions for operations
related to GNSS-products


it can be imported directly with:
from geodezyx import reffram

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
import itertools
#### Import the logger
import logging
import os
import re

import matplotlib
import matplotlib.pyplot as plt
import natsort
import numpy as np
import pandas as pd
import pyorbital.astronomy

#### geodeZYX modules
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import stats
from geodezyx import utils
from geodezyx.reffram import kepler_gzyx

### disabled and imported directly in the needed fct
## import geodezyx.reffram.sofa18 as sofa
log = logging.getLogger(__name__)

##########  END IMPORT  ##########

def compar_orbit(Data_inp_1,Data_inp_2,step_data = 900,
                 sats_used_list = ['G'],
                 name1='',name2='',use_name_1_2_for_table_name = False,
                 RTNoutput = True,convert_ECEF_ECI=True,
                 clean_null_values = True,
                 conv_coef=10**3,return_satNull = False):
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
        convert sp3 ECEF => ECI (Terrestrial => Celestrial)
        must be True in operational to avoid artifacts.

    name1 & name2 : str (optionals)
        optional custom names for the 2 orbits

    use_name_1_2_for_table_name : bool
        False : use name 1 and 2 for table name, use datafile instead

    clean_null_values : bool or str
        if True or "all" remove sat position in all X,Y,Z values
        are null (0.000000)
        if "any", remove sat position if X or Y or Z is null
        if False, keep everything
        
    conv_coef : int
        conversion coefficient, km to m 10**3, km to mm 10**6

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

    References
    ----------
    "Coordinate Systems", ASEN 3200 1/24/06 George H. Born

    """

    # selection of both used Constellations AND satellites
    sys_used_list = []
    prn_used_list  = []
    for e in sats_used_list:
        if len(e) == 1: ## it is a constellation
            sys_used_list.append(e)
        elif len(e) == 3: ## it is a satellite
            prn_used_list.append(e)
            if not e[0] in sys_used_list:
                sys_used_list.append(e[0])

    # Read the files or DataFrames
    # metadata attributes are not copied
    # Thus, manual copy ...
    # (Dirty way, should be impoved without so many lines ...)
    if type(Data_inp_1) is str:
        D1orig = files_rw.read_sp3(Data_inp_1,epoch_as_pd_index=True)
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
        D2orig = files_rw.read_sp3(Data_inp_2,epoch_as_pd_index=True)
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
            sat_nul = utils.join_improved(" " ,*list(set(D1orig[D1_null_bool]["prn"])))
            log.warning("Null values contained in SP3 files : ")
            log.warning("f1: %s %s", np.sum(D1_null_bool),
                        utils.join_improved(" ", *list(set(D1orig[D1_null_bool]["prn"]))))
            log.warning("f2: %s %s", np.sum(D2_null_bool),
                        utils.join_improved(" ", *list(set(D2orig[D2_null_bool]["prn"]))))
        else:
            sat_nul = []

    else:
        D1 = D1orig.copy()
        D2 = D2orig.copy()

    D1sys_grp = D1.groupby("sys")
    D2sys_grp = D2.groupby("sys")

    for sysuse in sys_used_list:
        D1sys = D1sys_grp.get_group(sysuse) # D1sys = D1[D1['sys'] == sysuse]
        D2sys = D2sys_grp.get_group(sysuse) # D2sys = D2[D2['sys'] == sysuse]

        # checking if the data correspond to the step
        bool_step1 = np.mod((D1sys.index - np.min(D1.index)).seconds,step_data) == 0
        bool_step2 = np.mod((D2sys.index - np.min(D2.index)).seconds,step_data) == 0

        D1win = D1sys[bool_step1]
        D2win = D2sys[bool_step2]
        
        # find common sats and common epochs
        prni_set = sorted(list(set(D1win['prni']).intersection(set(D2win['prni']))))
        epoc_set = sorted(list(set(D1win.index).intersection(set(D2win.index))))

        # if special selection of sats, then apply it
        # (it is late and this selection is incredibely complicated ...)
        if np.any([True if sysuse in e else False for e in prn_used_list]):
            # first find the selected sats for the good constellation
            prni_used_select_list = [int(e[1:]) for e in prn_used_list if sysuse in e]
            # and apply it
            prni_set = sorted(list(set(prni_set).intersection(set(prni_used_select_list))))

        D1win_prni_grp = D1win.groupby("prni")
        D2win_prni_grp = D2win.groupby("prni")

        for prni in prni_set:
            # First research : find corresponding epoch for the SV
            # this one is sufficent if there is no gaps (e.g. with 0.00000) i.e.
            # same nb of obs in the 2 files
            # NB : .reindex() is smart, it fills the DataFrame
            # with NaN
            try:
                D1prni_orig = D1win_prni_grp.get_group(prni).reindex(epoc_set)
                # D1win[D1win['prni'] == prni].reindex(epoc_set)
                D2prni_orig = D2win_prni_grp.get_group(prni).reindex(epoc_set)
                # D2win[D2win['prni'] == prni].reindex(epoc_set)
            except Exception as exce:
                log.info("ERR : Unable to re-index with an unique epoch")
                log.info("      are you sure there is no multiple-defined epochs for the same sat ?")
                log.info("      it happens e.g. when multiple ACs are in the same DataFrame ")
                log.info("TIP : Filter the input Dataframe before calling this fct with")
                log.info("      DF = DF[DF['AC'] == 'gbm']")
                
                Dtmp1 = D1orig[D1orig['prni'] == prni]
                Dtmp2 = D2orig[D2orig['prni'] == prni]
                
                dupli1 = np.sum(Dtmp1.duplicated(["epoch","prn"]))
                dupli2 = np.sum(Dtmp2.duplicated(["epoch","prn"]))
                
                log.info("FWIW: duplicated epoch/sat in DF1 & DF2: %s %s",dupli1,dupli2)

                raise exce

            # Second research, it is a security in case of gap
            # This step is useless, because .reindex() will fill the DataFrame
            # with NaN
            if len(D1prni_orig) != len(D2prni_orig):
                log.info("different epochs nbr for SV %s %s %s",prni, len(D1prni_orig), len(D2prni_orig))
                epoc_prni_set = sorted(list(set(D1prni_orig.index).intersection(set(D2prni_orig.index))))
                D1prni = D1prni_orig.loc[epoc_prni_set]
                D2prni = D2prni_orig.loc[epoc_prni_set]
            else:
                D1prni = D1prni_orig
                D2prni = D2prni_orig

            P1 = D1prni[['x','y','z']]
            P2 = D2prni[['x','y','z']]

            # Start ECEF => ECI
            if convert_ECEF_ECI:
                # Backup because the columns xyz will be reaffected
                #D1sv_bkp = D1prni.copy()
                #D2sv_bkp = D2prni.copy()
    
                P1b = conv.ECEF2ECI(np.array(P1),conv.dt_gpstime2dt_utc(P1.index.to_pydatetime(),out_array=True))
                P2b = conv.ECEF2ECI(np.array(P2),conv.dt_gpstime2dt_utc(P2.index.to_pydatetime(),out_array=True))

                D1prni[['x','y','z']] = P1b
                D2prni[['x','y','z']] = P2b

                P1 = D1prni[['x','y','z']]
                P2 = D2prni[['x','y','z']]
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

                Vx = utils.diff_pandas(D1prni,'x',use_np_diff=True)
                Vy = utils.diff_pandas(D1prni,'y',use_np_diff=True)
                Vz = utils.diff_pandas(D1prni,'z',use_np_diff=True)
                
                V = pd.concat((Vx , Vy , Vz),axis=1)
                V.columns = ['vx','vy','vz']

                R = P1.divide(rnorm,axis=0)
                R.columns = ['xnorm','ynorm','znorm']

                H = pd.DataFrame(np.cross(R,V),columns=['hx','hy','hz'])
                hnorm = np.linalg.norm(H,axis=1)

                C = H.divide(hnorm,axis=0)
                C.columns = ['hxnorm','hynorm','hznorm']

                I = pd.DataFrame(np.cross(C,R),columns=['ix','iy','iz'])

                R_ar = np.array(R)
                I_ar = np.array(I)
                C_ar = np.array(C)

                #R_ar[1]
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
                                        index = P1.index,
                                        columns=['dr','dt','dn'])

            Diff_sat = Diff_sat * conv_coef # metrer conversion

            Diff_sat['sys'] = [sysuse] * len(Diff_sat.index)
            Diff_sat['prni'] = [prni] * len(Diff_sat.index)
            Diff_sat['prn'] = [sysuse + str(prni).zfill(2)] * len(Diff_sat.index)

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
                                  ', doy', str(conv.dt2doy(Date))))
    
    if return_satNull:
        return Diff_sat_all, sat_nul
    else:
        return Diff_sat_all


def compar_orbit_plot(Diff_sat_all_df_in,
                      save_plot=False,
                      save_plot_dir="",
                      save_plot_name="auto",
                      save_plot_name_suffix=None,
                      save_plot_ext=(".pdf",".png",".svg"),
                      yaxis_limit=None,
                      yaxis_label_unit="m"):
    """
    General description

    Parameters
    ----------
    Diff_sat_all_df_in : DataFrame
        a DataFrame produced by compar_orbit
        
    yaxis_limit : 3-tuple iterable or 2-element tuple
        force the y axis limits. must look like 
        [(ymin_r,ymax_r),(ymin_t,ymax_t),(ymin_n,ymax_n)]
        to control all the axis independely
        OR
        (ymin,ymax)
        to set all th axis at the same limits 

    Returns
    -------
    the Figure and the 3 Axes if no save is asked
    export path (str) if save is asked
    but plot a plot anyway
    """

    fig,[axr,axt,axn] = plt.subplots(3,1,sharex='all')

    satdispo = natsort.natsorted(list(set(Diff_sat_all_df_in['prn'])))

    SymbStk = []

    cm = plt.get_cmap('viridis')
    NUM_COLORS = len(satdispo)
    Colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]

    # Pandas donesn't manage well iterable as attribute
    # So, it is separated
    try:
        col_name0 = Diff_sat_all_df_in.frame_col_name1
        col_name1 = Diff_sat_all_df_in.frame_col_name2
        col_name2 = Diff_sat_all_df_in.frame_col_name3
    except:
        col_name0 = Diff_sat_all_df_in.columns[0]
        col_name1 = Diff_sat_all_df_in.columns[1]
        col_name2 = Diff_sat_all_df_in.columns[2]

    for satuse,color in zip(satdispo,Colors):
        Diffuse = Diff_sat_all_df_in[Diff_sat_all_df_in['prn'] == satuse]

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


    ylabuni = " (" + yaxis_label_unit + ")"
    
    if Diff_sat_all_df_in.frame_type == 'RTN':
        axr.set_ylabel('Radial diff.'     + ylabuni)
        axt.set_ylabel('Transverse diff.' + ylabuni)
        axn.set_ylabel('Normal diff.'     + ylabuni)

    else:
        axr.set_ylabel(Diff_sat_all_df_in.frame_type + ' X diff.' + ylabuni)
        axt.set_ylabel(Diff_sat_all_df_in.frame_type + ' Y diff.' + ylabuni)
        axn.set_ylabel(Diff_sat_all_df_in.frame_type + ' Z diff.' + ylabuni)


    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    axr.yaxis.set_major_formatter(y_formatter)
    axt.yaxis.set_major_formatter(y_formatter)
    axn.yaxis.set_major_formatter(y_formatter)
    
    if yaxis_limit and len(yaxis_limit) == 3: ### indep. axis limit
        axr.set_ylim(yaxis_limit[0])
        axt.set_ylim(yaxis_limit[1])
        axn.set_ylim(yaxis_limit[2])
    elif yaxis_limit and len(yaxis_limit) == 2:
        axr.set_ylim(yaxis_limit)
        axt.set_ylim(yaxis_limit)
        axn.set_ylim(yaxis_limit)
    else:
        pass
        
        
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
            
        if save_plot_name_suffix:
            save_plot_name = save_plot_name + '_' + save_plot_name_suffix

        for ext in save_plot_ext:
            save_plot_path = os.path.join(save_plot_dir,save_plot_name)
            plt.savefig(save_plot_path + ext)
            return_val = save_plot_path
            
    else:
        return_val = fig,(axr,axt,axn)

    return return_val

def compar_orbit_table(Diff_sat_all_df_in,RMS_style = 'natural',
                       light_tab  = False):
    """
    Generate a table with statistical indicators for an orbit comparison
    (RMS mean, standard dev, ...)
    
    Parameters
    ----------
    Diff_sat_all_df_in : Pandas DataFrame
        a DataFrame produced by compar_orbit

    RMS_style : str
        'natural': use the natural definition of the RMS
        'GRGS': RMS calc based on the GRGS definition of the RMS (OV help)
        is actually the standard deviation
        'kouba': RMS as defined in Kouba et al. 1994, p75
        using the degree of freedom (3*Nobs - 7)

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

    sat_list = utils.uniq_and_sort(Diff_sat_all_df_in['prn'])

    # Pandas donesn't manage well iterable as attribute
    # So, it is separated
    try:
        col_name0 = Diff_sat_all_df_in.frame_col_name1
        col_name1 = Diff_sat_all_df_in.frame_col_name2
        col_name2 = Diff_sat_all_df_in.frame_col_name3
    except:
        col_name0 = Diff_sat_all_df_in.columns[0]
        col_name1 = Diff_sat_all_df_in.columns[1]
        col_name2 = Diff_sat_all_df_in.columns[2]

    rms_stk = []

    for sat in sat_list:
        Diffwork = utils.df_sel_val_in_col(Diff_sat_all_df_in,'prn',sat)

        if RMS_style == "natural":
            rms_A = stats.rms_mean(Diffwork[col_name0])
            rms_B = stats.rms_mean(Diffwork[col_name1])
            rms_C = stats.rms_mean(Diffwork[col_name2])
        elif RMS_style == "GRGS":
            rms_A = stats.rms_mean(Diffwork[col_name0] - Diffwork[col_name0].mean())
            rms_B = stats.rms_mean(Diffwork[col_name1] - Diffwork[col_name1].mean())
            rms_C = stats.rms_mean(Diffwork[col_name2] - Diffwork[col_name2].mean())
        elif RMS_style == "kouba":
            rms_A = stats.rms_mean_kouba(Diffwork[col_name0])
            rms_B = stats.rms_mean_kouba(Diffwork[col_name1])
            rms_C = stats.rms_mean_kouba(Diffwork[col_name2])
            
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
            rms_stk.append([rms_A,rms_B,rms_C,RMS3D])
        else:
            rms_stk.append([rms_A,rms_B,rms_C,RMS3D,
                            min_A,max_A,mean_A,
                            min_B,max_B,mean_B,
                            min_C,max_C,mean_C])


    #################################
             # ALL SATS
    if RMS_style == "natural":
        rms_A = stats.rms_mean(Diff_sat_all_df_in[col_name0])
        rms_B = stats.rms_mean(Diff_sat_all_df_in[col_name1])
        rms_C = stats.rms_mean(Diff_sat_all_df_in[col_name2])
        RMS3D = np.sqrt(rms_A**2 + rms_B**2 + rms_C**2)
    elif RMS_style == "GRGS":
        rms_A = stats.rms_mean(Diff_sat_all_df_in[col_name0] - Diff_sat_all_df_in[col_name0].mean())
        rms_B = stats.rms_mean(Diff_sat_all_df_in[col_name1] - Diff_sat_all_df_in[col_name1].mean())
        rms_C = stats.rms_mean(Diff_sat_all_df_in[col_name2] - Diff_sat_all_df_in[col_name2].mean())
        RMS3D = np.sqrt(rms_A**2 + rms_B**2 + rms_C**2)
    elif RMS_style == "kouba":
        rms_A = stats.rms_mean_kouba(Diff_sat_all_df_in[col_name0])
        rms_B = stats.rms_mean_kouba(Diff_sat_all_df_in[col_name1])
        rms_C = stats.rms_mean_kouba(Diff_sat_all_df_in[col_name2])
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
        rms_stk.append([rms_A,rms_B,rms_C,RMS3D])
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


def compar_orbit_frontend(DataDF1,DataDF2,ac1,ac2, sats_used_list = ['G']):
    K = compar_orbit(DataDF1[DataDF1["ac"] == ac1],
                     DataDF2[DataDF2["ac"] == ac2],
                     sats_used_list=sats_used_list)
    compar_orbit_plot(K)
    return K





def compar_clock(DFclk_inp_1,DFclk_inp_2,col_name = "name",bias_Col_name = "bias"):
    """
    Compares 2 GNSS clock bias DataFrames (from .clk), to a
    statistics table (with compar_clock_table)


    Parameters
    ----------
    DFclk_inp_1 & DFclk_inp_2 : DataFrame
        Clock DataFrame provided by files_rw.read_clk()

    Returns
    -------
    DFclk_diff : DataFrame
        Clock bias difference DataFrame
    """
    DF1idx = DFclk_inp_1.set_index([col_name,"epoch"])
    DF1idx.sort_index(inplace=True)
    
    DF2idx = DFclk_inp_2.set_index([col_name,"epoch"])
    DF2idx.sort_index(inplace=True)
    
    I1 = DF1idx.index
    I2 = DF2idx.index
    
    Iinter = I1.intersection(I2)
    Iinter = Iinter.sort_values()
    
    DF_diff_bias = DF1idx.loc[Iinter][bias_Col_name] - DF2idx.loc[Iinter][bias_Col_name]

    DFclk_diff = DF1idx.loc[Iinter].copy()
    DFclk_diff[bias_Col_name] = DF_diff_bias
    if "ac" in  DFclk_diff.columns:
        DFclk_diff.drop("ac",axis=1,inplace=True)
    else:
        DFclk_diff.drop("AC",axis=1,inplace=True)
    DFclk_diff.rename({bias_Col_name:bias_Col_name+"_diff"},inplace=True,axis=1)
    
    
    # Name definitions
    if 'ac' in DFclk_inp_1.columns:
        DFclk_diff['name1'] = DFclk_inp_1.ac.values[0]
    if 'AC' in DFclk_inp_1.columns:
        DFclk_diff['name1'] = DFclk_inp_1.AC.values[0]
    
    if 'ac' in DFclk_inp_1.columns:
        DFclk_diff['name2'] = DFclk_inp_2.ac.values[0]
    if 'AC' in DFclk_inp_1.columns:
        DFclk_diff['name2'] = DFclk_inp_2.AC.values[0]

    
    return DFclk_diff
    
def compar_clock_table(DFclk_diff_in,col_name = "name",bias_Col_name = "bias_diff"):
    """
    Generate a table with statistical indicators for a clock comparison
    (RMS mean, standard dev, ...)

    Parameters
    ----------
    DFclk_diff_in : DataFrame
        Clock bias difference DataFrame (from compar_clock)

    Returns
    -------
    DFcompar_out : DataFrame
        Statistical results of the comparison.

    """
    
    DF_diff_grp = DFclk_diff_in.groupby(col_name)[bias_Col_name]
    
    Smin  = DF_diff_grp.min().rename("min",inplace=True)
    Smax  = DF_diff_grp.max().rename("max",inplace=True)
    Smean = DF_diff_grp.mean().rename("mean",inplace=True)
    Sstd  = DF_diff_grp.std().rename("std",inplace=True)
    Srms  = DF_diff_grp.apply(stats.rms_mean).rename("rms",inplace=True)
    
    DFcompar_out = pd.concat([Smin,Smax,Smean,Sstd,Srms],axis=1)
    DFcompar_out.reset_index()
    
    return DFcompar_out

def compar_clk_plot(Diff_sat_all_df_in,
                    save_plot=False,
                    save_plot_dir="",
                    save_plot_name="auto",
                    save_plot_name_suffix=None,
                    save_plot_ext=(".pdf",".png",".svg"),
                    yaxis_limit=None,
                    yaxis_label_unit="psec",
                    col_name = 'name',
                    bias_Col_name = 'bias'):
    """
    General description

    Parameters
    ----------
    Diff_sat_all_df_in: DataFrame
        a DataFrame produced by compar_clk
        
    yaxis_limit: 3-tuple iterable or 2-element tuple
        force the y axis limits. must look like 
        (ymin,ymax) to set all th axis at the same limits 
    
    col_name: Normally the name of the column with the sat names
    bias_Col_name: The column with the clk values
    Default: 'bias'

    Returns
    -------
    the Figure and the 3 Axes if no save is asked
    export path (str) if save is asked
    but plot a plot anyway
    """

    fig,axr = plt.subplots(1,1,sharex='all')
    Diff_sat_all_df_in = Diff_sat_all_df_in.reset_index()
    satdispo = natsort.natsorted(list(set(Diff_sat_all_df_in[col_name])))
    # satdispo = natsort.natsorted(list(set(Diff_sat_all_df_in['sat'])))
       
    SymbStk = []
       
    cm = plt.get_cmap('viridis')
    NUM_COLORS = len(satdispo)
    Colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
       
    
    
    Date = conv.numpy_dt2dt(Diff_sat_all_df_in.epoch.values[0])
    Diff_sat_all_df_in.name = ' '.join(('Clock comparison  b/w',
                                  Diff_sat_all_df_in.name1.values[0] ,'(ref.) and',
                                  Diff_sat_all_df_in.name2.values[0] ,',',Date.strftime("%Y-%m-%d")))
    # Pandas donesn't manage well iterable as attribute
    # So, it is separated
    try:
        col_name0 = Diff_sat_all_df_in.frame_col_name1
        col_name1 = Diff_sat_all_df_in.frame_col_name2
        col_name2 = Diff_sat_all_df_in.frame_col_name3
    except:
        col_name0 = Diff_sat_all_df_in.columns[0]
        col_name1 = Diff_sat_all_df_in.columns[1]
        col_name2 = Diff_sat_all_df_in.columns[2]
       
    for satuse,color in zip(satdispo,Colors):
        Diffuse = Diff_sat_all_df_in[Diff_sat_all_df_in[col_name] == satuse]
       
        Time = Diffuse.epoch
        R    = Diffuse[bias_Col_name+'_diff']*10**12
       
       
        #fig.fmt_xdata = mdates.DateFormatter('%Y-%m-%d')
       
        Symb = axr.plot(Time,R,label=satuse,c=color)
       
       
        SymbStk.append(Symb[0])
       
        fig.autofmt_xdate()
       
       
    ylabuni = " (" + yaxis_label_unit + ")"
    
       
    axr.set_ylabel('Bias Diff.' + ylabuni)
       
       
       
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    axr.yaxis.set_major_formatter(y_formatter)
       
    
    if yaxis_limit and len(yaxis_limit) == 3: ### indep. axis limit
        axr.set_ylim(yaxis_limit[0])
       
    elif yaxis_limit and len(yaxis_limit) == 2:
        axr.set_ylim(yaxis_limit)
       
    else:
        pass
        
        
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
            save_plot_name = Diff_sat_all_df_in.name1.values[0]+"_"+Diff_sat_all_df_in.name2.values[0]+"_"+Date.strftime("%Y-%m-%d")
            
        if save_plot_name_suffix:
            save_plot_name = save_plot_name + '_' + save_plot_name_suffix
       
        for ext in save_plot_ext:
            save_plot_path = os.path.join(save_plot_dir,save_plot_name)
            plt.savefig(save_plot_path + ext)
            return_val = save_plot_path
            
    else:
        return_val = fig,(axr)

    return return_val






def compar_sinex(snx1 , snx2 , stat_select = None, invert_select=False,
                 out_means_summary=True,out_meta=True,out_dataframe = True,
                 manu_wwwwd=None):

    if type(snx1) is str:
        week1 = utils.split_improved(os.path.basename(snx1),"_",".")[:]
        week2 = utils.split_improved(os.path.basename(snx2),"_",".")[:]
        if week1 != week2:
            log.warning("Dates of 2 input files are differents !!! It might be very bad !!! %s %s",week1,week2)
        else:
            wwwwd = week1
        D1 = files_rw.read_sinex(snx1,True)
        D2 = files_rw.read_sinex(snx2,True)
    else:
        log.warning("WARN : you are giving the SINEX input as a DataFrame, wwwwd has to be given manually using manu_wwwwd")
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
        elif utils.is_iterable(stat_select):
            STATCommon = [sta for sta in STATCommon_init if select_fct(sta in stat_select) ]
        else:
            log.warning("WARN : check type of stat_select")

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
        enu   = conv.XYZ2ENU_2(l1["x"],l1["y"],l1["z"],l2["x"],l2["y"],l2["z"])
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

    #    E,N,U    = conv.XYZ2ENU_2((X,Y,Z,x0,y0,z0))
    #    E,N,U    = conv.XYZ2ENU_2((X,Y,Z,x0,y0,z0))

    if out_dataframe:
        out_meta = True


    if not out_means_summary:
        log.info("INFO : this is not used operationally and it can be improved")
        return Ddiff
    else:
        output = []

        col_names = ("x","y","z","d3D_xyz",
                     "e","n","u","d2D_enu","d3D_enu")

        for xyz in col_names:
            output.append(stats.rms_mean(Ddiff[xyz]))
        for xyz in col_names:
            output.append(np.nanmean(Ddiff[xyz]))
        for xyz in col_names:
            output.append(np.nanstd(Ddiff[xyz]))

        if out_meta:
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


#   ____       _     _ _     _____        _        ______                             
#  / __ \     | |   (_) |   |  __ \      | |      |  ____|                            
# | |  | |_ __| |__  _| |_  | |  | | __ _| |_ __ _| |__ _ __ __ _ _ __ ___   ___  ___ 
# | |  | | '__| '_ \| | __| | |  | |/ _` | __/ _` |  __| '__/ _` | '_ ` _ \ / _ \/ __|
# | |__| | |  | |_) | | |_  | |__| | (_| | || (_| | |  | | | (_| | | | | | |  __/\__ \
#  \____/|_|  |_.__/|_|\__| |_____/ \__,_|\__\__,_|_|  |_|  \__,_|_| |_| |_|\___||___/
#                                                                                     
       
### Orbit DataFrames   


def DFOrb_velocity_calc(DFOrb_in,
                        drop_nan=False):
    """
    Compute the velocity of satellites from a DFOrb dataframe
    (differentiate the position)

    Parameters
    ----------
    DFOrb_in : Pandas DataFrame
        an Orbit DataFrame.
    drop_nan : bool, optional
        Remove the nan values 
        (the first ones, since it is a numerical differentiation).
        The default is False.
        It is recommended to keep the NaN values, thus the output will keep
        same size and index as input

    Returns
    -------
    df_vel : Pandas DataFrame
        an Orbit DataFrame with velocities (vx, vy ,vz columns).

    """
    
    dfgrp = DFOrb_in.groupby('prn')
    
    dfprn_stk = []
    for prn,dfprn in dfgrp:    
        for coord in ['x','y','z']:
            dcoord = dfprn[coord].diff()
            dtime = (np.float64(dfprn['epoch'].diff()) * 10**-9) 
            #timedelta in ns per defalut
            dfprn['v' + coord] = dcoord / dtime 
        dfprn_stk.append(dfprn)
            
    df_vel = pd.concat(dfprn_stk)
    
    df_vel.sort_index(inplace=True)
    
    if drop_nan:
        df_vel.dropna(inplace=True)
    
    return df_vel 


def beta_sun_ra_dec(sun_dec,sun_ra,sat_i,sat_o_lan):
    """
    Compute beta angle based on Sun's right ascension and declination
    Angles are in radians

    Parameters
    ----------
    sun_dec : float
        Sun's declination.
    sun_ra : float
        Sun's right ascension.
    sat_i : float
        Satellite's inclination.
    sat_o_lan : float
        Satellite's Longitude of the ascending node.

    Returns
    -------
    beta : float
        beta angle (in radians).
        
    Source
    ------
    `Satellites: de Kepler au GPS`, Michel Capderou (2012)
    
    simpler formula here:
        https://www.fxsolver.com/browse/formulas/Beta+Angle
        
    
    Note
    ----
    you can use ``pyorbital.astronomy.sun_ra_dec()`` to compute ``sun_dec``
    and ``sun_ra``
    
    If so, Sun's right ascension and declination computation polynoms are 
    based on: `Astronomical algorithms`, Jean Meeus (1st edition, 1991)
    """
    beta = np.arcsin(np.cos(sun_dec)*np.sin(sat_i)*np.sin(sat_o_lan-sun_ra)+
                     np.sin(sun_dec)*np.cos(sat_i))
    return beta

def beta_sun_eclip_long(sun_ecl_long,sat_o_lan,sat_i,earth_i):
    """
    Compute beta angle based on Sun's Ecliptic longitude
    Angles are in radians

    Parameters
    ----------
    sun_ecl_long : float
        Sun's Ecliptic longitude.
    sat_i : float
        Satellite's inclination.
    sat_o_lan : float
        Satellite's Longitude of the ascending node.
    earth_i : float
        Earth's inclination.

    Returns
    -------
    beta : float
        beta angle (in radians).
        
    Source
    ------
    `Calculation of the Eclipse Factor for Elliptical Satellite Orbits` NASA (1962)
    https://ntrs.nasa.gov/citations/19630000622
    `Computation of Eclipse Time for Low-Earth Orbiting Small Satellites` Sumanth R M (2019)
    https://commons.erau.edu/ijaaa/vol6/iss5/15/

    Simpler formula:
        https://en.wikipedia.org/wiki/Beta_angle
    
    Note
    ----
    you can use ``pyorbital.astronomy.sun_ecliptic_longitude()``
    to compute ``sun_ecl_long``
    
    If so, Sun's right ascension and declination computation polynoms are 
    based on: `Astronomical algorithms`, Jean Meeus (1st edition, 1991)

    """
    p1 = np.cos(sun_ecl_long)*np.sin(sat_o_lan)*np.sin(sat_i)
    p2 = np.sin(sun_ecl_long)*np.cos(earth_i)*np.cos(sat_o_lan)*np.sin(sat_i)
    p3 = np.sin(sun_ecl_long)*np.sin(earth_i)*np.cos(sat_i)
    beta = np.arcsin(p1-p2+p3)
    return beta

def beta_angle_calc(DFOrb_in,
                    calc_beta_sun_ra_dec=True,
                    calc_beta_sun_eclip_long=True,
                    beta_rad2deg=True):
    """
    Compute beta angle for GNSS satellite's orbits stored in an orbit 
    DataFrame
    

    Parameters
    ----------
    DFOrb_in : Pandas DataFrame
        an Orbit DataFrame (ECEF frame).
    calc_beta_sun_ra_dec : bool, optional
        compute beta angle with Sun's right ascension and declination.
        The default is True.
    calc_beta_sun_eclip_long : bool, optional
        compute beta angle with Sun's Ecliptic longitude.
        The default is True.
    beta_rad2deg : bool, optional
        convert beta angle to degrees. The default is True.

    Returns
    -------
    df_out : Pandas DataFrame
        DFOrb_in with a new 'beta' column.
    df_wrk : Pandas DataFrame
        intermediate values dataframe for debug.
        here, coordinates are in ECI frame
    """
    
    #### convert in ECI
    df_eci = DFOrb_in.copy()
    df_eci[['x','y','z']] = conv.ECEF2ECI(DFOrb_in[['x','y','z']].values,
                                          DFOrb_in['epoch'].values)
    
    #### compute velocity
    df_wrk = DFOrb_velocity_calc(df_eci,drop_nan=False)
    ##keep drop nan False, thus the DF will keep same size and index as input
    
    #### compute Kepler's parameters
    p = df_wrk[['x','y','z']].values * 1000
    v = df_wrk[['vx','vy','vz']].values * 1000
    
    kep_col = ['a','ecc','i','o_peri','o_lan','m']
    kep_params = kepler_gzyx.ECI_2_kepler_elts(p,v,rad2deg=False)
    df_wrk[kep_col] = np.column_stack(kep_params)
    ######### COMPUTE BETA
    
    #### cosmetic changes
    if calc_beta_sun_ra_dec and calc_beta_sun_eclip_long:
        b1="1"
        b2="2"
    else:
        b1=""
        b2=""    
    
    if beta_rad2deg:
        r2dfct = np.rad2deg
    else:
        r2dfct = lambda x:x 
    
    ##### Beta computed based on sun declination / right ascension
    if calc_beta_sun_ra_dec:
        ##### sun_ra_dec output in RADIANS
        df_wrk[['sun_ra','sun_dec']] = np.column_stack(pyorbital.astronomy.sun_ra_dec(df_wrk['epoch']))
        df_wrk['beta'+b1] = beta_sun_ra_dec(df_wrk['sun_dec'], 
                                            df_wrk['sun_ra'],
                                            df_wrk['i'],
                                            df_wrk['o_lan']).apply(r2dfct)
        
    ##### Beta computed based on Ecliptic longitude of the sun  
    if calc_beta_sun_eclip_long:
        ### sun_ecliptic_longitude output in RADIANS but NO MODULO !!!
        df_wrk['sun_ecl_long'] = np.mod(pyorbital.astronomy.sun_ecliptic_longitude(df_wrk['epoch']),np.pi*2)
        df_wrk['beta'+b2] = beta_sun_eclip_long(df_wrk['sun_ecl_long'],
                                                df_wrk['o_lan'],
                                                df_wrk['i'],
                                                np.deg2rad(23.45)).apply(r2dfct)
    df_out = DFOrb_in.copy()
    df_out['beta'] = df_wrk['beta'+b1]

    return df_out, df_wrk



def OrbDF_lagrange_interpolate(DForb_in,Titrp,n=10,
                               append_to_input_DF = False,
                               plot=False):
    """
    High level function to interpolate an orbit DataFrame

    Parameters
    ----------
    DForb_in : DataFrame
        an Orbit DataFrame.
    Titrp : iterable of datetime
        Epochs of the wished points.
    n : int, optional
        degree of the polynom. Better if even. The default is 10.
    append_to_input_DF : bool, optional
        append the interpolated DF to the input DF. The default is False.
    plot : bool, optional
        Plot the values. For debug only. The default is False.

    Returns
    -------
    DForb_out : DataFrame
        Interpolated orbits.
        
    Tips
    ----
    Use conv.dt_range to generate the wished epochs range

    """
    DForb_stk = []
    
    for sat,ac in itertools.product(DForb_in.sat.unique(),DForb_in['ac'].unique()):
        
        log.info("process",ac,sat)
        
        DForb_use = DForb_in[(DForb_in.sat == sat) & (DForb_in['ac'] == ac)].copy()
            
        Tdata = DForb_use.epoch.dt.to_pydatetime()
        
        Xitrp = stats.lagrange_interpolate(Tdata,DForb_use.x,Titrp,n=n)
        Yitrp = stats.lagrange_interpolate(Tdata,DForb_use.y,Titrp,n=n)
        Zitrp = stats.lagrange_interpolate(Tdata,DForb_use.z,Titrp,n=n)
        
        ClkDummy = np.array([999999.999999] * len(Titrp))
        
        ARR = np.column_stack((Titrp,Xitrp,Yitrp,Zitrp,ClkDummy))
        
        DForb_tmp = pd.DataFrame(ARR,columns=["epoch","x","y","z","clk"])
        
        ### sometihng else must be tested o give the annex val directly in the col of DForb_tmp
        DFannex_vals = DForb_use.drop(["epoch","x","y","z","clk"],axis=1).drop_duplicates()
        DFannex_vals = pd.concat([DFannex_vals]*(len(Titrp)),ignore_index=True,axis=0)
        
        DForb_tmp = pd.concat((DForb_tmp,DFannex_vals),axis=1)
        DForb_stk.append(DForb_tmp)
        
        if plot:
            # plt.plot(Tdata,DForb_use.x,'o')
            # plt.plot(Titrp,Xitrp,'.')  
            ## GUS mod 220322
            fig,axr = plt.subplots(1,1,sharex='all')
            Symb = axr.plot(Tdata,DForb_use.x,'o')
            Symb = axr.plot(Titrp,Xitrp,'.')

        
    
    DForb_out = pd.concat(DForb_stk)
    
    if append_to_input_DF:
        DForb_out = pd.concat((DForb_in,DForb_out))
        
    DForb_out.reset_index(drop=True)
    DForb_out[["x","y","z","clk"]] = DForb_out[["x","y","z","clk"]].astype(float)
    return DForb_out

def OrbDF_crf2trf(DForb_inp,DF_EOP_inp,time_scale_inp="gps",
                  inv_trf2crf=False):
    """
    Convert an Orbit DataFrame from Celetrial Reference Frame to 
    Terrestrial Reference Frame (.
    
    Requires EOP to work. Cf. note below.

    Parameters
    ----------
    DForb_inp : DataFrame
        Input Orbit DataFrame in Celetrial Reference Frame.
    DF_EOP_inp : DataFrame
        EOP DataFrame  (C04 format).
    time_scale_inp : str, optional
        The time scale used in. manage 'utc', 'tai' and 'gps'.
        The default is "gps".
    inv_trf2crf : bool, optional
        Provide the inverse transformation TRF => CRF.
        The default is False.

    Returns
    -------
    DForb_out : DataFrame
        Output Orbit DataFrame in Terrestrial Reference Frame.
        (or Celestrial if inv_trf2crf is True)
        
    Note
    ----
    The EOP can be obtained from the IERS C04 products.
    e.g.
    https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt
    To get them as a Compatible DataFrame, use the function
    files_rw.read_eop_C04()
    """
    
    DForb = DForb_inp.copy()
    
    import geodezyx.reffram.sofa18 as sofa
    
    ### bring everything to UTC
    if time_scale_inp.lower() == "gps":
        DForb["epoch_utc"] = conv.dt_gpstime2dt_utc(DForb["epoch"])
    elif time_scale_inp.lower() == "tai":
        DForb["epoch_utc"] = conv.dt_tai2dt_utc(DForb["epoch"])
    elif time_scale_inp.lower() == "utc":
        DForb["epoch_utc"] = DForb["epoch"]
    ### TT and UT1 are not implemented (quite unlikely to have them as input)
    
    ### do the time scale's conversion
    DForb["epoch_tai"] = conv.dt_utc2dt_tai(DForb["epoch_utc"])
    DForb["epoch_tt"]  = conv.dt_tai2dt_tt(DForb["epoch_tai"])
    DForb["epoch_ut1"] = conv.dt_utc2dt_ut1_smart(DForb["epoch_utc"],DF_EOP_inp)
        
    ### Do the EOP interpolation 
    DF_EOP_intrp = eop_interpotate(DF_EOP_inp, DForb["epoch_utc"])
    ### bring the EOP to radians
    Xeop = np.deg2rad(conv.arcsec2deg(DF_EOP_intrp['x']))
    Yeop = np.deg2rad(conv.arcsec2deg(DF_EOP_intrp['y']))
    
    TRFstk = []
    
    for tt,ut1,xeop,yeop,x,y,z in zip(DForb["epoch_tt"],
                                      DForb["epoch_ut1"],
                                      Xeop,Yeop,
                                      DForb['x'],DForb['y'],DForb['z']):
    
        MatCRF22TRF = sofa.iau_c2t06a(2400000.5,
                                      conv.dt2MJD(tt),
                                      2400000.5,
                                      conv.dt2MJD(ut1),
                                      xeop,yeop)
        if inv_trf2crf:
            MatCRF22TRF = np.linalg.inv(MatCRF22TRF)
    
        CRF = np.array([x,y,z])
        TRF = np.dot(MatCRF22TRF,CRF)
    
        TRFstk.append(TRF)
    
    ### Final stack and replacement
    TRFall = np.vstack(TRFstk)
    DForb_out = DForb_inp.copy()
    DForb_out[["x","y","z"]] = TRFall
    
    return DForb_out
                                                                     

#### FCT DEF
def OrbDF_reg_2_multidx(OrbDFin,index_order=["prn","epoch"]):
    """
    From an regular Orbit DF generated by read_sp3(), set some columns 
    (typically ["prn","epoch"]) as indexes
    The outputed DF is then a Multi-index DF
    """
    OrbDFwrk = OrbDFin.reset_index()
    OrbDFwrk = OrbDFwrk.sort_values(index_order)
    OrbDFwrk = OrbDFwrk.set_index(index_order,inplace=False)
    return OrbDFwrk

def OrbDF_multidx_2_reg(OrbDFin,index_order=["prn","epoch"]):
    """
    Convert a Multi-index formatted OrbitDF to his original form
    """
    OrbDFwrk = OrbDFin.reset_index()
    OrbDFwrk = OrbDFwrk.sort_values(index_order)
    OrbDFwrk["sys"] = OrbDFwrk["prn"].apply(lambda x: x[0])
    OrbDFwrk["prni"] = OrbDFwrk["prn"].apply(lambda x: int(x[1:]))
    return OrbDFwrk

def OrbDF_common_epoch_finder(OrbDFa_in,OrbDFb_in,return_index=False,
                          supplementary_sort=False,order=["prn","epoch"],
                          skip_reg2multidx_OrbDFa=False,
                          skip_reg2multidx_OrbDFb=False):
    """
    This function finds common satellites and epochs in two Orbit DataFrames and outputs the corresponding Orbit DataFrames.

    Parameters
    ----------
    OrbDFa_in : DataFrame
        The first input Orbit DataFrame.
    OrbDFb_in : DataFrame
        The second input Orbit DataFrame.
    return_index : bool, optional
        If True, the function also returns the common index. Default is False.
    supplementary_sort : bool, optional
        If True, an additional sort is performed. This is useful for multi GNSS where the output DataFrame may not be well sorted. Default is False.
    order : list of str, optional
        The order of the index for the multi-index DataFrame. Default is ["prn","epoch"].
    skip_reg2multidx_OrbDFa : bool, optional
        If True, skips the conversion of the first input DataFrame to a multi-index DataFrame. Default is False.
        The inputs are assumed to be already in multi-index format to optimize execution speed.
        (For advanced use only)
    skip_reg2multidx_OrbDFb : bool, optional
        If True, skips the conversion of the second input DataFrame to a multi-index DataFrame. Default is False.
        The inputs are assumed to be already in multi-index format to optimize execution speed.
        (For advanced use only)

    Returns
    -------
    OrbDFa_out : DataFrame
        The first output Orbit DataFrame with common satellites and epochs.
    OrbDFb_out : DataFrame
        The second output Orbit DataFrame with common satellites and epochs.
    Iinter : Index, optional
        The common index. Only returned if return_index is True.

    Note
    ----
    designed for orbits/sp3 first with sat and epoch as order parmeter,
    but can be used also for instance for snx files with
    STAT and epoch as order parmeter
    """
    if not skip_reg2multidx_OrbDFa:
        OrbDFa = OrbDF_reg_2_multidx(OrbDFa_in,index_order = order)
    else:
        OrbDFa = OrbDFa_in
    if not skip_reg2multidx_OrbDFb:
        OrbDFb = OrbDF_reg_2_multidx(OrbDFb_in,index_order = order)
    else:
        OrbDFb = OrbDFb_in

    I1 = OrbDFa.index
    I2 = OrbDFb.index

    Iinter = I1.intersection(I2)
    Iinter = Iinter.sort_values()

    OrbDFa_out = OrbDFa.loc[Iinter]
    OrbDFb_out = OrbDFb.loc[Iinter]

    if supplementary_sort:
        OrbDFa_out = OrbDFa_out.sort_values(order)
        OrbDFb_out = OrbDFb_out.sort_values(order)

    if len(OrbDFa_out) != len(OrbDFb_out):
        log.warning("len(Orb/ClkDFa_out) != len(Orb/ClkDFb_out)")
        log.warning("TIPS : ClkDFa_in and/or ClkDFb_in might contain duplicates")

    if return_index:
        return OrbDFa_out , OrbDFb_out , Iinter
    else:
        return OrbDFa_out , OrbDFb_out


def OrbDF_const_sv_columns_maker(OrbDFin,inplace=True):
    """
    (re)generate the const and sv columns from the sat one
    """
    if inplace:
        OrbDFin['sys'] = OrbDFin['prn'].str[0]
        OrbDFin['prni']    = OrbDFin['prn'].apply(lambda x: int(x[1:]))
        return None
    else:
        OrbDFout = OrbDFin.copy()
        OrbDFout['sys'] = OrbDFout['prn'].str[0]
        OrbDFout['prni']    = OrbDFout['prn'].apply(lambda x: int(x[1:]))
        return OrbDFout

 #   _____ _            _      _____        _        ______                              
 #  / ____| |          | |    |  __ \      | |      |  ____|                             
 # | |    | | ___   ___| | __ | |  | | __ _| |_ __ _| |__ _ __ __ _ _ __ ___   ___  ___  
 # | |    | |/ _ \ / __| |/ / | |  | |/ _` | __/ _` |  __| '__/ _` | '_ ` _ \ / _ \/ __| 
 # | |____| | (_) | (__|   <  | |__| | (_| | || (_| | |  | | | (_| | | | | | |  __/\__ \ 
 #  \_____|_|\___/ \___|_|\_\ |_____/ \__,_|\__\__,_|_|  |_|  \__,_|_| |_| |_|\___||___/ 
                                                                                                                                                                        
### Clock DataFrames

def ClkDF_filter(ClkDF_in,
             typ=("AS","AR"),
             name=None,
             ac=None,
             epoch_strt=dt.datetime(1980,1,1),
             epoch_end=dt.datetime(2099,1,1),
             name_regex=False):
    """
    Filter the content of a Clock DataFrame

    Parameters
    ----------
    ClkDF_in : DataFrame
        Input Clock DataFrame 
        (a concatenation of DF generated by files_rw.read_clk.
    typ : iterable of str, optional
        List of the types of clocks: AS (satellite) or AR (receiver).
        The default is ("AS","AR").
    name : iterable of str, optional
        List of wished satellites/stations.
        Can be a regex (see also name_regex)
        The default is None.
    ac : iterable of str, optional
        List of wished ACs. The default is None.
    epoch_strt : datetime, optional
        Start epoch. The default is dt.datetime(1980,1,1).
    epoch_end : datetime, optional
        End epoch (not included). The default is dt.datetime(2099,1,1).
    name_regex : bool, optional
        the given names as 'name' arguments are regular expressions
        Some useful regex are given bellow
        The default is False

    Returns
    -------
    Clock DataFrame
        Output Clock DataFrame.
        
    Notes
    -----
    '^E[0-9]{2}': Galileo Satellites
    '^G[0-9]{2}': GPS Satellites

    """
    
    if type(ClkDF_in) is str:
        ClkDF_wrk = utils.pickle_loader(ClkDF_in)
    else:
        ClkDF_wrk = ClkDF_in
    
    BOOL = np.ones(len(ClkDF_wrk)).astype(bool)
    
    if typ:
        BOOLtmp = ClkDF_wrk.type.isin(typ)
        BOOL    = BOOL & np.array(BOOLtmp)

    if name:
        if not name_regex: ### full name mode
            BOOLtmp = ClkDF_wrk.name.isin(name)
            BOOL    = BOOL & np.array(BOOLtmp)   
        else: ### REGEX mode
            BOOLtmp = np.zeros(len(ClkDF_wrk.name)).astype(bool)
            for rgx in name:
                NamSerie = ClkDF_wrk.name
                BOOLtmp = BOOLtmp | np.array(NamSerie.str.contains(rgx))

            BOOL = BOOL & np.array(BOOLtmp)
                 
    if ac:
        BOOLtmp = ClkDF_wrk.ac.isin(ac)
        BOOL    = BOOL & np.array(BOOLtmp)    
        
    ##epoch
    BOOLtmp = (epoch_strt <= ClkDF_wrk.epoch) & (ClkDF_wrk.epoch < epoch_end)
    BOOL    = BOOL & np.array(BOOLtmp)    
    
    return ClkDF_wrk[BOOL]

def ClkDF_filter2(ClkDF_in,
             typ=("AS","AR"),
             name=None,
             ac=None,
             epoch_strt=dt.datetime(1980,1,1),
             epoch_end=dt.datetime(2099,1,1),
             name_regex=False):
    """
    attempt for a faster version of ClkDF_filter, but the original is faster
    """
    
    if type(ClkDF_in) is str:
        ClkDF_wrk = utils.pickle_loader(ClkDF_in)
    else:
        ClkDF_wrk = ClkDF_in
        
    clkdf_stk = []
    
    for (ityp, iname, iac), clkdf_grp in ClkDF_wrk.groupby(["type","name","ac"]):
        if typ:
            bool_typ = True if ityp in typ else False
        else:
            bool_typ = True
        
        if name:
            if not name_regex:
                bool_name = True if iname in name else False
            else:
                bool_name = any([re.search(n, iname) for n in name])
        else:
            bool_name = True
            
        if ac:
            bool_ac = True if iac in ac else False
        else:
            bool_ac = True
            
        
        if not (bool_typ and bool_name and bool_ac):
            continue
        else:
            if epoch_strt > dt.datetime(1980,1,1) or epoch_end < dt.datetime(2099,1,1):
                bool_epoc = (epoch_strt <= clkdf_grp["epoch"]) & (clkdf_grp["epoch"] < epoch_end)
                clkdf_stk.append(clkdf_grp[bool_epoc])
            else:
                clkdf_stk.append(clkdf_grp)

                
            
    clkdf_out = pd.concat(clkdf_stk)
    
    return clkdf_out
                
            


def ClkDF_reg_2_multidx(ClkDFin,index_order=["name","epoch"]):
    """
    From an regular Clock DF generated by read_clk(), set some columns 
    (typically ["name","epoch"]) as indexes
    The outputed DF is then a Multi-index DF
    
    It an adapted version of OrbDF_reg_2_multidx
    """
    
    return OrbDF_reg_2_multidx(ClkDFin,index_order)
    

def ClkDF_common_epoch_finder(ClkDFa_in,ClkDFb_in,return_index=False,
                              supplementary_sort=False,
                              order=["name","epoch"]):
    """
    Find common sats/station and epochs in to Clock DF, and output the
    corresponding Clock DFs
    
    Is an adapted version of OrbDF_common_epoch_finder
    """
    
    return OrbDF_common_epoch_finder(ClkDFa_in,ClkDFb_in,
                                     return_index=return_index,
                                     supplementary_sort=supplementary_sort,
                                     order=order)



def ClkDF_common_epoch_finder_multi(ClkDF_list_in,
                                    return_index=False,
                                    supplementary_sort=False,
                                    order=["name","epoch"]):
    
    """
    Find common sats/station and epochs in to Clock DF, and output the
    corresponding Clock DFs
    
    Is is the multi version of ClkDF_common_epoch_finder
    """
    
    ClkDFref = ClkDF_list_in[0]
    
    #### First loop: we find the common epochs
    for ClkDF in ClkDF_list_in[1:]:
        
        OUTTUP = OrbDF_common_epoch_finder(ClkDFref,ClkDF,
                                           return_index=True,
                                           supplementary_sort=supplementary_sort,
                                           order=order)
        
        ClkDFref , _ , Iinter = OUTTUP
        
    #### second loop: we use the common epochs found for the outputed ClkDF  
    ClkDF_list_out= []
    for ClkDF in ClkDF_list_in:
        ClkDFout = ClkDF.set_index(order).loc[Iinter]
        ClkDF_list_out.append(ClkDFout)
    
    if not return_index:
        return ClkDF_list_out
    else:
        return ClkDF_list_out,Iinter
        


 #   _____ _      _____   __      __   _ _     _       _   _             
 #  / ____| |    |  __ \  \ \    / /  | (_)   | |     | | (_)            
 # | (___ | |    | |__) |  \ \  / /_ _| |_  __| | __ _| |_ _  ___  _ __  
 #  \___ \| |    |  _  /    \ \/ / _` | | |/ _` |/ _` | __| |/ _ \| '_ \ 
 #  ____) | |____| | \ \     \  / (_| | | | (_| | (_| | |_| | (_) | | | |
 # |_____/|______|_|  \_\     \/ \__,_|_|_|\__,_|\__,_|\__|_|\___/|_| |_|
                                                                       


def svn_prn_equiv_DF(path_meta_snx):
    """
    generate a SVN <> PRN equivalent DataFrame

    Parameters
    ----------
    path_meta_snx : str
        path of the MGEX metadata sinex.
        last version avaiable here
        http://mgex.igs.org/IGS_MGEX_Metadata.php

    Returns
    -------
    DFfin : Pandas DataFrame
        SVN <> PRN equivalent DataFrame.

    """
    
    DFsvn  = files_rw.read_sinex_versatile(path_meta_snx,"SATELLITE/IDENTIFIER",
                                           header_line_idx=-2)    

    DFprn  = files_rw.read_sinex_versatile(path_meta_snx,"SATELLITE/PRN",
                                           header_line_idx=-2)
    
    DFsvn.drop(columns='Comment__________________________________',inplace=True)
    DFprn.drop(columns='Comment_________________________________',inplace=True)

    ## the next lines 1 and 3 seems like they have became useless
    DFsvn["SVN_"] = DFsvn["SVN_"].apply(lambda x:x[0] + x[1:])
    DFprn.replace(dt.datetime(1970,1,1),dt.datetime(2099,1,1),inplace=True)
    DFprn["SVN_"] = DFprn["SVN_"].apply(lambda x:x[0] + x[1:])
    
    
    DFstk = []
    
    for isat , sat in DFprn.iterrows():
        svn = sat["SVN_"]
        
        sat["Block"] = DFsvn[DFsvn["SVN_"] == svn]["Block__________"].values[0]
        DFstk.append(sat)
        
    DFfin = pd.concat(DFstk,axis=1).transpose()
    
    DFfin.rename(columns={"SVN_":"SVN",
                          "Valid_From____":"start",
                          "Valid_To______":"end"},inplace=True)
    
    
    DFfin["const"]   = DFfin["SVN"].apply(lambda x:x[0])
    DFfin["SVN_int"] = DFfin["SVN"].apply(lambda x:int(x[1:]))
    DFfin["PRN_int"] = DFfin["PRN"].apply(lambda x:int(x[1:]))    
    
    return DFfin
    

def svn_prn_equiv(sat_in,date_in,
                  svn_prn_equiv_DF,
                  mode="svn2prn",
                  full_output=False):
    """
    Get the equivalence SVN <> PRN for a given epoch
    
    Parameters
    ----------
    sat_in : str
        Satellite "ID", SVN or PRN.
    date_in : datetime
        wished epoch.
    svn_prn_equiv_DF : DataFrame
        Equivalence table generated by svn_prn_equiv_DF.
    mode : str, optional
        prn2svn: PRN > SVN
        svn2prn: SVN > PRN.
        The default is "svn2prn".
    full_output : bool, optional
        get the complete Equivalence table row. The default is False.

    Returns
    -------
    str or DataFrame
    """

    svnorprn1 = mode[:3].upper()    
    svnorprn2 = mode[-3:].upper()
    
    DFsat = svn_prn_equiv_DF[svn_prn_equiv_DF[svnorprn1] == sat_in]
    Bool_date = np.logical_and((DFsat.start <= date_in) , (date_in < DFsat.end))
    DFout = DFsat[Bool_date]
    
    if len(DFout) != 1:
        log.warning("several or no %s entries !!! %s %s",mode,sat_in,date_in)
        
    if full_output:
        return DFout
    else:
        return DFout[svnorprn2].values[0]
    
    
def get_block_svn(sat_in,svn_prn_equiv_DF):
    """
    Get the equivalence SVN block type
    
    Parameters
    ----------
    sat_in : str
        Satellite SVN.
    svn_prn_equiv_DF : DataFrame
        Equivalence table generated by svn_prn_equiv_DF.
    Returns
    -------
    str with the block name
    """
    DFsat = svn_prn_equiv_DF[svn_prn_equiv_DF['SVN'] == sat_in]
    # if DFsat.empty:
    #     print('SVN NOT FOUND')
    # else:        
    block = DFsat.Block.values[0]

    return block    


def stats_slr(DFin,grpby_keys = ['sat'],
              threshold = .5):
    """
    computes statistics for SLR Residuals
    
    Parameters
    ----------
    DFin : Pandas DataFrame
        Input residual Dataframe from read_pdm_res_slr.
    grpby_keys : list of str, optional
        The default is ['sat'].
        per day, per solution, per satellite: ['day','sol','sat']
        per day, per solution, per station: ['day','sol','sta']
        per day, per solution, per satellite, per station: ['day','sol','sta','sat']
    threshold : float
        apply a Threshold

    Returns
    -------
    DD : Output statistics DataFrame
        return the mean, the rms and the std.
    """
    
    DD = DFin[np.abs(DFin["res"]) < threshold]
    
    DD_grp  = DD.groupby(grpby_keys)
    DD_mean = DD_grp['res'].agg(np.mean).rename('mean') * 1000
    DD_rms  = DD_grp['res'].agg(stats.rms_mean).rename('rms')   * 1000
    DD_std  = DD_grp['res'].agg(np.std).rename('std')   * 1000
    DD = pd.concat([DD_mean,DD_std,DD_rms],axis=1)
    DD.reset_index(inplace = True)
    
    return DD    


 #  ______           _   _        ____       _            _        _   _               _____                               _                
 # |  ____|         | | | |      / __ \     (_)          | |      | | (_)             |  __ \                             | |               
 # | |__   __ _ _ __| |_| |__   | |  | |_ __ _  ___ _ __ | |_ __ _| |_ _  ___  _ __   | |__) |_ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___ 
 # |  __| / _` | '__| __| '_ \  | |  | | '__| |/ _ \ '_ \| __/ _` | __| |/ _ \| '_ \  |  ___/ _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
 # | |___| (_| | |  | |_| | | | | |__| | |  | |  __/ | | | || (_| | |_| | (_) | | | | | |  | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
 # |______\__,_|_|   \__|_| |_|  \____/|_|  |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_| |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/
                                                                                                                                         

### EOP / Earth Oreintation Parameters

def eop_interpotate(DF_EOP,Epochs_intrp,eop_params = ["x","y"]):
    """
    Interopolate the EOP provided in a C04-like DataFrame

    Parameters
    ----------
    DF_EOP : DataFrame
        Input EOP DataFrame (C04 format).
        Can be generated by files_rw.read_eop_C04
    Epochs_intrp : datetime of list of datetimes
        Wished epochs for the interpolation.
    eop_params : list of str, optional
        Wished EOP parameter to be interpolated.
        The default is ["x","y"].

    Returns
    -------
    OUT : DataFrame or Series
        Interpolated parameters.
        Series if onely one epoch is provided, DF_EOP elsewere
    """
    if not utils.is_iterable(Epochs_intrp):
        singleton = True
    else:
        singleton = False
    
    I_eop   = dict()
    Out_eop = dict()
    Out_eop["epoch"] = Epochs_intrp    
    
    for eoppar in eop_params:
        I = conv.interp1d_time(DF_EOP.epoch,DF_EOP[eoppar])
        I_eop[eoppar] = I
        try:
            Out_eop[eoppar] = I(Epochs_intrp)
        except ValueError as err:
            log.error("in EOP interpolation")
            log.error("param.: %s, epoch: %s",eoppar,Epochs_intrp)
            raise err
      
    if not singleton:
        OUT = pd.DataFrame(Out_eop)
    else:
        OUT = pd.Series(Out_eop)
        
    return OUT
