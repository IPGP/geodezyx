#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.operational generates RINEX timeline plot. 

it can be imported directly with:
from geodezyx import operational

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
#### Import the logger
import logging
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import tabulate

#### geodeZYX modules
from geodezyx import conv
from geodezyx import operational
from geodezyx import utils

log = logging.getLogger('geodezyx')

##########  END IMPORT  ##########


def rinex_lister(path, add_long_names=True):
    """
    find all rinex in a folder and his subfolders
    path can be a string or a tuple of string => manage multi paths :)

    is very similar with softs_runner.multi_finder_rinex,
    gins_runner.get_rinex_list and operational.rinex_finder

    Parameters
    ----------
    path : str
        archive path.
    add_long_names : bool, optional
        consider new names. The default is True.

    Returns
    -------
    rinexfilelist : list
        list of rinex files.

    Notes
    -----
    operational.rinex_finder must be used in priority !!! (July 2022)
    """

    log.warning("rinex_lister depreciated, use operational.rinex_finder instead!!")

    if type(path) is str:
        path = [path]

    paths_walked = []
    for p in path:
        paths_walked = paths_walked + list(os.walk(p))

    wholefilelist = []
    for tup in paths_walked:
        wholefilelist = wholefilelist + tup[-1]

    wholefilelist = list(set(wholefilelist))

    rinexfilelist = [
        fil
        for fil in wholefilelist
        if re.search(conv.rinex_regex(), os.path.basename(fil))
    ]
    if add_long_names:
        rinexfilelist_new_names = [
            fil
            for fil in wholefilelist
            if re.search(conv.rinex_regex_long_name(), os.path.basename(fil))
        ]
        rinexfilelist = rinexfilelist + rinexfilelist_new_names

    log.info("%s RINEXs found", len(rinexfilelist))

    return rinexfilelist


def rinex_timeline(
    inputlist_or_paths,
    start=dt.datetime(1980, 1, 1),
    end=dt.datetime(2099, 1, 1),
    use_rinex_lister=True,
    dots_plot=False,
    jul_date_plot=False,
    return_figure=False,
    xlim_start_end=False,
):
    """
    Frontend function to plot the aviable RINEXs

    Parameters
    ----------
    inputlist_or_paths : iterable
        list of rinex file paths or list of archive directories.
    start : datetime, optional
        start of the plot time span. The default is dt.datetime(1980,1,1).
    end : datetime, optional
        end of the plot time span. The default is dt.datetime(2099,1,1).
    use_rinex_lister : bool, optional
        if True, inputlist_or_paths is a list of archive directories.
        if False, inputlist_or_paths is a list of  rinex file paths.
        The default is True.
    dots_plot : bool, optional
        print dots instead of line.
        This mode should be avoided because it is slower.
        The default is False.
    jul_date_plot : bool, optional
        Plot the date in julian dates. The default is False.
    return_figure : bool, optional
        returns figure. The default is False.
    xlim_start_end : bool, optional
        Force start and end values to be the x axis limit
        The default is False.

    Returns
    -------
    datadico : dict
        datadico.
    fig : matplotlib figure
        figure

    """

    datadico = rinex_timeline_datadico(
        inputlist_or_paths, use_rinex_lister=use_rinex_lister
    )

    fig = timeline_plotter(
        datadico,
        start=start,
        end=end,
        dots_plot=dots_plot,
        jul_date_plot=jul_date_plot,
        xlim_start_end=xlim_start_end,
    )

    if not return_figure:
        return datadico
    else:
        return datadico, fig


def rinex_timeline_datadico(
    inputlist_or_paths, use_rinex_lister=True, optional_info=""
):
    """
    Generate a RINEX datadico

    Parameters
    ----------
    inputlist_or_paths : iterable
        list of rinex file paths or list of archive directories.
    use_rinex_lister : bool, optional
        if True, inputlist_or_paths is a list of archive directories.
        if False, inputlist_or_paths is a list of  rinex file paths.
        The default is True.
    optional_info : str, optional
        A addtional information for the rinexs found.
        Usually it is the archive name
        The default is ''.

    Returns
    -------
    datadico : dict
        datadico[stat] = [(rinexname1,optional1,date1) ... (rinexnameN,optionalN,dateN)]
    """

    if use_rinex_lister:
        filelist = rinex_lister(inputlist_or_paths)
    else:
        filelist = inputlist_or_paths

    # rinexfilelist = [fil for fil in filelist if re.search( '.*' + conv.rinex_regex() + '$', fil)]
    rinexfilelist_old = [fil for fil in filelist if re.search(conv.rinex_regex(), fil)]
    rinexfilelist_new = [
        fil for fil in filelist if re.search(conv.rinex_regex_long_name(), fil)
    ]
    
    rinexfilelist = rinexfilelist_old + rinexfilelist_new

    if not use_rinex_lister:
        rinexfilelist = [os.path.basename(e) for e in rinexfilelist]

    log.info("%s RINEXs will be ploted on the timeline", len(rinexfilelist))

    statname_lis = sorted(list(set([rin[0:4].lower() for rin in rinexfilelist])))

    log.info("%s stations will be ploted on the timeline", len(statname_lis))

    datadico = dict()

    for stat in statname_lis:
        datadico[stat] = []

    for rnx in rinexfilelist:
        try:
            datadico[rnx[0:4].lower()].append((rnx, optional_info, conv.rinexname2dt(rnx)))
        except:
            log.error("error with : %s", rnx)

    return datadico


def rinex_timeline_datadico_merge(datadico_list, priority_list=None):
    """
    Merge different RINEXs datadico, produced by rinex_timeline_datadico
    coming from different archives

    Parameters
    ----------
    datadico_list : list of dict
        list of RINEX datadico.
    priority_list : list, optional
        priority list of 'optional_info' (archive ID)
        it will erase optional_info of lower priority
        NB : it is not very useful, just sort datadico_list in the right order ....
        The default is None.

    Returns
    -------
    datadico_out : dict
        a merged datadico
    """

    datadico_out = dict()

    datadico_merged = utils.dicts_of_list_merge(*datadico_list)

    for k, dataval in datadico_merged.items():

        rnxname_list = [e[0] for e in dataval]
        archive_list = [e[1] for e in dataval]
        date_list = [e[-1] for e in dataval]

        out_date_list, out_all_list = [], []
        for r, a, d in zip(rnxname_list, archive_list, date_list):
            if d not in out_date_list:
                out_date_list.append(d)
                out_all_list.append((r, a, d))
            elif priority_list:
                ind_existing = out_date_list.index(d)
                archd_existing = out_all_list[ind_existing][1]
                if priority_list.index(a) < priority_list.index(archd_existing):
                    out_date_list.remove(d)
                    out_all_list.remove(out_all_list[ind_existing])
                    out_date_list.append(d)
                    out_all_list.append((r, a, d))
            else:
                pass

        datadico_out[k] = out_all_list

    return datadico_out


def timeline_plotter(
    datadico,
    start=dt.datetime(1980, 1, 1),
    end=dt.datetime(2099, 1, 1),
    dots_plot=False,
    jul_date_plot=False,
    datadico_anex_list=[],
    use_only_stats_of_main_datadico=False,
    colordico_for_main_datadico=None,
    xlim_start_end=False,
    stats_only_list=[],
):
    """


    Parameters
    ----------
    datadico : dict
        a RINEX datadico (see rinex_timeline_datadico for details).
    start : datetime, optional
        start of the plot time span. The default is dt.datetime(1980,1,1).
    end : datetime, optional
        end of the plot time span. The default is dt.datetime(2099,1,1).
    dots_plot : TYPE, optional
        print dots instead of line.
        This mode should be avoided because it is slower.
        The default is False.
    jul_date_plot : bool, optional
        Plot the date in julian dates. The default is False.
    datadico_anex_list : list, optional
        A list of secondary datadicos. The default is [].
    use_only_stats_of_main_datadico : bool, optional
        Use only stats of main datadico.
        Advanced usage.
        The default is False.
    colordico_for_main_datadico : dict, optional
        Color for the main datadico.
        Advanced usage.
        The default is None.
    xlim_start_end : bool, optional
        Force start and end values to be the x axis limit
        The default is False.
    stats_only_list : list, optional
        If given, keep only the stations of this list
        The default is [].
        

    Returns
    -------
    fig : matplotlib figure
        figure

    Notes
    -----
    A simpler version has been commited to geodezyx toolbox GitHub for archive
    on 20180118 15:59A

    """

    fig, ax = plt.subplots()
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)

    if not use_only_stats_of_main_datadico:
        stats_concat_list = list(datadico.keys()) + sum(
            [list(e.keys()) for e in datadico_anex_list], []
        )
        stats_concat_list = list(reversed(sorted(list(set(stats_concat_list)))))
    else:
        stats_concat_list = list(reversed(sorted(list(datadico.keys()))))
        
    stats_concat_list = [e.lower() for e in stats_concat_list]
        
    if stats_only_list:
        stats_only_list = [e.lower() for e in stats_only_list]
        stats_concat_list = [e for e in stats_concat_list if e in stats_only_list]
        

    # the plot has not the same behavior if it is the morning or not
    # (rinexs timelines wont be ploted if it is the morning)
    if dt.datetime.now().hour < 12:
        morning_shift = dt.timedelta(days=1)
    else:
        morning_shift = dt.timedelta(days=0)

    legend_list = []  # must be here before the loop
    for i, stat in enumerate(stats_concat_list):
                
        ## exclude stats
        if len(stats_only_list) > 0 and not stat in stats_only_list:
            continue
        
        # PART 1 : PLOT MAIN DATADICO
        if not stat in datadico.keys():
            print("AAAAAA",stat)
            continue

        # T = Time, O = Station name (Observation)
        Torig = [e[-1] for e in datadico[stat]]
        Oorig = [e[1] for e in datadico[stat]]

        # Time windowing
        T, O = [], []
        for t, o in zip(Torig, Oorig):
            if (start - morning_shift) <= t <= end:
                T.append(t)
                O.append(o)

        T, O = utils.sort_binom_list(T, O)

        TMJD = conv.dt2mjd(T)
        TGrpAll = utils.consecutive_groupIt(
            TMJD, True
        )  # Tuples (start,end) of continue period

        for tgrp in TGrpAll:
            color1 = "C0"
            color2 = ""
            extra_archive = False

            ###*** managing colors
            if colordico_for_main_datadico:
                igrpstart = TMJD.index(tgrp[0])
                igrpend = TMJD.index(tgrp[1])
                Ogrp = O[igrpstart : igrpend + 1]
                Tgrp = TMJD[
                    igrpstart : igrpend + 1
                ]  # all dates in the current continue period

                opt_set = list(set(Ogrp))

                if len(opt_set) == 1:  # Regular case : only one archive
                    if opt_set[0] in colordico_for_main_datadico:
                        color1 = colordico_for_main_datadico[opt_set[0]]
                else:  # several archives ... so the line has to be splited in segments
                    extra_archive = True
                    OSubgrp = utils.identical_groupIt(Ogrp)
                    TSubgrp = utils.sublistsIt(Tgrp, [len(e) for e in OSubgrp])

                    Tgrp_plt = [
                        (e[0], e[-1] + 1) for e in TSubgrp
                    ]  # +1 because the end boundary day is not included
                    Ogrp_plt = [list(set(e))[0] for e in OSubgrp]
                    Color_grp_plt = [
                        (
                            colordico_for_main_datadico[e]
                            if e in colordico_for_main_datadico
                            else color2
                        )
                        for e in Ogrp_plt
                    ]
            ###*** End of managing colors

            tgrp = (
                tgrp[0],
                tgrp[1] + 1,
            )  # +1 because the end boundary day is not included
            # must stay there, in case of not colordico_for_main_datadico
            if not jul_date_plot:
                tgrp = conv.mjd2dt(tgrp)
                if extra_archive:
                    Tgrp_plt = [
                        (conv.mjd2dt(e[0]), conv.mjd2dt(e[1])) for e in Tgrp_plt
                    ]

            # PLOT part
            if tgrp[0] == tgrp[1] + dt.timedelta(
                days=1
            ):  # CASE NO PERIOD, ONLY ONE DAY
                ax.plot(tgrp[0], i, "." + color1)
            else:  # CASE PERIOD
                if not extra_archive:  # ONE ARCHIVE REGULAR CASE
                    ax.plot(tgrp, [i] * 2, "-" + color1)
                else:  # SUBCASE "EXTRA" : SEVERAL ARCHIVE
                    for tt, cc in zip(Tgrp_plt, Color_grp_plt):
                        ax.plot(tt, [i] * 2, "-" + cc)

        # PART 2 : PLOT ANEX DATADICO
        for idatadico_anex, datadico_anex in enumerate(datadico_anex_list):
            if stat in list(datadico_anex.keys()):
                T = datadico_anex[stat]
                T = [t for t in T if start <= t <= end]
                T = sorted(T)

                if jul_date_plot:
                    T = conv.mjd2dt(T)

                (pale_blue_dot,) = ax.plot(
                    T, i * np.ones(len(T)), "o", color="skyblue", label="final SNX"
                )
                legend_list = [pale_blue_dot]

    #### LEGEND
    if colordico_for_main_datadico:
        import matplotlib.lines as mlines

        for arcnam, col in colordico_for_main_datadico.items():
            legend_list.append(mlines.Line2D([], [], color=col, label=arcnam))
            plt.legend(handles=legend_list, loc="upper left", ncol=3, columnspacing=1)

    #    ax.set_yticks(np.arange(0,len(plotstat_lis)-1),plotstat_lis)
    plt.yticks(np.arange(0, len(stats_concat_list)), stats_concat_list)
    fig.autofmt_xdate()
    koef = np.sqrt(2) * 1
    fig.set_size_inches(koef * 11.69, koef * len(stats_concat_list) * 0.28)  # 16.53
    ax.set_ylim([-1, len(stats_concat_list) + 1])
    ax.set_xlim([start, end])

    return fig


def rinex_check_epochs_availability(rinex_path_list):
    """
    Gives simple statistics about RINEX avaiability for all the listed stations

    Parameters
    ----------
    rinex_path_list : list
        A list of rinex paths

    Returns
    -------
    T : str
        a table with statistics.

    """

    results_stk = []

    for rinex_path in rinex_path_list:

        rinex_name = os.path.basename(rinex_path)

        QC = operational.teqc_qc(rinex_path)

        if not QC:
            continue

        epoc_all = int(
            utils.egrep_big_string(
                "Poss. # of obs epochs", QC, only_first_occur=True
            ).split()[-1]
        )
        epoc_disp = int(
            utils.egrep_big_string(
                "Epochs w/ observations", QC, only_first_occur=True
            ).split()[-1]
        )

        dt_rnx = conv.rinexname2dt(rinex_name)

        date_str = conv.dt2str(dt_rnx, "%F")

        percentage = (float(epoc_disp) / float(epoc_all)) * 100.0

        results = [rinex_name, date_str, epoc_disp, epoc_all, percentage]

        results_stk.append(results)

    header = ["RINEX", "date", "Avbl.", "Poss.", "%"]
    T = tabulate.tabulate(results_stk, headers=header)

    return T


def listing_gins_timeline(
    path, stat_strt, stat_end, date_strt, date_end, suffix_regex=""
):
    """find all gins listings in a folder and his subfolders
    and plot timeline of the avaiable listings

    stat_strt,stat_end,date_strt,stat_end : where to find
    in the name the statname and the date"""

    fig, ax = plt.subplots()

    aaa = list(os.walk(path))
    wholefilelist = []
    for tup in aaa:
        wholefilelist = wholefilelist + tup[-1]

    wholefilelist = list(set(wholefilelist))
    lifilelist = [
        fil for fil in wholefilelist if re.search(suffix_regex + r".*gins", fil)
    ]

    statname_lis = sorted(list(set([li[stat_strt:stat_end] for li in lifilelist])))

    list(set(statname_lis))

    datadico = dict()

    for stat in statname_lis:
        datadico[stat] = []

    for li in lifilelist:
        try:
            tup = (li, conv.jjul_cnes2dt(li[date_strt:date_end]))
            datadico[li[stat_strt:stat_end]].append(tup)
        except:
            log.error("error with : %s", li)
    plotstat_lis = []
    for i, stat in enumerate(sorted(datadico.keys())):
        log.info("%s %s", i, stat)
        T = [e[-1] for e in datadico[stat]]
        plotstat_lis.append(stat)
        ax.plot(T, i * np.ones(len(T)), ".")

    i_list = np.arange(0, len(plotstat_lis))

    plt.yticks(i_list, plotstat_lis)

    return fig


#  ______                _   _                _____                                         _
# |  ____|              | | (_)              / ____|                                       | |
# | |__ _   _ _ __   ___| |_ _  ___  _ __   | |  __ _ __ __ ___   _____ _   _  __ _ _ __ __| |
# |  __| | | | '_ \ / __| __| |/ _ \| '_ \  | | |_ | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
# | |  | |_| | | | | (__| |_| | (_) | | | | | |__| | | | (_| |\ v /  __/ |_| | (_| | | | (_| |
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \_____|_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
#                                                                        __/ |
#                                                                       |___/


# def rinex_timeline_datadico_merge_not_very_smart(datadico_list, priority_list):
#     """
#     Merge different RINEXs datadico, produced by rinex_timeline_datadico
#     coming from different archives
#     Args :
#         rinex_timeline_datadico : list of RINEX datadico
#         priority_list : priority list of 'optional_info' (archive ID)
#                         it will erase optional_info of lower priority
#     Returns :
#         datadico_out : a merged datadico
#     """

#     datadico_out = dict()

#     datadico_merged = utils.dicts_of_list_merge(*datadico_list)

#     for k, dataval in datadico_merged.items():

#         rnxname_list = [e[0] for e in dataval]
#         archive_list = [e[1] for e in dataval]
#         date_list = [e[-1] for e in dataval]

#         out_date_list, out_all_list = [], []
#         for r, a, d in zip(rnxname_list, archive_list, date_list):
#             if d not in out_date_list:
#                 out_date_list.append(d)
#                 out_all_list.append((r, a, d))
#             else:
#                 ind_existing = out_date_list.index(d)
#                 archd_existing = out_all_list[ind_existing][1]
#                 if priority_list.index(a) < priority_list.index(archd_existing):
#                     out_date_list.remove(d)
#                     out_all_list.remove(out_all_list[ind_existing])
#                     out_date_list.append(d)
#                     out_all_list.append((r, a, d))

#         datadico_out[k] = out_all_list

#     return datadico_out


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
#                 TMJD=dt2mjd(T)
#                 TGrpAll = utils.consecutive_groupIt(TMJD,True)
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
#                             OSubgrp = utils.identical_groupIt(Ogrp)
#                             TSubgrp = utils.sublistsIt(Tgrp , [len(e) for e in OSubgrp])
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
#                         tgrp = mjd2dt(tgrp)
#                         if extra_archive:
#                             Tgrp_plt = [ (mjd2dt(e[0]) , mjd2dt(e[1])) for e in Tgrp_plt ]
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
#                     T = mjd2dt(T)
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
#                            o_major = utils.most_common(Ogrp)
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
