#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.operational generates RINEX timeline plot. 

it can be imported directly with:
from geodezyx import operational

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU GPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (GFZ, pierre.sakic@gfz-postdam.de)
GitHub repository :
https://github.com/GeodeZYX/GeodeZYX-Toolbox_v4
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import os 
import re
import tabulate

#### geodeZYX modules
from geodezyx import conv
from geodezyx import operational
from geodezyx import utils

##########  END IMPORT  ##########

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
    
    rinexfilelist           = [fil for fil in wholefilelist if re.search( conv.rinex_regex() , os.path.basename(fil))]
    if add_new_names:
        rinexfilelist_new_names = [fil for fil in wholefilelist if re.search( conv.rinex_regex_new_name() , os.path.basename(fil))]
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
#    rinexfilelist = [fil for fil in filelist if re.match( '.*' + conv.rinex_regex() + '$', fil)]
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
#            TGrp = utils.consecutive_groupIt(dt2MJD(T),True)
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

    #rinexfilelist = [fil for fil in filelist if re.search( '.*' + conv.rinex_regex() + '$', fil)]
    rinexfilelist_old = [fil for fil in filelist if re.search( conv.rinex_regex() , fil)]
    rinexfilelist_new = [fil for fil in filelist if re.search( conv.rinex_regex_new_name() , fil)]

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
            datadico[rnx[0:4]].append((rnx,optional_info,conv.rinexname2dt(rnx)))
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

    datadico_merged = utils.dicts_of_list_merge(*datadico_list)

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

    datadico_merged = utils.dicts_of_list_merge(*datadico_list)

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

        T,O = utils.sort_binom_list(T,O)

        TMJD=conv.dt2MJD(T)
        TGrpAll = utils.consecutive_groupIt(TMJD,True) # Tuples (start,end) of continue period

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
                    OSubgrp = utils.identical_groupIt(Ogrp)
                    TSubgrp = utils.sublistsIt(Tgrp , [len(e) for e in OSubgrp])

                    Tgrp_plt      = [ (e[0] , e[-1] + 1)    for e in TSubgrp ] # +1 because the end boundary day is not included
                    Ogrp_plt      = [ list(set(e))[0] for e in OSubgrp     ]
                    Color_grp_plt = [ colordico_for_main_datadico[e] if e in colordico_for_main_datadico else color2 for e in Ogrp_plt ]
            ###*** End of managing colors

            tgrp = (tgrp[0] , tgrp[1] + 1) # +1 because the end boundary day is not included
                                           # must stay there, in case of not colordico_for_main_datadico

            if not jul_date_plot:
                tgrp = conv.MJD2dt(tgrp)
                if extra_archive:
                    Tgrp_plt = [ (conv.MJD2dt(e[0]) , conv.MJD2dt(e[1])) for e in Tgrp_plt ]

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
                    T = conv.MJD2dt(T)

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
            tup = (li , conv.jjulCNES2dt(li[date_strt:date_end] ))
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

        QC = operational.teqc_qc(rinex_path)

        if not QC:
            continue

        epoc_all  = int(utils.egrep_big_string("Poss. # of obs epochs" ,QC,only_first_occur=True).split()[-1])
        epoc_disp = int(utils.egrep_big_string("Epochs w/ observations",QC,only_first_occur=True).split()[-1])

        dt_rnx = conv.rinexname2dt(rinex_name)

        date_str = conv.dt2str(dt_rnx,"%F")

        percentage = (float(epoc_disp) / float(epoc_all)) * 100.

        results = [rinex_name,date_str,epoc_disp,epoc_all,percentage]

        results_stk.append(results)

    header = ['RINEX','date','Avbl.', 'Poss.', '%']
    T = tabulate.tabulate(results_stk,headers=header)

    return T
