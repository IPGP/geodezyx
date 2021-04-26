#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:18:10 2017

@author: adminuser

This script aims to compare two SP3 orbits
"""

from natsort import natsorted, ns

import geoclass as gcls

from tabulate import tabulate
import softs_runner

### SIMPLE COMPARISON
if True:
    p1="/home/psakicki/aaa_FOURBI/wk1989/com19896.sp3"
    p2="/home/psakicki/aaa_FOURBI/wk1989/grm19896.sp3"
    name1=''
    name2=''
    save_plot=True
    save_plot_dir=os.path.dirname(p1)

    RTNoutput         = True # False only for debug & personal interest !!
    convert_ECEF_ECI  = True # False only for debug !!
    clean_null_values = "any" # Can be "all" , "any", False or True (True => "all")

    sats_used_list = ["E"] # G,E,R,C ...

    Diff_all    = gcls.compar_orbit(p1,p2,
                                    sats_used_list=sats_used_list,
                                    RTNoutput=RTNoutput,
                                    convert_ECEF_ECI=convert_ECEF_ECI,
                                    name1=name1,name2=name2,
                                    clean_null_values = clean_null_values,
                                    step_data = 900)

    _           = gcls.compar_orbit_plot(Diff_all,
                                        save_plot=save_plot,
                                        save_plot_dir=save_plot_dir)

    ComparTable = gcls.compar_orbit_table(Diff_all)


### MORE COMPLEX CASE : LOOP FOR SEARCHING SAME ORBIT NAME IN DIFFERENTS FOLDERS
if 0:
    p2_dir_list = ['/home/adminuser/Documents/1709_compar_orbit/data/FINALorbs/FINALmod1',
    '/home/adminuser/Documents/1709_compar_orbit/data/FINALorbs/FINALatx1967ORI',
    '/home/adminuser/Documents/1709_compar_orbit/data/FINALorbs/FINAL',
    '/home/adminuser/Documents/1709_compar_orbit/data/FINALorbs/FINALatx1967MOD',
    '/home/adminuser/Documents/1709_compar_orbit/data/FINALorbs/FINALorig_good']

    orbnam = '/ANA/2017/2017_222/ORB/2017_222_orb_1d.sp3'


    p_ref = '/home/adminuser/Documents/1709_compar_orbit/data/FINALorbs/OFFICIAL/2017_222_orb_1d.sp3'

    for p2_dir in p2_dir_list:
        p2 = p2_dir + orbnam

        D1 = gcls.read_sp3(p_ref)
        D2 = gcls.read_sp3(p2)

        name1 = p_ref.split('/')[7]
        name2 = p2.split('/')[7]

        RTNoutput        = True # False only for debug & personal interest !!
        convert_ECEF_ECI = True # False only for debug !!

        sats_used_list = ['E'] # G,E,R,C ...

        Diff_all    = gcls.compar_orbit(p_ref,p2,
                                    sats_used_list=sats_used_list,
                                    RTNoutput=RTNoutput,
                                    convert_ECEF_ECI=convert_ECEF_ECI,
                                    name1=name1,name2=name2)


        _           = gcls.compar_orbit_plot(Diff_all)
        ComparTable = gcls.compar_orbit_table(Diff_all)
