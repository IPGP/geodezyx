# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 11:32:10 2015

@author: pierre
"""

import softs_runner
import glob
import itertools
import os
import collections

# MONO FILE MODE
if 1:
    # path of the generic configuration file
    generik_conf = "/home/psakicki/THESE/SOFTWARES/RTKLIB_Linux/conf_files/kine_diff_1.conf"
    #  working directory : where the temp file will be saved and the results
    # produced
    working_dir  = '/home/psakicki/THESE/DATA/1506_GEODESEA/GNSS_GEODESEA/REBOOT_4_CNFG2_16/PROCESSED/RTKLIB'
    # rinex base path
    rnx_base = '/home/psakicki/THESE/DATA/1506_GEODESEA/GNSS_ONSHORE/nica/1sec/nica/nica1730.15o'
    # rinex rover path
    rnx_rover = '/home/psakicki/THESE/DATA/1506_GEODESEA/GNSS_GEODESEA/REBOOT_4_CNFG2_16/RINEX/gspl173A.15d.Z'
    # name of the experience
    exp_prefix   = "GEODESEA_NICA_REBOOT"
    
    # run of the processing
    softs_runner.rtklib_run_from_rinex(rnx_rover,rnx_base,
                                       generik_conf,working_dir, 
                                       experience_prefix=exp_prefix,
                                       outtype='geo',calc_center='igs')

    

#MULTI MODE : severals days to process BUT only ONE base
if 0:
    ENSGrefpos = [4201575.85487 ,     189861.53838  ,   4779065.51962 ]
    XYZbase = ENSGrefpos
    
    # rinex rover directory path & List of the rover rinexs
    rnx_rover_dir = "/home/psakicki/THESE/DATA/1509_TOIT_DPTS/MLVL/mlvl"
    rnx_rover_lis = glob.glob(rnx_rover_dir + '/*')
    
    # rinex rover directory path & List of the rover rinexs
    rnx_base_dir  = "/home/psakicki/THESE/DATA/1509_TOIT_DPTS/CRINEXs_Daily1Hz"
    rnx_base_lis  = glob.glob(rnx_base_dir  + '/*')
    
    #  working directory : where the temp files will be saved and the results
    # produced
    working_dir   = '/home/psakicki/THESE/DATA/1509_TOIT_DPTS/PROCESSING'
    
    # path of the generic configuration file
    generik_conf  = "/home/psakicki/THESE/SOFTWARES/RTKLIB_Linux/conf_files/kine_diff_1.conf"
    
    # name of the experience
    exp_prefix    = "TOITDPTS_ensg_mlvl"
    
    # run of the processing
    for (rnx_base,rnx_rover) in itertools.product(rnx_base_lis,rnx_rover_lis):
        if softs_runner.same_day_rinex_check(rnx_base,rnx_rover):
            
            softs_runner.rtklib_run_from_rinex(rnx_rover,rnx_base,
                                               generik_conf,
                                               working_dir, 
                                               experience_prefix=exp_prefix,
                                               XYZbase=XYZbase,outtype='deg')  
    
#MULTI MODE : severals days to process with several bases
if 0:
    # reference positions of the severals base
    # stocked in a dictionnary
    refposdic     = collections.OrderedDict([('ANGL', (4404491.463891803, -108110.74519999999, 4596491.369008197)), 
                                             ('AUNI', (4429451.702603825, -73344.22005901638 , 4573322.050031694)),
                                             ('BRES', (4370725.68305082 , -36125.112323497255, 4629768.617508197)),
                                             ('CHPH', (4236233.126345355, 110998.20365901639, 4751117.440314207)), 
                                             ('ILDX', (4436670.968227869, -91137.98939453551, 4566018.180908743)), 
                                             ('MAN2', (4274275.839450819, 11584.458500000006, 4718386.110908197)), 
                                             ('ROYA', (4466458.820191803, -79862.86325300546, 4537304.7706551915)),
                                             ('SMNE', (4201791.951862842, 177945.60513551912, 4779286.986814207))])
                                         
    #  working directory : where the temp files will be saved and the results
    # produced
    working_dir       = '/home/psakicki/THESE/DATA/1604_BOUEES/PROCESSING/RTKLIB'
    # name of the experience
    experience_prefix = "RTKLIB_BOUEES_AIX"
    
    # rinex rover directory path & List of the rover rinexs
    rnx_base_dir  = "/home/psakicki/THESE/DATA/1604_BOUEES/AIX/STATIONSSS_RGP/"
    rnx_base_lis  = glob.glob(rnx_base_dir  + '/*o')
    
    # rinex base directory path & List of the base rinexs
    rnx_rover_dir = "/home/psakicki/THESE/DATA/1604_BOUEES/AIX/UMRB/MERGED/"
    rnx_rover_lis = glob.glob(rnx_rover_dir + '/*o')

    # run of the processing
    for (rnx_base,rnx_rover) in itertools.product(rnx_base_lis,rnx_rover_lis):
        if softs_runner.same_day_rinex_check(rnx_base,rnx_rover):
            #the position of the base is seeked in the dictionnary defined
            # above
            base_CODEname = os.path.basename(rnx_base)[0:4].upper()
            refpos = refposdic[base_CODEname]
            print("working with :" ,rnx_base,rnx_rover)
            softs_runner.rtklib_run_from_rinex(rnx_rover,rnx_base,
                                               generik_conf,working_dir, 
                                               experience_prefix=exp_prefix,
                                               XYZbase=refpos,outtype='xyz')
