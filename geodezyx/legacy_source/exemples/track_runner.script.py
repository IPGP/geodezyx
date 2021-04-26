# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 14:08:14 2015

@author: psakicki
"""

import softs_runner
import glob , itertools , os
import collections

# SINGLE MODE
if 0:
    # rinex rover path
    rnx_rover    = "/home/psakicki/THESE/DATA/1506_TEST_TOIT_3/OPERA_RINEX/ROOF/4h/roof161a.15o"
    # rinex base path
    rnx_base     = "/home/psakicki/THESE/DATA/1506_TEST_TOIT_3/OPERA_RINEX/ILEL/4h/ilel161a.15o"
    
    #  working directory : where the temp file will be saved and the results
    # produced
    working_dir       = '/home/psakicki/THESE/RUNNING_EXP/track/TESTOIT3'
    # name of the experience
    experience_prefix = "TESTOIT3_air"
    
    # XYZ position of the reference
    ILELrefpos        = [ 0.442604443913433E+07 , 
                         -0.894254996071955E+05 ,
                          0.457629676666969E+07 ]

    XYZbase  = ILELrefpos

    # calc mode of track : air , short or long (air is better for kinematics)
    mode     = 'air'
    
    # output : XYZ or FLH (geographic , lat lon h)
    outtype  = 'XYZ'

    # run of the processing
    softs_runner.track_runner(rnx_rover,rnx_base,working_dir,
                              experience_prefix,XYZbase,mode=mode)
        
                  
#MULTI MODE : severals days to process BUT only ONE base
if 0:  
    # XYZ position of the reference                   
    NICArefpos       = [4581808.8798  ,   581032.2302   ,   4384493.0375  ]
    XYZbase  = NICArefpos
    
    # rinex rover directory path & List of the rover rinexs
    rnx_rover_dir = "/home/psakicki/THESE/DATA/1506_GEODESEA/GNSS_GEODESEA/PRINCIPAL_GSP1/RINEX_GSPL/1Hz/splited_1h/"
    rnx_rover_lis = glob.glob(rnx_rover_dir + '/*')

    # rinex base directory path & List of the base rinexs
    rnx_base_dir  = "/home/psakicki/THESE/DATA/1506_GEODESEA/GNSS_ONSHORE/nica/1sec/nica"
    rnx_base_lis  = glob.glob(rnx_base_dir  + '/*')

    #  working directory : where the temp file will be saved and the results
    # produced
    working_dir   = '/home/psakicki/THESE/DATA/1509_TOIT_DPTS/PROCESSING/TRACK'
    # name of the experience
    experience_prefix = "GEODESEA_NICA"
    
    # calc mode of track : air , short or long (air is better for kinematics)
    mode     = 'air'
    # output : XYZ or FLH (geographic , lat lon h)
    outtype  = 'FLH'
    
    # run of the processing
    for (rnx_base,rnx_rover) in itertools.product(rnx_base_lis,rnx_rover_lis):
        if softs_runner.same_day_rinex_check(rnx_base,rnx_rover):
            softs_runner.track_runner(rnx_rover,rnx_base,working_dir,
                                      experience_prefix,XYZbase,mode=mode)

#MULTI MODE : severals days to process with several bases
if 0: 
    # reference positions of the severals base
    # stocked in a dictionnary
    refposdic = collections.OrderedDict([('CHPH', (4236233.08156, 110998.26463599999, 4751117.477176)), 
                                         ('MAN2', (4274275.799144, 11584.521544, 4718386.149148)), 
                                         ('MLVL', (4201576.823332, 189860.284372, 4779064.903872)), 
                                         ('SIRT', (4213550.772464001, 162494.69744, 4769661.8907079995)),
                                         ('SMNE', (4201791.9088, 177945.66576799998, 4779287.023676001))])
    
    # rinex rover directory path & List of the rover rinexs
    rnx_rover_dir = "/home/psakicki/THESE/DATA/1509_TOIT_DPTS/CRINEXs_Daily1Hz/"
    rnx_rover_lis = glob.glob(rnx_rover_dir + '/*249*') + glob.glob(rnx_rover_dir + '/*251*')

    # rinex base directory path & List of the base rinexs
    rnx_base_dir  = "/home/psakicki/THESE/DATA/1509_TOIT_DPTS/STATIONSSS_RGP/telechargement_RGP_*/recherche_1"
    rnx_base_lis  = glob.glob(rnx_base_dir  + '/*o')

    #  working directory : where the temp files will be saved and the results
    # produced
    working_dir       = '/home/psakicki/THESE/DATA/1509_TOIT_DPTS/PROCESSING/TRACK_longmode'
    # name of the experience
    experience_prefix = "TRACK_TOITDPTS1603_longmode"
    
    # calc mode of track : air , short or long (air is better for kinematics)
    mode     = 'air'
    mode     = 'short'
    mode     = 'long'

    # output : XYZ or FLH (geographic , lat lon h)
    outtype  = 'FLH'
    
    # run of the processing
    for (rnx_base,rnx_rover) in itertools.product(rnx_base_lis,rnx_rover_lis):
        if softs_runner.same_day_rinex_check(rnx_base,rnx_rover):            
            #the position of the base is seeked in the dictionnary defined
            # above
            base_CODEname = os.path.basename(rnx_base)[0:4].upper()
            XYZbase = refposdic[base_CODEname]
            print("working with :" ,rnx_base,rnx_rover)
            softs_runner.track_runner(rnx_rover,rnx_base,working_dir,
                                      experience_prefix,XYZbase,mode=mode,
                                      calc_center='igs',XYZbase=XYZbase,
                                      XYZrover=XYZrover)
                              


                             
            