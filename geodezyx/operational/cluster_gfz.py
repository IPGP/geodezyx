#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:55:51 2019

@author: psakicki
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
import numpy as np
import subprocess
import re
import time


#### geodeZYX modules
from geodezyx import utils


#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

##########  END IMPORT  ##########
def cluster_GFZ_run(commands_list,
                    bunch_on_off = True,
                    bunch_job_nbr = 10,
                    bunch_wait_time = 600,
                    bj_check_on_off = True,
                    bj_check_mini_nbr = 2,
                    bj_check_wait_time = 120,
                    bj_check_user="auto"):


    history_file_path = None
    wait_sleeping_before_launch=5

    i_bunch = 0

    if bj_check_user == "auto":
        bj_check_user=utils.get_username()

    log_path = "/home/" + bj_check_user + "/test_tmp.log"
    LOGobj = open(log_path , 'w+')

    print ("****** JOBS THAT WILL BE LAUNCHED ******")
    print ('Number of jobs : ' + str(len(commands_list)))
    print ("****************************************")

    for kommand in commands_list:

        ########## LOG/PRINT command
        print(kommand)
        LOGobj.write(kommand + '\n')

        ########## LOG/PRINT sleep
        info_sleep = "INFO : script is sleeping for " +str(wait_sleeping_before_launch) +"sec (so you can cancel it) "
        print(info_sleep)
        LOGobj.write(info_sleep + '\n')
        time.sleep(wait_sleeping_before_launch)

        ########## LOG/PRINT start
        info_start = "INFO : script starts @ " + str(dt.datetime.now())
        print(info_start)
        LOGobj.write(info_start + '\n')

        ########## RUN command here !!
        p = subprocess.Popen('',executable='/bin/csh', stdin=subprocess.PIPE , stdout=subprocess.PIPE , stderr=subprocess.PIPE)
        stdout,stderr = p.communicate( kommand.encode() )

        if history_file_path:
            with open(history_file_path , "a") as myfile:
                myfile.write(kommand + '\n')

        i_bunch += 1

        ########## Bunch On/Off : check if a bunch of job has been launched
        if bunch_on_off and np.mod(i_bunch , bunch_job_nbr) == 0:
            info_bunch = "INFO : sleeping @ " + str(dt.datetime.now()) + " for " + str(bunch_wait_time) + "s b.c. a bunch of " + str(bunch_job_nbr) + " jobs has been launched"
            print(info_bunch)
            time.sleep(bunch_wait_time)
            LOGobj.write(info_bunch + '\n')

            ########## Bunch On/Off Check : check if the bunch is finished (experimental but on is better)
            ###### THIS PART MUST BE MERGE WITH THE SMALLER FCT BELLOW
            if bj_check_on_off:
                print("INFO : BJ Check : All jobs should be finished now, let's see if there is some latecomers")
                bj_check_tigger = False
                bj_command = "perl /dsk/igs2/soft_wrk/psakicki/SOFT_EPOS8_BIN_TOOLS/SCRIPTS/e8_bjobs_local.pl"
            while not bj_check_tigger:
                bj_list = subprocess.check_output(bj_command,shell='/bin/csh')
                bj_pattern_checked =   bj_check_user + ' *' +  bj_check_user

                bj_list_checked = [re.search(bj_pattern_checked,l) for l in bj_list]
                bj_list_checked_sum = np.sum(bj_list_checked)

                if bj_list_checked_sum > bj_check_mini_nbr:
                    print("INFO : sleeping @ " + str(dt.datetime.now()) + " for " + str(bj_check_wait_time) + "s b.c." + str(bj_list_checked_sum) + "job(s) match pattern " + bj_pattern_checked)
                    print(bj_list)
                    time.sleep(bj_check_wait_time)

                else:
                    bj_check_tigger = True
                    print("INFO : let's continue, no job matchs the pattern " + bj_pattern_checked)

def number_job_user(bj_check_user=None,verbose=True):
    if not bj_check_user:
        bj_check_user=utils.get_username()
        
    bj_command = "perl /dsk/igs2/soft_wrk/psakicki/SOFT_EPOS8_BIN_TOOLS/SCRIPTS/e8_bjobs_local.pl"

    bj_list = subprocess.check_output(bj_command,shell='/bin/csh')
    bj_list = bj_list.decode("utf-8")
    bj_list = bj_list.split("\n")
        
    bj_pattern_checked =   bj_check_user + ' *' +  bj_check_user

    bj_list_checked = [bool(re.search(bj_pattern_checked,l)) for l in bj_list]
    bj_list_checked_sum = np.sum(bj_list_checked)
    
    if verbose:
        print("INFO: ",bj_list_checked_sum,"running jobs found for",bj_check_user)
    
    return bj_list_checked_sum

def sleep_job_user(bj_check_user=None,maxjob=20,bj_check_wait_time=20):
    if not bj_check_user:
        bj_check_user=utils.get_username()
    
    n_job = number_job_user(bj_check_user)
    
    if n_job >= maxjob:
        print("INFO : sleeping @ " + str(dt.datetime.now()) + " for " + str(bj_check_wait_time) + "s b.c." + str(n_job) + " jobs are runing")
        time.sleep(bj_check_wait_time)
        
    return None
    





