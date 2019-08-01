# -*- coding: utf-8 -*-
import softs_runner
import datetime as dt

# STATDICO 
# a statdico is a dictionary associating Archives Centers to list of stations 
# EX :
# statdico['archive center 1'] = ['STA1','STA2','STA3', ...]
# statdico['archive center 2'] = ['STA2','STA1','STA4', ...]

# ========= END OF EXEMPLES =========

statdico = dict()
#POUR LA GWADA
statdico['rgp']      = ['ABMF',"LMMF","FFT2","PPTG"]
statdico['igs']      = ['ZIMM',"POTS"]
#statdico['unavco'] = ['CN40' , 'CN19' , 'CN38' , 'SAN0' , 'CN35' , 'ROA0' , 'CN29' , 'CN18' , 'CBMD' , 'LCSB' , 'GCFS' , 'GCEA' , 'CN10' , 'CN12' , 'PRCG' , 'SMRT' , 'BARA']

# START & END DATES OF DATA TO DOWNLOAD
s,e = dt.datetime(2016,1,1) , dt.datetime(2016,12,31)
s,e = softs_runner.start_end_date_easy(2015,1,2018,310)

# PLACE OF THE ARCHIVE WHERE THE RINEX WILL BE SAVED
archive_dir = "/home/psakicki/aaa_FOURBI/DL_RINEX_TEST"

# ARCHTYPE : structure of the archive like
#archtype ='stat/year/doy'
#archtype ='stat'
#archtype ='week/dow/stat'
archtype ='stat/year'

parallel_download=1

user=""
passwd=""


if 1:
    url_list = softs_runner.multi_downloader_rinex(statdico,archive_dir,s,e,archtype,
                                                   parallel_download,sorted_mode=0,
                                                   user=user,passwd=passwd)

# ================= ORBITS ==============================
if 0:
    orblis = softs_runner.multi_downloader_orbs_clks( archive_dir , s , e ,
                                                     archtype='/',
                                                     calc_center = 'igs')
                                               
