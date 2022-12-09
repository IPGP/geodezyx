# -*- coding: utf-8 -*-
import datetime as dt

#### Import star style
from geodezyx import operational

# STATDICO 
# a statdico is a dictionary associating Archives Centers to list of stations 
# EX :
# statdico['archive center 1'] = ['STA1','STA2','STA3', ...]
# statdico['archive center 2'] = ['STA2','STA1','STA4', ...]

# EXEMPLEs
#statdico['igs'] = ['scub' , 'guat' , 'ssia' , 'slor' , 'mana' , 'cro1']
#statdico['rgp'] = []
#statdico['ovsg'] = ['abd0']
#statdico['unavco'] = ['CN40' , 'CN19' , 'CN38' , 'SAN0' , 'CN35' , 'ROA0' , 'CN29' , 'CN18' , 'CBMD' , 'LCSB' , 'GCFS' , 'GCEA' , 'CN10' , 'CN12' , 'PRCG' , 'SMRT' , 'BARA']
#statdico['igs'] = ['scub' , 'guat' , 'ssia' , 'slor' , 'mana' , 'cro1']
#statdico['igs'] = ['SFER','WARN','GANP','UZHL','PTBB','TITZ','FFMJ','WSRT','KOSG','DLFT']

# if you want you can download the same stat. in several servers (assuming is the same data)
#statdico['igs'] = ['TLSE']
#statdico['rgp'] = ['TLSE']
#statdico['renag'] = ['TLSE']

#statdico['igs'] = ['LROC']
#statdico['rgp'] = ['LROC']
#statdico['renag'] = ['LROC']

# It can also download broadcasts
# statdico['brdc'] = ['BRDC']
#
# ========= END OF EXEMPLES =========

statdico = dict()
statdico['igs'] = ['scub' , 'guat' , 'ssia' , 'slor' , 'mana' , 'cro1']


# START & END DATES OF DATA TO DOWNLOAD
s,e = conv.start_end_date_easy(2020,300,2020,301)


# PLACE OF THE ARCHIVE WHERE THE RINEX WILL BE SAVED
archive_dir = "/wrk/psakicki/PLAYGROUND/psakicki/0000_PROJECTS_OTHERS/2112_jura21/02_data_rinexs/orpheon"


# ARCHTYPE : structure of the archive like
#archtype ='stat/year/doy'
#archtype ='stat'
#archtype ='week/dow/stat'
#archtype ='stat/year'
#archtype ='year/doy'
from datetime import datetime
parallel_download=1

## login for some servers
user   = 'garcia'
passwd = 'testpass' #(dummy value here ;) )


statdico = dict()
statdico['ens_fr'] = ['abel','mesa']
start = datetime.strptime('19/09/20', '%d/%m/%y')
end = datetime.strptime('20/09/20', '%d/%m/%y')
archtype ='stat/year/doy'
archive_dir = "/dsk/mansur/tst_rnx/"
parallel_download=1

if True:
    urllist,savedirlist = operational.multi_downloader_rinex(statdico,
                                                             archive_dir,
                                                             start,end,
                                                             archtype,
                                                             parallel_download,
                                                             sorted_mode=0,
                                                             user=user,
                                                             passwd=passwd)