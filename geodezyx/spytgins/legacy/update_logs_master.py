"""
See also update_master_station.py
(PS250317)
"""

from MasterStationFile import *
import os
import shutil
from datetime import datetime as dt

### VARIABLES A ADAPTER ###
#SPOTGINS logfiles directory
spotgins_logdir = 'spotgins/metadata/logfiles'
#Local directory with up-to-date logfile
local_logdir = 'S:/data/meta/gpslog'
# AC code
ac = 'ULR'
###########################

#Master station file
master = MasterStationFile.from_file('spotgins/metadata/stations/station_master_file.dat')

#List of updated logfiles
out = open('new_logfiles.list','w')

for (acro9,log) in master.data[master.data['AC'] == ac][['NAME','LOGFILE_NAME']].values:

    if str(log) == 'nan':
        spotgins_logdate = dt(1900,1,1)
    else:
        try:
            spotgins_logdate = dt.strptime(log.split('_')[1].split('.')[0],'%Y%m%d')
        except:
            print (f"Unable to extract log date for {acro9}")
            continue

    #new logfiles?
    for local_log in os.listdir(local_logdir):
        if local_log.lower().startswith(acro9.lower()):
            local_logdate = dt.strptime(local_log.split('_')[1].split('.')[0], '%Y%m%d')

            if local_logdate > spotgins_logdate:
                print (f"New logfile for {acro9}: {log} -> {local_log}")

                #Update in master
                index = master.data[master.data['NAME'] == acro9].index[0]
                master.data.loc[index, 'LOGFILE_NAME'] = local_log

                #Copy in SPOTGINS logdir
                shutil.copyfile(f"{local_logdir}/{local_log}",f"{spotgins_logdir}/{local_log}")

                out.write(f"{acro9}\n")

#writing the updated master station file. Default output: station_master_file.dat.new
master.write()
out.close()





