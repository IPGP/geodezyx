"""
See also create_stationfile_sitelogs.py which is 100% Python and does not use bash commands.
(PS2503xx)
"""

from StationFile import *
from MasterStationFile import *
import argparse
import pandas as pd
import subprocess
from pathlib import Path
import os

#Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-list', dest='stalist',help='list of 9char acronyms of the stations to process',required=False)
args = parser.parse_args()

sf = 'spotgins/metadata/stations/station_file.dat'
ma = 'spotgins/metadata/stations/station_master_file.dat'

master = MasterStationFile.from_file(ma)
sta = StationFile.from_file(sf,ma)

miss = open('notinmaster.list','w')

datalist = open(args.stalist)
stalist = []
for line in datalist:
	acro9 = line.strip().split()[0]
	print (acro9)

	#in master?
	inmaster = master.data[master.data['NAME'] == acro9]
	if len(inmaster) == 0:
		print ('WARNING: ' + acro9 + ' not in master file')
		miss.write(acro9 + '\n')
	else:

		stalist.append(acro9)
		subprocess.call('bash /Calculs/utils/create_station_file_from_log.sh -name=' + acro9,shell=True)
miss.close()

out = Path(args.stalist).name.split('.')[0] + '.stafile'

stafile = open(out,'w')
stafile.write('\n')

problems = open('logs_issues.list','w')
for acro9 in stalist:
	try:
		stafnew = open(acro9 + '_STATION.sta')
		for line in stafnew:
			if line[0:1] != '#':
				stafile.write(line)
	
		stafile.write('\n')
		os.remove(acro9 + '_STATION.sta')

	except:
		print ('Unable to open ' + acro9 + '_STATION.sta')
		problems.write(acro9 + '\n')

problems.close()
stafile.close()

