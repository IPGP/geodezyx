#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from spytgins.spotgins_utils import functions
import sys
import re
import datetime
import numpy as np
import matplotlib.pyplot as plt
import os

plt.style.use('perso')


sol_file=sys.argv[1]
name=sol_file.split("/")[-1][:9]

station_infos = os.environ["station_infos"]
jjj = np.array([datetime.datetime.strptime(s, " %Y%m%d").timetuple().tm_yday for s in re.split("-o|\n", station_infos)[1::2]])
yyyy= np.array([datetime.datetime.strptime(s, " %Y%m%d").timetuple().tm_year for s in re.split("-o|\n", station_infos)[1::2]])
last= np.array([datetime.datetime.strptime(s[0:5]+"1231", " %Y%m%d").timetuple().tm_yday for s in re.split("-o|\n", station_infos)[1::2]])
dates_changes=yyyy + jjj/last

offset_file="/".join(sol_file.split("/")[:-1])+"/../"+name+"_OFFSET.sort"
if os.path.exists(offset_file):
    offsets = np.loadtxt("/".join(sol_file.split("/")[:-1])+"/../"+name+"_OFFSET.sort",usecols=[0])
else:
    offsets = np.array([])

TS = np.loadtxt(sol_file)

fig = plt.figure(figsize = (14,8))
ax1 = plt.subplot(311)
cursor1 = functions.SnaptoCursor(ax1, TS[:,1], TS[:,3]-np.median(TS[:,3]))
cid1 = plt.connect('motion_notify_event', cursor1.mouse_move)
ax1.plot(TS[:,1],TS[:,3]-np.median(TS[:,3]),linestyle='-',color='tab:blue')
for off in dates_changes:
    ax1.axvline(x=off,linestyle='--',color='k')
for off2 in offsets:
    ax1.axvline(x=off2,linestyle='--',color='r')
ax1.set_title("East", loc='right',fontsize=23, va='center')
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel("mm")
ax1.set_xlim(2000,2023)
plt.title(name+" -- SPOTGINS")

ax2 = plt.subplot(312, sharex=ax1)
cursor2 = functions.SnaptoCursor(ax2, TS[:,1], TS[:,4]-np.median(TS[:,4]))
cid2 = plt.connect('motion_notify_event', cursor2.mouse_move)
ax2.plot(TS[:,1],TS[:,4]-np.median(TS[:,4]),linestyle='-',color='tab:orange')
for off in dates_changes:
    ax2.axvline(x=off,linestyle='--',color='k')
for off2 in offsets:
    ax2.axvline(x=off2,linestyle='--',color='r')
ax2.set_title("Nort", loc='right',fontsize=23, va='center')
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("mm")
ax2.set_xlim(2000,2023)

ax3 = plt.subplot(313, sharex=ax1)
cursor3 = functions.SnaptoCursor(ax3, TS[:,1], TS[:,5]-np.median(TS[:,5]))
cid3 = plt.connect('motion_notify_event', cursor3.mouse_move)
ax3.plot(TS[:,1],TS[:,5]-np.median(TS[:,5]),linestyle='-',color='tab:green')
for off in dates_changes:
    ax3.axvline(x=off,linestyle='--',color='k')
for off2 in offsets:
    ax3.axvline(x=off2,linestyle='--',color='r')
ax3.set_title("Vert", loc='right',fontsize=23, va='center')
ax3.set_xlabel("years")
ax3.set_ylabel("mm")
ax3.set_xlim(2000,2023)

plt.show()

