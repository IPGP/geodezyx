#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#import matplotlib as mpl
#mpl.use('TkAgg')
from spytgins.files_classes import MasterStationFile
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import numpy as np

print (f"Working directory: {os.getcwd()}")
ma = 'metadata/stations/station_master_file.dat'

mast = MasterStationFile.from_file(ma)

acs = np.unique(mast.data['AC'].values)

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.Robinson())
ax.coastlines(color="black")
gl = ax.gridlines(draw_labels=True,linestyle="dashed",zorder=1)
ax.add_feature(cfeature.BORDERS,color='lightgrey',linewidth=0.5)

#taille des points
sp=10
lw = 1
colors = {"ALL":"black","ULR":"blue","OMP":"gold","EOST":"red","GRGS":"deepskyblue","ESGT":"chartreuse", "IPGP":"green","UGA":"pink"}

for ac in acs:
    print (ac)
    data = mast.data[mast.data['AC'] == ac][['LONGITUDE','LATITUDE']]

    zo= 2 ; ma = 'o' ; s = sp
    if ac == 'ALL':
        zo = 3 ; ma = '*'; s = 50

    netmap = ax.scatter(data['LONGITUDE'].values,
           data['LATITUDE'].values,
           c=colors[ac],
            marker=ma,
           s=s,edgecolor='black',
           transform=ccrs.PlateCarree(),zorder=zo,linewidth=0.3,alpha=0.8,label=ac + "(#" + str(len(data)) + ")")

plt.legend(bbox_to_anchor=(0.5, -0.18),fancybox=True, shadow=True,loc= 'lower center',ncol = 4,markerscale = 2)

plt.savefig('www/SPOTGINS_map.png',dpi = 300, bbox_inches='tight')

print ("Fichiers presents dans www")
for f in os.listdir("www"):
    print (f)
