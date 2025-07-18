#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import regionmask
import sys

# #print the regional mask used here
# print(regionmask.defined_regions.ar6.all)

# https://regionmask.readthedocs.io/en/stable/defined_scientific.html

def get_indiv_mask(lon,lat):

    #create the mask
    lon=np.array([lon])
    lat=np.array([lat])
    mask = regionmask.defined_regions.ar6.all.mask(lon, lat)
    value = mask.values[0][0]

    #region que l'on souhaite créer :
    # 0 = **JOKER**
    # 1 = Europe (tout ce qui est à l’ouest de la Russie sauf Turquie).
    # 2 = Asie (inclus la Turquie, le proche orient, et tout ce qui est à l’Est de la Mer Rouge)
    # 3 = Afrique
    # 4 = Amérique du Nord (USA, Mexique, Canada et Alaska inclu)
    # 5 = Océanie (Australie, Nouvelle Zélande, Tonga, Samoa, Fidji, …) --> rugby power
    # 6 = Antarctique
    # 7 = Amérique du Sud (commence au Sud du Mexique)
    # 8 = Petits térritoires / Iles (inclassables ou pas évident que ça rentre dans les autres catégories, terres australes, antilles, ...)
    # 9 = Terres du Nord (Svalbard, Groenland, Islande)
    EUROPE=1
    ASIA=2
    AFRICA=3
    NORTH_AMERICA=4
    OCEANIA=5
    ANTARTICA=6
    SOUTH_AMERICA=7
    OTHERS=8
    NORTH_LANDS=9

    if (value==16) | (value==17) | (value==19):
        spot = EUROPE
    elif (value==18) | ((value>=28) & (value<=38)):
        spot = ASIA
    elif (value>=20) & (value<=27):
        spot = AFRICA
    elif (value>=1) & (value<=6):
        spot = NORTH_AMERICA
    elif (value>=39) & (value<=43):
        spot = OCEANIA
    elif (value>=44) & (value<=45):
        spot = ANTARTICA
    elif (value>=7) & (value<=15):
        spot = SOUTH_AMERICA
    elif (value==0) | (value==46):
        spot = NORTH_LANDS
    else:
        spot = OTHERS

    return spot




def create_mask():

    #create the mask
    lon = np.arange(-179.5, 180)
    lat = np.arange(-89.5, 90)
    mask = regionmask.defined_regions.ar6.all.mask(lon, lat)
    mask_arr = mask.values

    #region que l'on souhaite créer :
    # 0 = **JOKER**
    # 1 = Europe (tout ce qui est à l’ouest de la Russie sauf Turquie).
    # 2 = Asie (inclus la Turquie, le proche orient, et tout ce qui est à l’Est de la Mer Rouge)
    # 3 = Afrique
    # 4 = Amérique du Nord (USA, Mexique, Canada et Alaska inclu)
    # 5 = Océanie (Australie, Nouvelle Zélande, Tonga, Samoa, Fidji, …) --> rugby power
    # 6 = Antarctique
    # 7 = Amérique du Sud (commence au Sud du Mexique)
    # 8 = Petits térritoires / Iles (inclassables ou pas évident que ça rentre dans les autres catégories, terres australes, antilles, ...)
    # 9 = Terres du Nord (Svalbard, Groenland, Islande)
    EUROPE=1
    ASIA=2
    AFRICA=3
    NORTH_AMERICA=4
    OCEANIA=5
    ANTARTICA=6
    SOUTH_AMERICA=7
    OTHERS=8
    NORTH_LANDS=9

    mask_new = mask_arr.copy()
    mask_new[(mask_arr==16) | (mask_arr==17) | (mask_arr==19)] = EUROPE
    mask_new[(mask_arr==18) | ((mask_arr>=28) & (mask_arr<=38))] = ASIA
    mask_new[(mask_arr>=20) & (mask_arr<=27)] = AFRICA
    mask_new[(mask_arr>=1) & (mask_arr<=6)] = NORTH_AMERICA
    mask_new[(mask_arr>=39) & (mask_arr<=43)] = OCEANIA
    mask_new[(mask_arr>=44) & (mask_arr<=45)] = ANTARTICA
    mask_new[(mask_arr>=7) & (mask_arr<=15)] = SOUTH_AMERICA
    mask_new[(mask_arr==0) | (mask_arr==46)] = NORTH_LANDS
    mask_new[(mask_new>9)] = OTHERS

    return mask_new
    # #afficher la carte rapidement
    # plt.imshow(mask_new,origin='lower')
    # plt.colorbar()
    # plt.show()


if __name__ == "__main__":
    lon = float(sys.argv[1])
    lat = float(sys.argv[2])

    if lon > 180.:
        lon += -360.

    print(get_indiv_mask(lon,lat))

