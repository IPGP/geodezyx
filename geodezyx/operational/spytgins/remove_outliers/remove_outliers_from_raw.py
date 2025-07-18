#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from spytgins.spotgins_utils import functions
import sys


if __name__ == "__main__":

    ##############################
    # Load cats time series file
    ##############################
    file=sys.argv[1]
    tol=sys.argv[2]
    TS_final = functions.remove_outliers_raw(file,tol=float(tol))

    #################
    # Print results
    #################
    for i in range(len(TS_final[:,0])):
        print(" {0:8.0f} {1:12.7f}  {2:8.2f} {3:14.6f} {4:14.6f} {5:14.6f} {6:14.6f} {7:14.6f} {8:14.6f}".format(*TS_final[i,:]))



