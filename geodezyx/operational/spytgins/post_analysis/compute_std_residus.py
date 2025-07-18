#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from spytgins.spotgins_utils.functions import *
import sys
import re


if __name__ == "__main__":
    time_series_file=sys.argv[1]
    fit=sys.argv[2]

    points = os.environ["points"]
    li = [float(s) for s in points.split()]
    POINTS = np.reshape(np.array(li),(int(len(li)/4),4))

    compare_std2sol(time_series_file,fit,POINTS)

    # if re.match(r"^([0-9]+\.?[0-9]*[e,E]?[\+,-]?[0-9]*)(,[0-9]+\.?[0-9]*[e,E]?[\+,-]?[0-9]*)*$",dates):
    #     list_dates = np.array([float(s) for s in dates.split(",")])
    # else:
    #     exit()
    # if re.match(r"^([0-9]+\.?[0-9]*[e,E]?[\+,-]?[0-9]*)(,[0-9]+\.?[0-9]*[e,E]?[\+,-]?[0-9]*)*$",east):
    #     list_east = np.array([float(s) for s in east.split(",")])
    # else:
    #     exit()
    # if re.match(r"^([0-9]+\.?[0-9]*[e,E]?[\+,-]?[0-9]*)(,[0-9]+\.?[0-9]*[e,E]?[\+,-]?[0-9]*)*$",nort):
    #     list_nort = np.array([float(s) for s in nort.split(",")])
    # else:
    #     exit()
    # if re.match(r"^([0-9]+\.?[0-9]*[e,E]?[\+,-]?[0-9]*)(,[0-9]+\.?[0-9]*[e,E]?[\+,-]?[0-9]*)*$",vert):
    #     list_vert = np.array([float(s) for s in vert.split(",")])
    # else:
    #     exit()
    # points=np.column_stack((list_dates,list_east,list_nort,list_vert))


