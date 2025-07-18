#!/usr/bin/python

import numpy as np
import sys


def resample_time_series(time_series_file,sampling="7d"):

    TS = np.loadtxt(time_series_file)

    if ".cats" in time_series_file:
        time=0
    elif ".spotgins" in time_series_file:
        time=1
    else:
        print("ERROR : time series file format is not recognized.")
        print("Only use .cats & .spotgins files !")
        exit()

    #############################################
    # Calculate the new time range in yyyy.yyyy 
    #############################################
    sampling = sampling.replace("d","*(1./365.25)")
    sampling = sampling.replace("w","*(7./365.25)")
    sampling = sampling.replace("m","*(1./12.)")
    sampling = sampling.replace("y","*1.")
    step=eval(sampling)

    #############################################
    # Make week-median series to detect ouliers
    #############################################
    TS_resampled=np.zeros((1,len(TS[0,:])))
    ranges=np.arange(TS[0,time]+step,TS[-1,time],step)
    ranges[-1]=TS[-1,time]+0.0001
    BEG=TS[0,time]
    for END in ranges:
        cond=(TS[:,time] >= BEG) & (TS[:,time] <= END)
        if np.any(cond):
            med = np.median(TS[cond,:],axis=0)
            TS_resampled = np.concatenate((TS_resampled,med[np.newaxis,:]))
        BEG=END

    return np.delete(TS_resampled,0,0)



if __name__ == "__main__":


    ##############################
    # Load cats time series file
    ##############################
    file=sys.argv[1]
    sampling=sys.argv[2]
    TS_resampled = resample_time_series(file,sampling=sampling)

    #################
    # Print results
    #################
    for i in range(len(TS_resampled[:,0])):
        print(" {0:8.0f} {1:12.7f}  {2:8.2f} {3:14.6f} {4:14.6f} {5:14.6f} {6:14.6f} {7:14.6f} {8:14.6f}".format(*TS_resampled[i,:]))
