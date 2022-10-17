# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 09:45:05 2022

@author: iant

VENSPEC-H OPERATIONS PLAN

HIGH DATA RATE:
    DAYSIDE:
        MINIMAL ACROSS-TRACK BINNING
        SHORT INTEGRATION TIME - MAY NEED MULTIPLE ACCUMULATIONS?
        DARK FRAMES EVERY N MEASUREMENTS
    NIGHTSIDE:
        BINNED ACROSS-TRACK
        LONG INTEGRATION TIME - SINGLE ACCUMULATION
        DARK FRAMES EVERY N MEASUREMENTS

"""

# import os
import spiceypy as sp
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

from vs_obs.other.load_spice_kernels import load_spice_kernels
load_spice_kernels()

from vs_obs.config.constants import CHANNEL_ID
from vs_obs.planning.make_orbit_plan import make_orbit_plan
from vs_obs.io.plot_footprint import plot_footprint



def equator_altitude(orbit_plan):
    orbit_ixs = []
    for orbit in orbit_plan:
        if "footprint_s" not in orbit[1].keys():
            continue
        
        if orbit[1]["semiorbit"] in orbit_ixs:
            continue
        
        
        footprint = orbit[1]["footprint_s"]
        lat = np.mean(footprint[:, 1])
        if lat > -1 and lat < 1:
            print(orbit[1]["daynight"], orbit[1]["semiorbit"], orbit[1]["alt"])

            footprint_s = orbit[1]["footprint_s"]
            footprint_e = orbit[1]["footprint_e"]
            
            #make footprint 
            footprint = np.zeros((5, 2))
            footprint[0, :] = footprint_s[1, :]
            footprint[1, :] = footprint_e[2, :]
            footprint[2, :] = footprint_e[3, :]
            footprint[3, :] = footprint_s[0, :]
            footprint[4, :] = footprint_s[1, :]

            plt.plot(footprint[:, 0], footprint[:, 1])
            print(footprint)
            

    
            orbit_ixs.append(orbit[1]["semiorbit"])


[channelShape, name, fov_centre, nvectors, fov_vectors] = sp.getfov(CHANNEL_ID, 4)


# for days in range(0, 300, 5):
#     utc_start_time = datetime(2035, 1, 1) + timedelta(days=days)

#     orbit_plan, h = make_orbit_plan(utc_start_time, fov_vectors)
#     print("### ", utc_start_time)
#     equator_altitude(orbit_plan)



# utc_start_time = datetime(2034, 11, 26) #mission start approx
# utc_start_time = datetime(2035, 5, 16) #max altitude 451km
# utc_start_time = datetime(2035, 7, 20) #equatorial nightside max 395km
utc_start_time = datetime(2035, 7, 25) #equatorial nightside min 288km
print(utc_start_time)

orbit_plan, h = make_orbit_plan(utc_start_time, fov_vectors)
equator_altitude(orbit_plan)


# plot_footprint(orbit_plan)


#footprint statistics
#make HR lon lat grid
#check if centre of each grid point is within each footprint
#record each filter separately
#select colour based on number of filters
#plot output for 1 orbit sequence, 14, 28 etc.











