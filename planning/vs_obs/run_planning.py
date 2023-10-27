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
import time

from matplotlib.backends.backend_pdf import PdfPages

from vs_obs.other.load_spice_kernels import load_spice_kernels
load_spice_kernels()

from vs_obs.config.constants import CHANNEL_ID
from vs_obs.planning.make_orbit_plan import make_orbit_plan
from vs_obs.planning.make_coverage_grids import make_coverage_grids
from vs_obs.io.plot_footprint import plot_footprint
from vs_obs.io.plot_latitude_coverage import plot_latitude_coverage
# from vs_obs.planning.print_equator_altitude_footprint import print_equator_altitude_footprint




[channelShape, name, fov_centre, nvectors, fov_vectors] = sp.getfov(CHANNEL_ID, 4)


# for days in range(0, 300, 5):
#     utc_start_time = datetime(2035, 1, 1) + timedelta(days=days)

#     orbit_plan, h = make_orbit_plan(utc_start_time, fov_vectors)
#     print("### ", utc_start_time)
#     equator_altitude(orbit_plan)



# utc_start_time = datetime(2034, 11, 26) #mission start approx
utc_start_time = datetime(2035, 3, 18) #mission start approx
# utc_start_time = datetime(2035, 5, 16) #max altitude 451km
# utc_start_time = datetime(2035, 7, 20) #equatorial nightside max 395km
# utc_start_time = datetime(2035, 7, 25) #equatorial nightside min 288km
print("Ops start time", utc_start_time)

print("Making orbit plan")
start = time.time()
orbit_plan, h = make_orbit_plan(utc_start_time, fov_vectors)
end = time.time()
print("Orbit plan made in %0.1fs" %(end - start))        

print("Plotting footprints")
start = time.time()
plot_footprint(orbit_plan)
end = time.time()
print("Footprints plotted in %0.1fs" %(end - start))        


# daynight = "n"
# print("Making coverage grid")
# start = time.time()
# coverage_grids, lats_hr, lons_hr = make_coverage_grids(orbit_plan, daynight)
# end = time.time()
# print("Coverage grid made in %0.1fs" %(end - start))        

# fig = plot_latitude_coverage(coverage_grids, lats_hr, lons_hr, "Filters 1, 2, 3 and dark")
# fig = plot_latitude_coverage(coverage_grids, lats_hr, lons_hr, "Filter 1, 1, 1 and dark on orbit 2")
# fig = plot_latitude_coverage(coverage_grids, lats_hr, lons_hr, "Filter 1, 1, 1 and dark on orbits 2 and 4")



# """plot increasing latitudinal coverage with orbit numbers"""
# daynight = "n"
# with PdfPages("latitudinal_coverage.pdf") as pdf: #open pdf
#     for final_orbit in range(30, 3001, 30):
#         print("Making coverage maps orbits 0-%i" %final_orbit)
#         start = time.time()
#         coverage_grids, lats_hr, lons_hr = make_coverage_grids(orbit_plan, daynight, semiorbits=range(0, final_orbit, 1))
#         end = time.time()
#         print("Coverage map made in %0.1fs" %(end - start))        
#         fig = plot_latitude_coverage(coverage_grids, lats_hr, lons_hr, title="Orbits 0-%i" %final_orbit)
#         pdf.savefig()
#         plt.close()


#plot coverage on map for filter 1 (red)
# plot_footprint(orbit_plan)
# masked_grid = np.ma.masked_where(coverage_grids["1"] < 0.1, coverage_grids["1"])
# plt.imshow(masked_grid, extent=(-180, 180, -60, 60), cmap="Set1")



    
