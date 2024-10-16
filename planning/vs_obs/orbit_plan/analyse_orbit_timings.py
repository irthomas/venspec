# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 11:57:45 2024

@author: iant

ANALYSE ORBIT TIMINGS
"""

from scipy.spatial import ConvexHull, convex_hull_plot_2d
import os
import spiceypy as sp
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
# import time
# from matplotlib.backends.backend_pdf import PdfPages

# from planning.vs_obs.config.constants import CHANNEL_ID
# from planning.vs_obs.orbit_plan.make_orbit_plan import make_orbit_plan
# from planning.vs_obs.io.plot_footprint import plot_footprint

from planning.vs_obs.spice.find_lat_crossings import get_et_crossing_lat, get_et_crossing_lats
from planning.vs_obs.spice.find_sza_crossings import get_et_crossing_sza
# from tools.general.progress_bar import progress
from planning.vs_obs.spice.spice_functions import et2dt, get_centre_lonlat, get_corners_lonlats
from planning.vs_obs.instrument.filter_selection import filter_sequences
from planning.vs_obs.instrument.filter_wheel_functions import calculate_sequences_duration

from planning.vs_obs.config.constants import PRECOOLING_DURATION, STANDBY_DURATION, CCU_RECONFIGURATION_TIME


utc_start_time = datetime(2035, 6, 1)
# et_start = sp.utc2et(str(utc_start_time))
nightside_sza_limit = 95.0  # for all points in FOV
dayside_sza_limit = 95.0  # for centre of FOV


# get first nightside start time
et_nightside_start1 = get_et_crossing_sza(utc_start_time, nightside_sza_limit, "gt")
dt_nightside_start1 = et2dt(et_nightside_start1)
lon_nightside_start1, lat_nightside_start1 = get_centre_lonlat(et_nightside_start1)  # Latitude of FOV centre at the SZA crossing time

# start real observation 15 seconds earlier
et_nightside_obs_start1 = et_nightside_start1 - 15.0  # start 15 seconds early for dark
dt_nightside_obs_start1 = et2dt(et_nightside_obs_start1)
lon_nightside_obs_start1, lat_nightside_obs_start = get_centre_lonlat(et_nightside_obs_start1)  # Latitude of FOV centre at the SZA crossing time

# calculate the theoretical nightside start and end times for the 4 orbits, considering only sza/lat constraints

# once the real nightside starting latitude is known, we can calculate all starting times for the 4 orbits
# here the nightsides are on the descending latitude side
ascending_crossing_ets, descending_crossing_ets = get_et_crossing_lats(
    dt_nightside_obs_start1, dt_nightside_obs_start1 + timedelta(minutes=390), lat_nightside_obs_start)
