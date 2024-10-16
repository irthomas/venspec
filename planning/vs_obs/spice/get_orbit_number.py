# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 21:34:07 2024

@author: iant

OLD METHOD - CALCULATE IT
NEW METHOD - GET FILE FROM ESA DATALABS
"""

import os
# import numpy as np
# import spiceypy as sp
from datetime import datetime
# import matplotlib.pyplot as plt
from planning.vs_obs.config.paths import paths
from planning.vs_obs.spice.spice_functions import et2dt
# from planning.vs_obs.spice.find_sza_crossings import get_corners_sza


def read_orbit_numbers():
    """new version: read orbit numbers from datalabs file, make dictionary for quick searching"""

    print("Reading in orbit number data")
    with open(os.path.join(paths["REFERENCE_DIRECTORY"], "quindec_cycle_orbit_table.txt"), "r") as f:
        lines = f.readlines()
    lines_split = [line.split() for line in lines]

    # make dictionary of type {(dt_start, dt_end): orbit_number}
    orbit_number_dts = {}
    for i in range(len(lines_split) - 1):
        dt_start = datetime.strptime(lines_split[i][0], "%Y-%m-%dT%H:%M:%S")
        dt_end = datetime.strptime(lines_split[i + 1][0], "%Y-%m-%dT%H:%M:%S")
        orbit_number = int(lines_split[i][1])
        orbit_number_dts[(dt_start, dt_end)] = orbit_number

    return orbit_number_dts


def get_orbit_number(orbit_number_dts, etdt):
    """find orbit number given an input et or dt"""

    # convert et to dt if et
    if isinstance(etdt, float):
        etdt = et2dt(etdt)

    # loop through until matching orbit time found
    for (dt_start, dt_end), orbit_number in orbit_number_dts.items():
        if etdt >= dt_start and etdt < dt_end:
            return orbit_number  # return as soon as found


# for testing
# etdt = 1117844370.0
# orbit_number = get_orbit_number(etdt)


# def get_local_minima(values):
#     """find indices of all local minima"""
#     indices = (np.diff(np.sign(np.diff(values))) > 1).nonzero()[0] + 1
#     return indices


# def get_local_maxima(values):
#     """find indices of all local minima"""
#     indices = (np.diff(np.sign(np.diff(values))) < -1).nonzero()[0] + 1
#     return indices


# def get_orbit_times(dt_before, dt_after, orbit_delimiter="north", resolution=5):
#     """make dictionary of orbit start and end times with orbit number"""

#     et_before = sp.utc2et(str(dt_before))
#     et_after = sp.utc2et(str(dt_after))

#     ets = np.arange(et_before, et_after, resolution)
#     lats = np.asarray([get_centre_lonlat(et)[1] for et in ets])

#     if orbit_delimiter == "north":
#         ixs = get_local_maxima(lats)
#     elif orbit_delimiter == "south":
#         ixs = get_local_minima(lats)

#     orbit_dts = {}
#     for i in range(len(ixs) - 1):
#         dt_start = et2dt(ets[ixs[i]])
#         dt_end = et2dt(ets[ixs[i+1]])
#         orbit_dts[(dt_start, dt_end)] = i

#     return orbit_dts

# orbit_dts = get_orbit_times(datetime(2035, 6, 1), datetime(2035, 6, 5))
