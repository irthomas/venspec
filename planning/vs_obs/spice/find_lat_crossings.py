# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 12:25:30 2024

@author: iant
"""

import numpy as np
import spiceypy as sp
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
# from tools.general.progress_bar import progress
from planning.vs_obs.spice.spice_functions import et2dt, get_centre_lonlat
from planning.vs_obs.spice.find_sza_crossings import get_corners_sza


def get_et_crossing_lat(dt_before, latitude, direction):
    """get the et when the centre of the FOV crosses a given latitude"""

    dt_after = dt_before + timedelta(minutes=5)

    et_before = sp.utc2et(str(dt_before))
    et_after = sp.utc2et(str(dt_after))

    ets = np.arange(et_before, et_after, 1)
    lats = np.asarray([get_centre_lonlat(et)[1] for et in ets])

    # # try the interpolation, use polyfit as could be multiple crossings
    # x = np.arange(len(lats)) - len(lats) / 2.0
    # polyfit = np.polyfit(x, lats, 2)

    # check all increasing or decreasing
    diff = np.diff(lats)
    if not np.all(diff > 0.0) or not np.all(diff < 0.0):
        # if partially increasing or decreasing, assume transition is at the start
        if direction == "a":  # keep ascending values
            ixs = np.where(diff > 0.0)[0]
        if direction == "d":  # keep descending values and reverse indices
            ixs = np.where(diff < 0.0)[0][::-1]

        if len(ixs) < 3:
            print("Error: insufficient monotonically increasing or decreasing points (%i points)" % len(ixs), lats)
            return 0.0

    if latitude < min(lats) or latitude > max(lats):
        print("Error: latitude %0.2f not in given range" % latitude, lats)
        return 0.0

    et_crossing = np.interp(latitude, lats[ixs], ets[ixs])

    print("Latitude of FOV centre at crossing point", get_centre_lonlat(et_crossing))
    return et_crossing


def get_local_minima_maxima_or_equals(values):
    """find indices of all local minima and maxima, including those with equal adjacent values """

    indices_min = (np.diff(np.sign(np.diff(values))) > 0).nonzero()[0] + 1
    indices_max = (np.diff(np.sign(np.diff(values))) < 0).nonzero()[0] + 1
    return np.sort(np.concatenate((indices_min, indices_max)))


def get_et_crossing_lats(dt_before, dt_after, latitudes, plot=True):
    """get all latitude crossings within a given time range for both ascending and descending """

    # convert to et and get latitudes for every second
    et_before = sp.utc2et(str(dt_before))
    et_after = sp.utc2et(str(dt_after))

    ets = np.arange(et_before, et_after, 1)
    lats = np.asarray([get_centre_lonlat(et)[1] for et in ets])

    # find local minima and maxima including last point
    min_max_ixs = list(get_local_minima_maxima_or_equals(lats)) + [len(lats)]

    ets_split = [ets[min_max_ixs[i]:min_max_ixs[i+1]] for i in range(len(min_max_ixs)-1)]
    lats_split = [lats[min_max_ixs[i]:min_max_ixs[i+1]] for i in range(len(min_max_ixs)-1)]

    ascending_crossing_ets = []
    descending_crossing_ets = []
    for latitude in latitudes:
        for ets_orbit, lats_orbit in zip(ets_split, lats_split):
            if lats_orbit[5] > lats_orbit[0]:  # lats ascending
                et_crossing = np.interp(latitude, lats_orbit, ets_orbit)
                # only save if not equal to last point in et range (i.e. interpolation error)
                if et_crossing != ets_orbit[-1]:
                    ascending_crossing_ets.append(et_crossing)
            else:  # lats descending
                et_crossing = np.interp(latitude, lats_orbit[::-1], ets_orbit[::-1])
                # only save if not equal to last point in et range (i.e. interpolation error)
                if et_crossing != ets_orbit[-1]:
                    descending_crossing_ets.append(et_crossing)

    ascending_crossing_ets = sorted(ascending_crossing_ets)
    descending_crossing_ets = sorted(descending_crossing_ets)

    crossings = {"ascending": {}, "descending": {}}

    for crossing_et in ascending_crossing_ets:
        corner_szas = get_corners_sza([crossing_et])
        dt = et2dt(crossing_et)
        lon, lat = get_centre_lonlat(crossing_et)
        crossings["ascending"][crossing_et] = {"dt": dt, "corner szas": corner_szas, "centre lon": lon, "centre lat": lat}

    for crossing_et in descending_crossing_ets:
        corner_szas = get_corners_sza([crossing_et])
        dt = et2dt(crossing_et)
        lon, lat = get_centre_lonlat(crossing_et)
        crossings["descending"][crossing_et] = {"dt": dt, "corner szas": corner_szas, "centre lon": lon, "centre lat": lat}

    if plot:
        plt.figure()
        plt.plot(ets, lats, alpha=0.5)

        for latitude in latitudes:
            plt.axhline(y=latitude, c="k", linestyle="--")
        for et in crossings["ascending"].keys():
            plt.scatter(et, crossings["ascending"][et]["centre lat"], s=20, c="k")
        for et in crossings["descending"].keys():
            plt.scatter(et, crossings["descending"][et]["centre lat"], s=20, c="r")

    # all_ets = sorted(ascending_crossing_ets + descending_crossing_ets)
    if descending_crossing_ets[0] < ascending_crossing_ets[0]:
        crossings["type"] = "Latitudes descending"

        if descending_crossing_ets[1] > ascending_crossing_ets[0]:
            print("Error: incomplete first orbit")
    else:
        crossings["type"] = "Latitudes ascending"

        if ascending_crossing_ets[1] > descending_crossing_ets[0]:
            print("Error: incomplete first orbit")

    return crossings
