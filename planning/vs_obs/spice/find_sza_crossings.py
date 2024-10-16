# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 11:22:19 2024

@author: iant

GET TIME OF FIRST NIGHTSIDE WHEN SZA OF FOV CORNERS MEETS A CERTAIN CONDITION

"""


# from planning.vs_obs.spice.find_lat_crossings import get_et_crossing_lat
import numpy as np
import spiceypy as sp
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
# from tools.general.progress_bar import progress
from planning.vs_obs.spice.spice_functions import et2dt, get_centre_lonlat
from planning.vs_obs.config.constants import SPICE_METHOD, SPICE_TARGET, \
    SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER, \
    SPICE_SHAPE_MODEL_METHOD, SP_DPR, SPICE_DREF, SPICE_VENUS_RADIUS, \
    SPICE_FORMATSTR, SPICE_PRECISION, SPICE_DATETIME_FMT


def get_corners_sza(ets):
    """get the solar zenith angles for VSH FOV corners at the given ephemeris times"""

    if isinstance(ets, list):
        ets = np.asarray(ets)

    n_ets = ets.shape[0]

    channel_id = sp.bods2c(SPICE_DREF)  # find channel id number
    [channel_shape, name, boresight_vector, n_corners, boresight_corners] = sp.getfov(channel_id, 4)

    fov_corners = [list(boresight_corner[:]) for boresight_corner in boresight_corners]

    surf_points = np.zeros((n_ets, n_corners, 3)) + np.nan  # use nan here, not -999
    surf_szas = np.zeros((n_ets, n_corners)) + np.nan  # use nan here, not -999

    # surface intercepts
    # run in loop - need to catch off-planet pointing errors
    for j, et in enumerate(ets):
        for i, fov_corner in enumerate(fov_corners):
            try:
                sincpt = sp.sincpt(SPICE_SHAPE_MODEL_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME,
                                   SPICE_ABCORR, SPICE_OBSERVER, SPICE_DREF, fov_corner)
            except sp.stypes.SpiceyError:  # error thrown if lat and lon point not on planet
                continue

            surf_points[j, i, :] = sincpt[0]

            if not np.any(np.isnan(surf_points[j, i, :])):

                surf_ilumin = sp.ilumin(
                    SPICE_SHAPE_MODEL_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER,
                    surf_points[j, i, :])[2:5]

                surf_szas[j, i] = surf_ilumin[1]*SP_DPR
    return surf_szas


def get_et_crossing_sza(dt_before, sza_limit, limit_type):
    """get the et when the maximum SZA of all corners crosses the given sza limit
    use gt for nightside start times, when SZA>given value"""

    dt_after = dt_before + timedelta(minutes=5)

    et_before = sp.utc2et(str(dt_before))
    et_after = sp.utc2et(str(dt_after))

    ets = np.arange(et_before, et_after, 1)
    surf_szas = get_corners_sza(ets)
    max_szas = np.max(surf_szas, axis=1)

    if limit_type == "gt":
        if max_szas[0] > sza_limit:
            print("Error: given time was already beyond the SZA limit")
    elif limit_type == "lt":
        if max_szas[0] < sza_limit:
            print("Error: given time was already beyond the SZA limit")
    else:
        print("Error: limit type %s must be lt or gt" % limit_type)

    extend_range = False
    if limit_type == "gt" and max_szas[-1] < sza_limit:
        extend_range = True
    if limit_type == "lt" and max_szas[-1] > sza_limit:
        extend_range = True

    if extend_range:
        # if insufficient points, extend to whole half orbit
        print("Warning: transition not found, extending range")

        dt_after = dt_before + timedelta(minutes=40)
        et_before = sp.utc2et(str(dt_before))
        et_after = sp.utc2et(str(dt_after))

        ets = np.arange(et_before, et_after, 1)
        surf_szas = get_corners_sza(ets)
        max_szas = np.max(surf_szas, axis=1)

    if limit_type == "gt":
        et_transition = np.interp(sza_limit, max_szas, ets)
    elif limit_type == "lt":  # reverse arrays
        et_transition = np.interp(sza_limit, max_szas[::-1], ets[::-1])

    # print(list(zip(max_szas, ets-ets[0])))
    # print(et_transition-ets[0])

    print("Corner SZAs at transition:", get_corners_sza([et_transition]), "crossing point", get_centre_lonlat(et_transition))
    return et_transition


def get_local_minima_maxima_or_equals(values):
    """find indices of all local minima and maxima, including those with equal adjacent values """

    indices_min = (np.diff(np.sign(np.diff(values))) > 0).nonzero()[0] + 1
    indices_max = (np.diff(np.sign(np.diff(values))) < 0).nonzero()[0] + 1
    return np.sort(np.concatenate((indices_min, indices_max)))


def get_et_crossing_szas(dt_before, dt_after, sza, limit_type="all", plot=True):
    """get all SZAs FOV crossings within a given time range for both ascending and descending orbits
    limit type refers to whether any FOV corner value should be greater than the limit or all should be"""

    # # for testing
    # dt_before = datetime(2035, 6, 1)
    # dt_after = datetime(2035, 6, 1, 6, 31, 45, 856000)
    # sza = 95.0
    # limit_type = "all"

    et_before = sp.utc2et(str(dt_before))
    et_after = sp.utc2et(str(dt_after))

    ets = np.arange(et_before, et_after, 1)
    surf_szas = get_corners_sza(ets)

    if limit_type == "any":
        szas = np.max(surf_szas, axis=1)
    elif limit_type == "all":
        szas = np.min(surf_szas, axis=1)
    elif limit_type == "mean":
        szas = np.mean(surf_szas, axis=1)
    else:
        print("Error: limit type incorrectly specified")

    min_max_ixs = np.concatenate(([0], get_local_minima_maxima_or_equals(szas)))

    ets_split = [ets[min_max_ixs[i]:min_max_ixs[i+1]] for i in range(len(min_max_ixs)-1)]
    szas_split = [szas[min_max_ixs[i]:min_max_ixs[i+1]] for i in range(len(min_max_ixs)-1)]

    ascending_crossing_ets = []
    descending_crossing_ets = []
    for ets_orbit, szas_orbit in zip(ets_split, szas_split):
        if szas_orbit[5] > szas_orbit[0]:  # szas ascending
            et_crossing = np.interp(sza, szas_orbit, ets_orbit)
            ascending_crossing_ets.append(et_crossing)
        if szas_orbit[5] < szas_orbit[0]:  # szas descending
            et_crossing = np.interp(sza, szas_orbit[::-1], ets_orbit[::-1])
            descending_crossing_ets.append(et_crossing)

    if plot:
        plt.figure()
        for i in range(4):
            plt.plot(ets, surf_szas[:, i])
        plt.axhline(y=sza, c="k", linestyle="--")
        for ix in min_max_ixs:
            plt.axvline(x=ets[ix])

        for et in ascending_crossing_ets:
            plt.scatter(et, sza, s=20, c="k")
        for et in descending_crossing_ets:
            plt.scatter(et, sza, s=20, c="r")

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

    if ascending_crossing_ets[0] < descending_crossing_ets[0]:
        crossings["type"] = "Nightside first"
    else:
        crossings["type"] = "Dayside first"

    return crossings
