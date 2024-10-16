# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 21:24:22 2024

@author: iant

ASSUMING THAT NIGHTSIDES WILL BE LIMITED TO A FEW DEGREES AWAY FROM TERMINATOR, CHECK LATITUDE COVERAGE NEAR POLES

"""


import spiceypy as sp
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
# import time

# from matplotlib.backends.backend_pdf import PdfPages

from tools.general.progress import progress
from tools.file.read_write_hdf5 import write_hdf5_from_dict, read_hdf5_to_dict

from tools.spice.load_spice_kernels import load_spice_kernels

from tools.spice.load_spice_kernels import SPICE_TARGET, \
    SPICE_PLANET_REFERENCE_FRAME, SPICE_OBSERVER, \
    SPICE_SHAPE_MODEL_METHOD, SP_DPR, SPICE_DREF, \
    SPICE_FORMATSTR, SPICE_PRECISION, SPICE_DATETIME_FMT, \
    SPICE_INTERCEPT_METHOD, SPICE_ABERRATION_CORRECTION, \
    SPICE_REFERENCE_FRAME, KILOMETRES_TO_AU, SPICE_PLANET_ID, SPICE_LONGITUDE_FORM


NA_VALUE = -999


# load spice kernels for nominal mission
load_spice_kernels("envision_study_et1_north_voi_ml008_v050.tm")


def make_coverage_h5(utc_start_time, ets):

    et_start = sp.utc2et(str(utc_start_time))

    # assume 1 second IT here, not used for plotting
    ets_all = np.asfarray([[et + et_start, et + et_start + 13] for et in ets]).T

    et_end = ets_all[0, -1]
    utc_end_time = sp.et2utc(et_end, SPICE_FORMATSTR, SPICE_PRECISION)
    utc_start_time = sp.et2utc(et_start, SPICE_FORMATSTR, SPICE_PRECISION)

    nSpectra = ets_all.shape[1]

    N_OBS_DATETIMES = 1  # ets_all.shape[0]

    CHANNEL_ID = sp.bods2c(SPICE_DREF)  # find channel id number
    # SPICE_VENUS_AXES = sp.bodvrd("VENUS", "RADII", 3)[1]  # no ellipsoid for Venus -> circular so only one radius value
    # SPICE_VENUS_RADIUS = SPICE_VENUS_AXES[0]

    dref = SPICE_DREF
    [channelShape, name, boresightVector, nvectors, boresightVectorboundsCorners] = sp.getfov(CHANNEL_ID, 4)

    # boresightVectorbounds = np.vstack((np.asarray([0., 0., 1.0]), boresightVectorboundsCorners))
    boresightVectorbounds = boresightVectorboundsCorners

    nPoints = boresightVectorbounds.shape[0]

    # make non-point dictionary
    # d = {}
    d_tmp = {}

    # make empty point dictionaries
    dp = {}
    dp_tmp = {}
    for point in range(nPoints):
        dp[point] = {}
        dp_tmp[point] = {}

    d_tmp["fov_corners"] = [list(boresightVectorbound[:]) for boresightVectorbound in boresightVectorbounds]

    print("Calculating surface point geometry")
    # Loop through array of times for a first time, storing surface intercept points
    # Valid=times when nadir pointed to Mars otherwise nan
    """make empty array nspectra x npoints x [start,end] x nsurface points"""
    surf_points = [np.zeros((nSpectra, nPoints, 3)) + np.nan for k in range(N_OBS_DATETIMES)]  # use nan here, not -999

    # surface intercepts
    # run in loop - need to catch off-Mars pointing errors
    for j in progress(range(nSpectra)):

        ets = [ets_all[k][j] for k in range(N_OBS_DATETIMES)]
        for i in range(nPoints):
            fov_corner = d_tmp["fov_corners"][i]
            try:
                sincpt = [sp.sincpt(SPICE_SHAPE_MODEL_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME,
                                    SPICE_ABERRATION_CORRECTION, SPICE_OBSERVER, dref, fov_corner) for et in ets]
            except sp.stypes.SpiceyError:  # error thrown if lat and lon point not on planet
                continue

            for k in range(N_OBS_DATETIMES):
                surf_points[k][j, i, :] = sincpt[k][0]

    nans_found = False
    all_nans = False

    n_nans = 0
    n_points = 0

    for i in range(nPoints):
        dp[i]["surf_xyz"] = [surf_points[k][:, i, :] for k in range(N_OBS_DATETIMES)]
        # dp[i]["point_xy"] = points[i] #save FOV points i.e. [0,0], [1,1], etc.

        surf_xyzs = np.asarray(dp[i]["surf_xyz"])

        # count the number of nans
        n_nans += np.count_nonzero(np.isnan(surf_xyzs))
        # count number of elements
        n_points += surf_xyzs.size

        # check if any nans found
        # nans = np.isnan(np.asarray(dp[i]["surf_xyz"]))
        if n_nans > 0:
            nans_found = True

        # check if all nans
        if n_nans == n_points:
            all_nans = True

    """code to check FOV pointing"""

    # check if FOV always on planet
    if not nans_found:  # if always on planet
        print("Nadir observation always points to planet")

    elif all_nans:  # if always off planet
        print("Warning: never points towards planet")

    else:  # if mixed (e.g. normal occs)
        print("Warning: only sometimes points towards planet (%i/%i points)", n_nans, n_points)

    print("Adding nadir point geometry")

    fovCornersAll = [a for a in d_tmp["fov_corners"]]

    # initialise empty arrays
    surfCoords = [np.zeros((nSpectra, nPoints, 3)) + NA_VALUE for k in range(N_OBS_DATETIMES)]
    surfRadius = [np.zeros((nSpectra, nPoints)) + NA_VALUE for k in range(N_OBS_DATETIMES)]
    surfLons = [np.zeros((nSpectra, nPoints)) + NA_VALUE for k in range(N_OBS_DATETIMES)]
    surfLats = [np.zeros((nSpectra, nPoints)) + NA_VALUE for k in range(N_OBS_DATETIMES)]
    surfLSTs = [np.zeros((nSpectra, nPoints)) + NA_VALUE for k in range(N_OBS_DATETIMES)]
    surfLSTHMSs = [np.zeros((nSpectra, nPoints, 3)) + NA_VALUE for k in range(N_OBS_DATETIMES)]
    surfIllumins = [np.zeros((nSpectra, nPoints, 3)) + NA_VALUE for k in range(N_OBS_DATETIMES)]
    surfSunSZAs = [np.zeros((nSpectra, nPoints)) + NA_VALUE for k in range(N_OBS_DATETIMES)]
    surfIncidenceAngles = [np.zeros((nSpectra, nPoints)) + NA_VALUE for k in range(N_OBS_DATETIMES)]
    surfEmissionAngles = [np.zeros((nSpectra, nPoints)) + NA_VALUE for k in range(N_OBS_DATETIMES)]
    surfPhaseAngles = [np.zeros((nSpectra, nPoints)) + NA_VALUE for k in range(N_OBS_DATETIMES)]

    # loop through spectra and points, adding to arrays
    for rowIndex in progress(range(nSpectra)):
        for k in range(N_OBS_DATETIMES):
            et = ets_all[k][rowIndex]

            # Points Datasets
            for pointIndex, fovCorner in enumerate(fovCornersAll):
                if not np.any(np.isnan(surf_points[k][rowIndex, pointIndex, :])):

                    surfCoords[k][rowIndex, pointIndex, :] = sp.reclat(surf_points[k][rowIndex, pointIndex, :])
                    surfLons[k][rowIndex, pointIndex] = surfCoords[k][rowIndex, pointIndex, 1] * SP_DPR
                    surfLats[k][rowIndex, pointIndex] = surfCoords[k][rowIndex, pointIndex, 2] * SP_DPR

                    surfRadius[k][rowIndex, pointIndex] = surfCoords[k][rowIndex, pointIndex, 0]

                    surfLSTHMSs[k][rowIndex, pointIndex, :] = sp.et2lst(
                        et, SPICE_PLANET_ID, surfCoords[k][rowIndex, pointIndex, 1], SPICE_LONGITUDE_FORM)[0: 3]

                    surfLSTs[k][rowIndex, pointIndex] = (surfLSTHMSs[k][rowIndex, pointIndex, 0] + surfLSTHMSs[k]
                                                         [rowIndex, pointIndex, 1]/60.0 + surfLSTHMSs[k][rowIndex, pointIndex, 2]/3600.0)

                    surfIllumins[k][rowIndex, pointIndex, :] = sp.ilumin(
                        SPICE_SHAPE_MODEL_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABERRATION_CORRECTION, SPICE_OBSERVER,
                        surf_points[k][rowIndex, pointIndex, :])[2: 5]

                    surfSunSZAs[k][rowIndex, pointIndex] = surfIllumins[k][rowIndex, pointIndex, 1]*SP_DPR
                    surfIncidenceAngles[k][rowIndex, pointIndex] = surfIllumins[k][rowIndex, pointIndex, 1]*SP_DPR
                    surfEmissionAngles[k][rowIndex, pointIndex] = surfIllumins[k][rowIndex, pointIndex, 2]*SP_DPR
                    surfPhaseAngles[k][rowIndex, pointIndex] = surfIllumins[k][rowIndex, pointIndex, 0]*SP_DPR

    """convert to arrays, add to dictionary"""
    for i in range(nPoints):
        dp[i]["surf_lon"] = np.asarray(surfLons)[:, :, i].T
        dp[i]["surf_lat"] = np.asarray(surfLats)[:, :, i].T
        dp[i]["surf_lst"] = np.asarray(surfLSTs)[:, :, i].T
        dp[i]["surf_sza"] = np.asarray(surfSunSZAs)[:, :, i].T
        dp[i]["surf_inc"] = np.asarray(surfIncidenceAngles)[:, :, i].T
        dp[i]["surf_emi"] = np.asarray(surfEmissionAngles)[:, :, i].T
        dp[i]["surf_pha"] = np.asarray(surfPhaseAngles)[:, :, i].T
        dp[i]["surf_dsk_dist"] = np.asarray(surfRadius)[:, :, i].T

    # just take acquisition start times
    lonlats = np.asarray([[dp[i]["surf_lon"][:, 0] for i in range(nPoints)], [dp[i]["surf_lat"][:, 0] for i in range(nPoints)]])
    szas = np.asarray([dp[i]["surf_sza"][:, 0] for i in range(nPoints)])

    h5_name = r"C:\Users\iant\Documents\DATA\sza_maps_utc_%s-%s" % (utc_start_time[0:11].replace(" ", "_"), utc_end_time[0:11].replace(" ", "_"))
    write_hdf5_from_dict(h5_name, {"et": ets_all[0, :], "lon": lonlats[0, :, :], "lat": lonlats[1, :, :], "sza": szas}, {}, {}, {})


# utc_start_time = datetime(2035, 3, 17, 1, 0, 0)
# ets = np.arange(0, 4*365*24*3600, 15)
# et_start = sp.utc2et(str(utc_start_time))
# et_end = max(ets) + et_start
# utc_end_time = sp.et2utc(et_end, SPICE_FORMATSTR, SPICE_PRECISION)
# utc_start_time = sp.et2utc(et_start, SPICE_FORMATSTR, SPICE_PRECISION)
# print(utc_start_time, utc_end_time)
# make_coverage_h5(utc_start_time, ets)


print("Loading data")
h5_name = r"C:/Users/iant/Documents/DATA/sza_maps_utc_2035_MAR_17-2039_MAR_16"
d = read_hdf5_to_dict(h5_name)[0]

utc_str = h5_name[-23:]

ets = d["et"]
lons = d["lon"]
lats = d["lat"]
szas = d["sza"]


NIGHTSIDE_SZA_LIMITS = [91., 95., 100.]

# plt.figure(figsize=(9, 9))
# for sza_limit in NIGHTSIDE_SZA_LIMITS:
#     print(sza_limit)

#     sza_max = np.max(szas, axis=0)  # max of the FOV corners for each acquisition
#     ixs_night_good = np.where(sza_max > sza_limit)[0]

#     lat_bins = np.arange(-90, 91, 1.0)

#     lons_night = lons[:, ixs_night_good]
#     lats_night = lats[:, ixs_night_good]
#     szas_night = szas[:, ixs_night_good]

#     freq, bins = np.histogram(np.ravel(lats[:, ixs_night_good]), bins=lat_bins)

#     plt.barh(bins[:-1], freq, height=1.0, alpha=0.5, label="SZA of whole FOV>%0.f degrees" % sza_limit)
#     plt.ylim((-90, 90))


# plt.title("Nightside latitude coverage histogram for different solar zenith angle limits\n2035 SEP 17 - 2037 SEP 15")
# plt.grid()
# plt.legend(loc="center left")
# plt.xlabel("Number of observation points at each latitude")
# plt.ylabel("Latitude (degrees)")

# plt.savefig("latitude_coverage_sza_histograms_%s.png" % utc_str)
# plt.close()

# for sza_limit in NIGHTSIDE_SZA_LIMITS:
#     print(sza_limit)

#     sza_max = np.max(szas, axis=0)  # max of the FOV corners for each acquisition
#     ixs_night_good = np.where(sza_max > sza_limit)[0]

#     lons_night = lons[:, ixs_night_good]
#     lats_night = lats[:, ixs_night_good]
#     szas_night = szas[:, ixs_night_good]

#     fig1, ax1 = plt.subplots(figsize=(14, 7), constrained_layout=True)
#     plot = ax1.scatter(lons_night, lats_night, c=szas_night, s=2, cmap="gnuplot_r", alpha=1.0)
#     ax1.set_xlabel("Longitude")
#     ax1.set_ylabel("Latitude")
#     ax1.set_xlim((-181, 181))
#     ax1.set_ylim((-90, 90))
#     cb1 = fig1.colorbar(plot)
#     cb1.set_label("SZA (degrees)", rotation=270, labelpad=10)
#     ax1.grid()
#     ax1.set_title("Solar zenith angle maps where entire VenSpec-H FOV>%.0f degrees\n%s" % (sza_limit, utc_str))

#     plt.savefig("sza_maps_utc_%s_sza_gt_%i.png" % (utc_str, sza_limit))
#     plt.close()


# take a slice through a lat/lon range and plot the szas
lat_range = [-1., 1.]


sza_limit = 90.

sza_max = np.max(szas, axis=0)  # max of the FOV corners for each acquisition
ixs_night_good = np.where((sza_max > sza_limit) & (lats[0, :] > lat_range[0]) & (lats[0, :] < lat_range[1]))[0]

lats_night = lats[:, ixs_night_good]
lons_night = lons[:, ixs_night_good]
szas_night = szas[:, ixs_night_good]
ets_night = ets[ixs_night_good]


plt.figure(figsize=(12, 6))
plt.title("VenSpec-H Nightside Observations")
plt.scatter((ets_night-ets_night[0])/(3600*24*365), lons_night[0, :], c=szas_night[0, :], cmap="gnuplot_r")
plt.xlabel("Years since start of mission")
plt.ylabel("Longitude of equator crossing (degrees)")
cb1 = plt.colorbar()
cb1.set_label("Solar zenith angle (degrees)", rotation=270, labelpad=10)
plt.grid()
plt.savefig("longitude_sza_nightside_equator_crossings.png", dpi=100)

# plt.figure(figsize=(9, 9))
# sza_bins = np.arange(90.0, 181.0, 1.0)

# freq, bins = np.histogram(np.ravel(szas[0, ixs_night_good]), bins=sza_bins)

# plt.bar(bins[:-1], freq, width=1.0)

# plt.title("Solar zenith angle distribution for latitudes between %0.0f and %0.0f" % (lat_range[0], lat_range[1]))
# plt.xlabel("Solar zenith angles")
# plt.ylabel("Number of equatorial observations at this solar zenith angle")
