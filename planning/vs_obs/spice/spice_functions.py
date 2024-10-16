# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:35:16 2022

@author: iant
"""

import numpy as np
import spiceypy as sp
from datetime import datetime

from planning.vs_obs.config.constants import SPICE_METHOD, SPICE_TARGET, \
    SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER, \
    SPICE_SHAPE_MODEL_METHOD, SP_DPR, SPICE_DREF, SPICE_VENUS_RADIUS, \
    SPICE_FORMATSTR, SPICE_PRECISION, SPICE_DATETIME_FMT


def get_fov_vectors(channel_id):
    """return vectors for FOV centre and corners in the VSH reference frame"""
    [channelShape, name, fov_centre, nvectors, fov_vectors] = sp.getfov(channel_id, 4)
    return fov_centre, fov_vectors


def get_centre_sza(et):
    """get solar zenith angle"""
    subpnt = sp.subpnt(SPICE_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER)
    subpnt_xyz = subpnt[0]
    # incidence angles
    surf_ilumin = sp.ilumin(SPICE_SHAPE_MODEL_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER, subpnt_xyz)
    incidence_angle = surf_ilumin[3] * sp.dpr()
    return incidence_angle


def get_sc_latlon(et):
    """get longitude and latitude of S/C"""
    # sub point data
    subpnt = sp.subpnt(SPICE_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER)
    subpnt_xyz = subpnt[0]

    # convert to lat/lons
    reclat = sp.reclat(subpnt_xyz)
    lon = reclat[1] * SP_DPR
    lat = reclat[2] * SP_DPR

    return lon, lat


def get_corner_lonlats(et, fov_vectors):
    """get lon/lat pairs for field of view vectors"""

    lonlats = np.zeros((len(fov_vectors), 2))
    for i, fov_vector in enumerate(fov_vectors):
        sincpt = sp.sincpt(SPICE_SHAPE_MODEL_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER, SPICE_DREF, fov_vector)[0]
        reclat = sp.reclat(sincpt)

        lonlats[i, :] = [reclat[1] * SP_DPR, reclat[2] * SP_DPR]

    return lonlats


def get_corners_lonlats(ets):
    """get lon/lat pairs for field of view vectors"""

    if isinstance(ets, list):
        ets = np.asarray(ets)

    n_ets = ets.shape[0]

    channel_id = sp.bods2c(SPICE_DREF)  # find channel id number
    [channel_shape, name, boresight_vector, n_corners, boresight_corners] = sp.getfov(channel_id, 4)

    fov_corners = [list(boresight_corner[:]) for boresight_corner in boresight_corners]

    surf_points = np.zeros((n_ets, n_corners, 3)) + np.nan  # use nan here, not -999
    surf_latlons = np.zeros((n_ets, n_corners, 2)) + np.nan  # use nan here, not -999

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

                reclat = sp.reclat(surf_points[j, i, :])

                surf_latlons[j, i, :] = [reclat[1] * SP_DPR, reclat[2] * SP_DPR]
    return surf_latlons


def get_centre_lonlat(et):
    """get lon/lat for field of view centre"""

    fov_vector = np.asarray([0., 0., 1.0])

    sincpt = sp.sincpt(SPICE_SHAPE_MODEL_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER, SPICE_DREF, fov_vector)[0]
    reclat = sp.reclat(sincpt)

    return [reclat[1] * SP_DPR, reclat[2] * SP_DPR]


def get_sc_alt(et):
    """get S/C orbit altitude"""
    # find obs position/velocity rel to mars in J2000
    obs2venus_spkezr = sp.spkezr(SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER)

    # height of observer above Mars centre
    alt = sp.vnorm(obs2venus_spkezr[0][0:3]) - SPICE_VENUS_RADIUS

    return alt


def et2dt(et):
    """convert ephemeris time to datetime"""
    return datetime.strptime(sp.et2utc(et, SPICE_FORMATSTR, SPICE_PRECISION), SPICE_DATETIME_FMT)


def dt2et(dt):
    """convert datetime to ephemeris time"""
    return sp.utc2et(datetime.strftime(dt, SPICE_DATETIME_FMT))
