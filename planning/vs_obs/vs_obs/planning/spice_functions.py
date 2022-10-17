# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:35:16 2022

@author: iant
"""

import numpy as np
import spiceypy as sp
from datetime import datetime


from vs_obs.config.constants import SPICE_METHOD, SPICE_TARGET, \
    SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER, \
    SPICE_SHAPE_MODEL_METHOD, SP_DPR, SPICE_DREF, SPICE_VENUS_RADIUS, \
    SPICE_FORMATSTR, SPICE_PRECISION, SPICE_DATETIME_FMT




def sza(et):
    """get solar zenith angle"""
    subpnt = sp.subpnt(SPICE_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER)
    subpnt_xyz = subpnt[0]
    #incidence angles
    surf_ilumin = sp.ilumin(SPICE_SHAPE_MODEL_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER, subpnt_xyz)
    incidence_angle = surf_ilumin[3] * sp.dpr()
    return incidence_angle



def sc_latlon(et):
    """get longitude and latitude of S/C"""
    #sub point data
    subpnt = sp.subpnt(SPICE_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER)
    subpnt_xyz = subpnt[0]
                
    #convert to lat/lons
    reclat = sp.reclat(subpnt_xyz)
    lon = reclat[1] * SP_DPR
    lat = reclat[2] * SP_DPR
    
    return lon, lat
    


def fov_lonlat(et, fov_vectors):
    """get lon/lat pairs for field of view vectors"""
    
    lonlats = np.zeros((len(fov_vectors), 2))
    for i, fov_vector in enumerate(fov_vectors):
        sincpt = sp.sincpt(SPICE_SHAPE_MODEL_METHOD, SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER, SPICE_DREF, fov_vector)[0]
        reclat = sp.reclat(sincpt)
        
        lonlats[i, :] = [reclat[1] * SP_DPR, reclat[2] * SP_DPR]
        
    return lonlats



def sc_alt(et):
    """get S/C orbit altitude"""
    # find obs position/velocity rel to mars in J2000
    obs2venus_spkezr = sp.spkezr(SPICE_TARGET, et, SPICE_PLANET_REFERENCE_FRAME, SPICE_ABCORR, SPICE_OBSERVER)
    
    #height of observer above Mars centre
    alt = sp.vnorm(obs2venus_spkezr[0][0:3]) - SPICE_VENUS_RADIUS

    return alt






def et2dt(et):
    """convert ephemeris time to datetime"""
    return datetime.strptime(sp.et2utc(et, SPICE_FORMATSTR, SPICE_PRECISION), SPICE_DATETIME_FMT)




def dt2et(dt):
    """convert datetime to ephemeris time"""
    return sp.utc2et(datetime.strftime(dt, SPICE_DATETIME_FMT))




def next_d2n_terminator(dt):
    """find first day-to-night terminator"""
    et = dt2et(dt)
    incidence_angle = sza(et)
    
    #if already on nightside
    if incidence_angle > 90:
        while incidence_angle > 90:
            et += 1. * 60.
            incidence_angle = sza(et)
            # print(incidence_angle)
    # print("**")
    while incidence_angle < 90:
        et += 1
        incidence_angle = sza(et)
    # print(incidence_angle)
    
    return et2dt(et)





def next_n2d_terminator(dt):
    """find first day-to-night terminator"""
    et = dt2et(dt)
    incidence_angle = sza(et)
    
    #if already on dayside
    if incidence_angle < 90:
        while incidence_angle < 90:
            et += 1. * 60.
            incidence_angle = sza(et)
            # print(incidence_angle)
    # print("**")
    while incidence_angle > 90:
        et += 1
        incidence_angle = sza(et)
    # print(incidence_angle)
    
    return et2dt(et)
