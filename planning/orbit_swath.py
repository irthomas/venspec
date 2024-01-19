# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 09:29:24 2023

@author: iant

VENSPEC-H LONGITUDE SWATH OVERLAP
"""


import os
import numpy as np
import spiceypy as sp
from datetime import datetime, timedelta

import matplotlib.pyplot as plt

from tools.spice.load_spice_kernels import load_spice_kernels, SPICE_DATETIME_FMT, \
    SPICE_METHOD, SPICE_ABCORR, SPICE_SHAPE_MODEL_METHOD

from tools.general.progress import progress
from tools.general.get_nearest_index import get_nearest_index

#load spice kernels
orbit_name = "EnVision_ET1_2031_NorthVOI"
orbit_dict = load_spice_kernels("%s.bsp" %orbit_name)

SPICE_OBSERVER = "-668"
HALF_ORBIT_LENGTH = 47*60.0
N_ORBITS = 3000


# START_DATETIME = datetime(2035, 3, 18, 0, 0)
START_DATETIME = datetime(2038, 4, 18, 0, 0)

VENUS_RADIUS = 6051.8 #km


def get_spice_info(ets, out):

    d = {}
    
    #check if single value => convert to 1 element array
    if isinstance(ets, np.generic):
        ets = [ets]
        
    

    #sub-observer point data
    subpnts = [sp.subpnt(SPICE_METHOD, "VENUS", et, "IAU_VENUS", SPICE_ABCORR, SPICE_OBSERVER) for et in ets]
    
    if any(x in out for x in ["lat", "lon", "sza"]):
        subpnts_xyz = [subpnt[0] for subpnt in subpnts]
    if any(x in out for x in ["lat", "lon"]):
        reclats = [sp.reclat(subpnt_xyz) for subpnt_xyz in subpnts_xyz]
    if "lon" in out:
        lons_rad = [reclat[1] for reclat in reclats]
        lons = np.asfarray(lons_rad) * sp.dpr()
        d["lon"] = lons
    if "lat" in out:
        lats_rad = [reclat[2] for reclat in reclats]
        lats = np.asfarray(lats_rad) * sp.dpr()
        d["lat"] = lats
    if "sza" in out:
        #incidence angles
        surf_ilumin = [sp.ilumin(SPICE_SHAPE_MODEL_METHOD, "VENUS", et, "IAU_VENUS", SPICE_ABCORR, SPICE_OBSERVER, subpnt_xyz) for et, subpnt_xyz in zip(ets, subpnts_xyz)]
        #trgepc, srfvec, sslphs, sslsol, sslemi
        szas = [ilumin[3] * sp.dpr() for ilumin in surf_ilumin]
        
        d["sza"] = szas
    
    return d



    


def get_data():
    print("SPICE calculations")
    szas = {"ets":[], "lons":[], "szas":[]}
    
    for i in range(N_ORBITS):
        
        if np.mod(i, int(N_ORBITS/20)) == 0:
            print(i)
    
        if i == 0:
            #get science phase start/end ephemeris times
            dt_start = START_DATETIME + timedelta(days=1)
            et_start = sp.utc2et(datetime.strftime(dt_start, SPICE_DATETIME_FMT))
    
            #find first equator crossing
            # 10 second resolution
            ets = np.arange(et_start, et_start + HALF_ORBIT_LENGTH, 1)
        else:
            ets = np.arange(et_start, et_start + 240.0, 0.1)
            
        
        lats = get_spice_info(ets, ["lat"])["lat"]
        
        #index at equator
        eq_ix = get_nearest_index(0.0, lats)
        
        #params at equator
        et_eq = ets[eq_ix]
        # print(i, et_eq - et_start)
        info_eq = get_spice_info(et_eq, ["lon", "sza"])
        lon_eq = info_eq["lon"]
        sza_eq = info_eq["sza"]
        
        szas["ets"].append(et_eq)
        szas["lons"].append(lon_eq)
        szas["szas"].append(sza_eq)
        
        #find next equator
        et_start = et_eq + HALF_ORBIT_LENGTH - 120.0
            
    return szas


if "geom_d" not in globals():
    geom_d = get_data()

print("Plot")
for key in geom_d:
    geom_d[key] = np.asarray(geom_d[key])

# plt.figure()
# plt.scatter(szas["ets"], szas["szas"], c=szas["lons"])
# plt.figure()
# plt.scatter(szas["lons"], szas["szas"], c=szas["ets"], s=1)

day_ixs = np.where(geom_d["szas"] < 90.0)[0][5:]
night_ixs = np.where(geom_d["szas"] > 90.0)[0][5:]

eq_crossing_lons = np.squeeze(geom_d["lons"][night_ixs])
# eq_crossing_szas = geom_d["szas"][night_ixs]

eq_crossing_days = (np.squeeze(geom_d["ets"][night_ixs]) - geom_d["ets"][night_ixs[0]])/3600.0 / 24.0

lons_diff = np.diff(eq_crossing_lons)

# plt.scatter(eq_crossing_lons, np.zeros_like(eq_crossing_lons), c=eq_crossing_szas)

km_diff = (2 * np.pi * VENUS_RADIUS) * lons_diff / 360.0

#remove transitions
trans_ixs = np.where(np.abs(km_diff) > 100.0)[0]
for trans_ix in trans_ixs:
    km_diff[trans_ix] = km_diff[trans_ix-1]


plt.figure(figsize=(12, 8))
plt.title("Longitudinal difference between consecutive orbits when crossing the equator")
plt.scatter(eq_crossing_days[:-1], km_diff)
plt.xlabel("Days elapsed starting %s" %str(START_DATETIME))
plt.ylabel("Distance between consecutive orbits (km)")
plt.grid()
plt.savefig("distance_between_consecutive_orbits_%s.png" %(str(START_DATETIME)[0:10]))


