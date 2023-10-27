# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 12:43:13 2023

@author: iant


VENSPEC-H SZA AT THE EQUATOR
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
    
    for i in range(44000):
        
        if np.mod(i, 1000) == 0:
            print(i)
    
        if i == 0:
            #get science phase start/end ephemeris times
            dt_start = orbit_dict["science_start"] + timedelta(days=1)
            et_start = sp.utc2et(datetime.strftime(dt_start, SPICE_DATETIME_FMT))
    
            #find first equator crossing
            # 10 second resolution
            ets = np.arange(et_start, et_start + HALF_ORBIT_LENGTH, 10)
        else:
            ets = np.arange(et_start, et_start + 240.0, 10)
            
        
            lats = get_spice_info(ets, ["lat"])["lat"]
            
            #index at equator
            eq_ix = get_nearest_index(0.0, lats)
            
            #params at equator
            et_eq = ets[eq_ix]
            info_eq = get_spice_info(et_eq, ["lon", "sza"])
            lon_eq = info_eq["lon"]
            sza_eq = info_eq["sza"]
            
            szas["ets"].append(et_eq)
            szas["lons"].append(lon_eq)
            szas["szas"].append(sza_eq)
            
            #find next equator
            et_start = et_eq + HALF_ORBIT_LENGTH - 120.0
            
    return szas


if "szas" not in globals():
    szas = get_data()

print("Plot")
for key in szas:
    szas[key] = np.asarray(szas[key])

# plt.figure()
# plt.scatter(szas["ets"], szas["szas"], c=szas["lons"])
# plt.figure()
# plt.scatter(szas["lons"], szas["szas"], c=szas["ets"], s=1)

day_ixs = np.where(szas["szas"] < 90.0)[0]
night_ixs = np.where(szas["szas"] > 90.0)[0]

plt.figure(figsize=(12, 8))
plt.title("Longitude of equator crossing: dayside observations (colour = time)")
plt.scatter(szas["lons"][day_ixs], szas["szas"][day_ixs], c=szas["ets"][day_ixs], s=1)
plt.xlabel("Longitude")
plt.ylabel("Solar zenith angle")
plt.xlim((-180, 180))
plt.grid()
plt.xticks(np.arange(-180, 181, 30))
plt.savefig("dayside_sza_coverage.png")

plt.figure(figsize=(12, 8))
plt.title("Longitude of equator crossing: nightside observations (colour = time)")
plt.scatter(szas["lons"][night_ixs], szas["szas"][night_ixs], c=szas["ets"][night_ixs], s=1)
plt.xlabel("Longitude")
plt.ylabel("Solar zenith angle")
plt.xlim((-180, 180))
plt.grid()
plt.xticks(np.arange(-180, 181, 30))
plt.savefig("nightside_sza_coverage.png")


# #make strings
# dts = [sp.et2utc(et, "C", 0) for et in ets]

# n_points = len(ets)

