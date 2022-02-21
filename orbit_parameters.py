# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 10:13:18 2022

@author: iant

CHECK ENVISION KERNELS FOR
1. EARTH/VENUS LIMB SCAN OPPORTUNITIES
2. ORBIT ALTITUDES
3. LOCAL SOLAR TIME

"""

import os
import spiceypy as sp
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

from matplotlib.backends.backend_pdf import PdfPages


from tools.file.paths import paths
# from tools.plotting.get_colours import get_colours

# SPICE_FORMATSTR = "C"
# SPICE_PRECISION = 0
SPICE_DATETIME_FMT = "%Y %b %d %H:%M:%S.%f"
SHORT_DATETIME_FMT = "%d %b %Y"

SPICE_METHOD = "INTERCEPT/ELLIPSOID"
SPICE_ABCORR = "NONE"
# SPICE_OBSRVR = "-668"

SPICE_LONGITUDE_FORM = "PLANETOCENTRIC"
SPICE_PLANET_ID = 299 # venus

kernel_root_dir = paths["KERNEL_ROOT_DIRECTORY"]

#load each kernel manually

kernels = [
    "fk/envision_v00.tf",
    "fk/envision_sc_venus_npo_v00.tf",
    "lsk/naif0012.tls",
    "pck/pck00010.tpc",
    "pck/de-403-masses.tpc",
    "spk/de432s.bsp",
    ]

orbit_d ={
    "EnVision_ALT_T1_2032_NorthVOI.bsp":{
        "title":"ALT_T1_2032_NorthVOI (2022/01/18)\nAlternative backup scenario with Venus Orbit Insertion over the Northern hemisphere",
        "science_start":datetime(2034, 12, 20),
        "science_end":datetime(2038, 12, 18),
        "spice_observer":"-999", #error in kernel production => should be -668
        },
    "EnVision_ALT_T4_2032_SouthVOI.bsp":{
        "title":"ALT_T4_2032_SouthVOI (2022/01/18)\nAlternative baseline scenario with Venus Orbit Insertion over the Southern hemisphere",
        "science_start":datetime(2035, 3, 27),
        "science_end":datetime(2039, 3, 25),
        "spice_observer":"-999", #error in kernel production => should be -668
        },
    "EnVision_ESC_T4_2032_NorthVOI.bsp":{
        "title":"ESC_T4_2032_NorthVOI (2020/06/10)\nLaunch into direct escape at -5.0 deg declination and Venus Orbit Insertion over the Northern hemisphere.",
        "science_start":datetime(2035, 6, 15),
        # "science_end":datetime(2035, 9, 14),
        "science_end":datetime(2039, 6, 14),
        "spice_observer":"-668",
        },
    "EnVision_ESC_T2_2032_SouthVOI.bsp":{
        "title":"ESC_T2_2032_SouthVOI (2019/07/26)\nLaunch into direct escape at -3.0 deg declination and Venus Orbit Insertion over the Southern hemisphere",
        "science_start":datetime(2034, 11, 26),
        "science_end":datetime(2040, 3, 24),
        "spice_observer":"-668",
        },
    "EnVision_HEO_T2_2032_NorthVOI.bsp":{
        "title":"HEO_T2_2032_NorthVOI (2019/07/26)\nLaunch into HEO (310x100,000 km) and Venus Orbit Insertion over the Northern hemisphere",
        "science_start":datetime(2034, 11, 26),
        "science_end":datetime(2040, 3, 24),
        "spice_observer":"-668",
        },
    }



for key in orbit_d.keys():
    for kernel in kernels:
        sp.furnsh(os.path.join(kernel_root_dir, os.path.normcase(kernel)))
        
    sp.furnsh(os.path.join(kernel_root_dir, "spk", key))
    
    min_alt = 300
    max_alt = 300
        
    
    with PdfPages("%s.pdf" %key) as pdf:
    
        #plot for each MTP i.e. 28 days at a time
        delta_days = 28
        
        #get start datetime of each MTP for whole of mission
        dt_starts = [orbit_d[key]["science_start"] + timedelta(days=delta_days*i) for i in np.arange(np.floor((orbit_d[key]["science_end"] - orbit_d[key]["science_start"]).days/delta_days))]
    
        #loop through MTPs
        for mtp, dt_start in  enumerate(dt_starts):
        
            #get science phase start/end ephemeris times
            dt_end = dt_start + timedelta(days=delta_days)
            et_start = sp.utc2et(datetime.strftime(dt_start, SPICE_DATETIME_FMT))
            et_end = sp.utc2et(datetime.strftime(dt_end, SPICE_DATETIME_FMT))
            
            #make list of ephemeris times        
            ets = np.arange(et_start, et_end, 60.0) #one value per 60 seconds
            
            
            #spice calculations
            venus_radius = sp.bodvrd("VENUS", "RADII", 3)[1][0] #no ellipsoid for Venus -> circular so only one radius value
    
            #sub point data
            subpnts = [sp.subpnt(SPICE_METHOD, "VENUS", et, "IAU_VENUS", SPICE_ABCORR, orbit_d[key]["spice_observer"]) for et in ets]
            subpnts_xyz = [subpnt[0] for subpnt in subpnts]
            
            
            #convert to lat/lons
            reclats = [sp.reclat(subpnt_xyz) for subpnt_xyz in subpnts_xyz]
            lons_rad = [reclat[1] for reclat in reclats]
            lats_rad = [reclat[2] for reclat in reclats]
            lons = np.asfarray(lons_rad) * sp.dpr()
            lats = np.asfarray(lats_rad) * sp.dpr()
            
            #get orbit altitude
            # find obs position/velocity rel to mars in J2000
            obs2venus_spkezrs = [sp.spkezr("VENUS", et, "IAU_VENUS", SPICE_ABCORR, orbit_d[key]["spice_observer"]) for et in ets]
            
            #height of observer above Mars centre
            alts = [sp.vnorm(spkezr[0][0:3]) for spkezr in obs2venus_spkezrs] - venus_radius
            
            #store minimum and maximum altitudes
            if min_alt > np.min(alts):
                min_alt = np.min(alts)
            if max_alt < np.max(alts):
                max_alt = np.max(alts)

            #get LST - H, m, s
            lsts_hms = [sp.et2lst(et, SPICE_PLANET_ID, lon_rad, SPICE_LONGITUDE_FORM)[0:3] for et, lon_rad in zip(ets, lons_rad)]
            #get LST in hours
            lsts = [lst[0] + lst[1]/60.0 + lst[2]/3600.0 for lst in lsts_hms]
    
    
            
            fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(14, 7), constrained_layout=True)
            scatter1 = ax1.scatter(lons, lats, c=alts, vmin=220, rasterized=True, s=3)
            cbar1 = plt.colorbar(scatter1, ax=ax1)
            cbar1.set_label("Altitude", rotation=270, labelpad=10)

            scatter2 = ax2.scatter(lons, lats, c=lsts, vmin=0, vmax=24, s=3, rasterized=True, cmap=plt.cm.get_cmap("twilight_shifted"))
            cbar2 = plt.colorbar(scatter2, ax=ax2)
            cbar2.set_label("Local solar time", rotation=270, labelpad=10)
            
            title = "%s\nMTP%03i: %s - %s" %(orbit_d[key]["title"], mtp+1, datetime.strftime(dt_start, SHORT_DATETIME_FMT), datetime.strftime(dt_end, SHORT_DATETIME_FMT))
            fig.suptitle(title)
            ax1.set_ylabel("Latitude")
            ax1.set_xlabel("Longitude")
            ax2.set_xlabel("Longitude")
            ax1.grid()
            ax2.grid()
            
            # plt.figure()
            # plt.hist(alts)
            
            pdf.savefig()
            plt.close()
            

        
    sp.kclear()
    
    print(key)
    print("Altitude range=", min_alt, max_alt)