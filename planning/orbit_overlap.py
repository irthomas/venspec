# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 13:37:37 2022

@author: iant
"""
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 10:13:18 2022

@author: iant

CHECK ENVISION KERNELS FOR
ORBIT OVERLAP
"""

import os
import spiceypy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime, timedelta

from matplotlib.backends.backend_pdf import PdfPages


from tools.file.paths import paths
from tools.plotting.plot_coloured_line import plot_coloured_line

alpha = 1.0


# SPICE_FORMATSTR = "C"
# SPICE_PRECISION = 0
SPICE_DATETIME_FMT = "%Y %b %d %H:%M:%S.%f"
SHORT_DATETIME_FMT = "%d %b %Y"

LONG_DATETIME_FMT = "%Y %b %d %H:%M:%S"

SPICE_METHOD = "INTERCEPT/ELLIPSOID"
SPICE_SHAPE_MODEL_METHOD = "Ellipsoid"
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
        "title":"ALT_T1_2032_NorthVOI (2022/01/18)\nAlternative backup scenario with Venus Orbit Insertion\nover the Northern hemisphere",
        "science_start":datetime(2034, 12, 20, 0, 44, 15),
        "science_end":datetime(2038, 12, 18),
        "spice_observer":"-999", #error in kernel production => should be -668
        },
    # "EnVision_ALT_T4_2032_SouthVOI.bsp":{
    #     "title":"ALT_T4_2032_SouthVOI (2022/01/18)\nAlternative baseline scenario with Venus Orbit Insertion over the Southern hemisphere",
    #     "science_start":datetime(2035, 3, 27),
    #     "science_end":datetime(2039, 3, 25),
    #     "spice_observer":"-999", #error in kernel production => should be -668
    #     },
    # "EnVision_ESC_T4_2032_NorthVOI.bsp":{
    #     "title":"ESC_T4_2032_NorthVOI (2020/06/10)\nLaunch into direct escape at -5.0 deg declination and Venus Orbit Insertion over the Northern hemisphere.",
    #     "science_start":datetime(2035, 6, 15),
    #     # "science_end":datetime(2035, 9, 14),
    #     "science_end":datetime(2039, 6, 14),
    #     "spice_observer":"-668",
    #     },
    # "EnVision_ESC_T2_2032_SouthVOI.bsp":{
    #     "title":"ESC_T2_2032_SouthVOI (2019/07/26)\nLaunch into direct escape at -3.0 deg declination and Venus Orbit Insertion over the Southern hemisphere",
    #     "science_start":datetime(2034, 11, 26),
    #     "science_end":datetime(2040, 3, 24),
    #     "spice_observer":"-668",
    #     },
    # "EnVision_HEO_T2_2032_NorthVOI.bsp":{
    #     "title":"HEO_T2_2032_NorthVOI (2019/07/26)\nLaunch into HEO (310x100,000 km) and Venus Orbit Insertion over the Northern hemisphere",
    #     "science_start":datetime(2034, 11, 26),
    #     "science_end":datetime(2040, 3, 24),
    #     "spice_observer":"-668",
    #     },
    }



for key in orbit_d.keys():
    for kernel in kernels:
        sp.furnsh(os.path.join(kernel_root_dir, os.path.normcase(kernel)))
        
    sp.furnsh(os.path.join(kernel_root_dir, "spk", key))
    
    
    
    obs_length = 92. * 4.5 #length of observation sequence (minutes)
    obs_start_delta = 24. * 60. #time between observation sequence start times (minutes)
    n_starts = 5. #number of observation sequences

    with PdfPages("%s.pdf" %key) as pdf:
        
        for days in range(45, 56+45, 1):
            science_start = orbit_d[key]["science_start"] + timedelta(days=days)
        
            
            dt_starts = [science_start + timedelta(minutes=obs_start_delta*i) for i in np.arange(n_starts)]
        
            fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(8, 6), constrained_layout=True)
        
            lons_all = []
            lats_all = []
            incidence_all = []
        
            #loop through orbits
            for i, dt_start in  enumerate(dt_starts):
                
            
                #get science phase start/end ephemeris times
                dt_end = dt_start + timedelta(minutes=obs_length)
                et_start = sp.utc2et(datetime.strftime(dt_start, SPICE_DATETIME_FMT))
                et_end = sp.utc2et(datetime.strftime(dt_end, SPICE_DATETIME_FMT))
                
                #make list of ephemeris times        
                ets = np.arange(et_start, et_end, 2.0) #one value per 1 second
                
                
                #spice calculations
                #sub point data
                subpnts = [sp.subpnt(SPICE_METHOD, "VENUS", et, "IAU_VENUS", SPICE_ABCORR, orbit_d[key]["spice_observer"]) for et in ets]
                subpnts_xyz = [subpnt[0] for subpnt in subpnts]
                
                
                #convert to lat/lons
                reclats = [sp.reclat(subpnt_xyz) for subpnt_xyz in subpnts_xyz]
                lons_rad = [reclat[1] for reclat in reclats]
                lats_rad = [reclat[2] for reclat in reclats]
                lons = np.asfarray(lons_rad) * sp.dpr()
                lats = np.asfarray(lats_rad) * sp.dpr()
                
                #incidence angles
                surf_ilumin = [sp.ilumin(SPICE_SHAPE_MODEL_METHOD, "VENUS", et, "IAU_VENUS", SPICE_ABCORR, orbit_d[key]["spice_observer"], subpnt_xyz) for et, subpnt_xyz in zip(ets, subpnts_xyz)]
                incidence_angles = np.asarray([ilumin[3] * sp.dpr() for ilumin in surf_ilumin])
        
                # incidence_angles_scaled = (incidence_angles - np.min(incidence_angles)) / (np.max(incidence_angles) - np.min(incidence_angles))
        
                #plot only when contiguous i.e. avoid longitude wrapping
                lons_diff = np.abs(np.diff(lons))
                wrap_ixs = np.where(lons_diff > 150)[0]
                
                previous_wrap_ix = 0
                #loop through contiguous data
                for wrap_ix in wrap_ixs:
                    lons_all.append(lons[previous_wrap_ix:wrap_ix])
                    lats_all.append(lats[previous_wrap_ix:wrap_ix])
                    incidence_all.append(incidence_angles[previous_wrap_ix:wrap_ix])
                    
                    previous_wrap_ix = wrap_ix + 1
        
        
                # plot end of file after final wrap
                lons_all.append(lons[previous_wrap_ix:])
                lats_all.append(lats[previous_wrap_ix:])
                incidence_all.append(incidence_angles[previous_wrap_ix:])
                
            # vmin = np.min([np.min(i) for i in incidence_all])
            # vmax = np.max([np.max(i) for i in incidence_all])
            
            vmin = 0.
            vmax = 180.
                
            for lons, lats, incidences in zip(lons_all, lats_all, incidence_all):
                plot_coloured_line(ax1, lons, lats, (incidences-vmin)/(vmax-vmin))
                
            norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
            sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=norm)
            sm.set_array([])
        
            cbar = fig.colorbar(sm)
            cbar.set_label("Solar Incidence Angle", rotation=270, labelpad=10)
            
            # title = "%s\n%s - %s" %(orbit_d[key]["title"], datetime.strftime(dt_starts[0], LONG_DATETIME_FMT), datetime.strftime(dt_end, LONG_DATETIME_FMT))
            title = "%s - %s" %(datetime.strftime(dt_starts[0], LONG_DATETIME_FMT), datetime.strftime(dt_end, LONG_DATETIME_FMT))
            fig.suptitle(title)
            ax1.set_ylabel("Latitude")
            ax1.set_xlabel("Longitude")
            ax1.grid()
            
            stop()
            
            # fig.savefig("EnVision_orbit_incidence_angle_day%i.png" %days)
            pdf.savefig()
            plt.close()
            

    
    sp.kclear()
    
    # print(key)
    # print("Altitude range=", min_alt, max_alt)