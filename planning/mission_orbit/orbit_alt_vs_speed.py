# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 10:27:56 2022

@author: iant
"""


import os
import numpy as np
import spiceypy as sp
from datetime import datetime, timedelta

import matplotlib.pyplot as plt

# from tools.file.write_log import write_log

# from tools.general.progress import progress
from tools.file.paths import paths

# # SPICE_FORMATSTR = "C"
# # SPICE_PRECISION = 0
# SPICE_DATETIME_FMT = "%Y %b %d %H:%M:%S.%f"
# SHORT_DATETIME_FMT = "%d %b %Y"

# LONG_DATETIME_FMT = "%Y %b %d %H:%M:%S"

# SPICE_METHOD = "INTERCEPT/ELLIPSOID"
# SPICE_SHAPE_MODEL_METHOD = "Ellipsoid"
# SPICE_ABCORR = "NONE"
# # SPICE_OBSRVR = "-668"

# SPICE_LONGITUDE_FORM = "PLANETOCENTRIC"
# SPICE_PLANET_ID = 299 # venus

# kernel_root_dir = paths["KERNEL_ROOT_DIRECTORY"]

# #load each kernel manually

# kernels = [
#     "fk/envision_v00.tf",
#     "fk/envision_sc_venus_npo_v00.tf",
#     "lsk/naif0012.tls",
#     "pck/pck00010.tpc",
#     "pck/de-403-masses.tpc",
#     "spk/de432s.bsp",
#     ]

# orbit_d ={
#     "EnVision_ALT_T1_2032_NorthVOI.bsp":{
#         "title":"ALT_T1_2032_NorthVOI (2022/01/18)\nAlternative backup scenario with Venus Orbit Insertion\nover the Northern hemisphere",
#         "science_start":datetime(2034, 12, 20, 0, 44, 15),
#         "science_end":datetime(2038, 12, 18),
#         "spice_observer":"-999", #error in kernel production => should be -668
#         },
#     "EnVision_ALT_T4_2032_SouthVOI.bsp":{
#         "title":"ALT_T4_2032_SouthVOI (2022/01/18)\nAlternative baseline scenario with Venus Orbit Insertion over the Southern hemisphere",
#         "science_start":datetime(2035, 3, 27),
#         "science_end":datetime(2039, 3, 25),
#         "spice_observer":"-999", #error in kernel production => should be -668
#         },
#     "EnVision_ESC_T4_2032_NorthVOI.bsp":{
#         "title":"ESC_T4_2032_NorthVOI (2020/06/10)\nLaunch into direct escape at -5.0 deg declination and Venus Orbit Insertion over the Northern hemisphere.",
#         "science_start":datetime(2035, 6, 15),
#         # "science_end":datetime(2035, 9, 14),
#         "science_end":datetime(2039, 6, 14),
#         "spice_observer":"-668",
#         },
#     "EnVision_ESC_T2_2032_SouthVOI.bsp":{
#         "title":"ESC_T2_2032_SouthVOI (2019/07/26)\nLaunch into direct escape at -3.0 deg declination and Venus Orbit Insertion over the Southern hemisphere",
#         "science_start":datetime(2034, 11, 26),
#         "science_end":datetime(2040, 3, 24),
#         "spice_observer":"-668",
#         },
#     "EnVision_HEO_T2_2032_NorthVOI.bsp":{
#         "title":"HEO_T2_2032_NorthVOI (2019/07/26)\nLaunch into HEO (310x100,000 km) and Venus Orbit Insertion over the Northern hemisphere",
#         "science_start":datetime(2034, 11, 26),
#         "science_end":datetime(2040, 3, 24),
#         "spice_observer":"-668",
#         },
#     }


# plt.figure(figsize=(10, 8))
# plt.title("EnVision altitude vs velocity for different SPICE kernels")

# for key in orbit_d.keys():
#     for kernel in kernels:
#         sp.furnsh(os.path.join(kernel_root_dir, os.path.normcase(kernel)))

#     sp.furnsh(os.path.join(kernel_root_dir, "spk", key))


#     print("SPICE calculations", key)

#     #get science phase start/end ephemeris times
#     dt_start = orbit_d[key]["science_start"]
#     dt_end = dt_start + timedelta(days=200)
#     et_start = sp.utc2et(datetime.strftime(dt_start, SPICE_DATETIME_FMT))
#     et_end = sp.utc2et(datetime.strftime(dt_end, SPICE_DATETIME_FMT))

#     #make list of ephemeris times
#     ets = np.arange(et_start, et_end, 1200.0) #one value per X seconds

#     #make strings
#     dts = [sp.et2utc(et, "C", 0) for et in ets]

#     n_points = len(ets)

#     #spice calculations
#     venus_radius = sp.bodvrd("VENUS", "RADII", 3)[1][0] #no ellipsoid for Venus -> circular so only one radius value

#     #sub-observer point data
#     observer = orbit_d[key]["spice_observer"]
#     subpnts = [sp.subpnt(SPICE_METHOD, "VENUS", et, "IAU_VENUS", SPICE_ABCORR, observer) for et in ets]
#     subpnts_xyz = [subpnt[0] for subpnt in subpnts]


#     #convert to lat/lons in degrees
#     reclats = [sp.reclat(subpnt_xyz) for subpnt_xyz in subpnts_xyz]
#     lons_rad = [reclat[1] for reclat in reclats]
#     lats_rad = [reclat[2] for reclat in reclats]
#     lons = np.asfarray(lons_rad) * sp.dpr()
#     lats = np.asfarray(lats_rad) * sp.dpr()

#     #get orbit altitude
#     # find obs position/velocity rel to venus in venus frame
#     obs2venus_spkezrs = [sp.spkezr("VENUS", et, "IAU_VENUS", SPICE_ABCORR, observer) for et in ets]

#     #height of observer above Mars centre
#     altitudes = [sp.vnorm(spkezr[0][0:3]) for spkezr in obs2venus_spkezrs] - venus_radius

#     #velocities
#     speeds = [sp.vnorm(spkezr[0][3:6]) for spkezr in obs2venus_spkezrs]


#     if key == "EnVision_ESC_T2_2032_SouthVOI.bsp":
#         with open("altitude_vs_velocity.tsv", "w") as f:
#             line = "Spacecraft Altitude (km)\tSpacecraft Velocity (km/s)\tGenerated with SPICE kernel %s\n" %key
#             f.write(line)
#             for altitude, speed in zip(altitudes, speeds):
#                 line = "%0.2f\t%0.3f\n" %(altitude, speed)
#                 f.write(line)

#     plt.scatter(altitudes, speeds, label=key)

# OIP circular orbit calc
gsd = 100e3  # m
r = 6.051e6  # m
h = 470.0e3  # m
# e = np.arange(250, 500, 10) * 1e3 #m
g = 6.674e-11  # m3 kg−1 s−2
m = 4.87e24

# it = (gsd / r) * np.sqrt(((r + h)**3) / (g * m))

# v = gsd / it / 1e3
# print(v)


# v = np.sqrt(g * m / (r + h))
circum_v = 2.0 * np.pi * (r + 0)
circum_e = 2.0 * np.pi * (r + h)

v = 6.98  # km/s
print(v * 14.401)

ground_track = v * 14.401 / circum_e * circum_v

print(ground_track)


# plt.plot(h / 1e3, v, color="k", label="OIP formula")


# plt.xlabel("Altitude above surface (km)")
# plt.ylabel("Spacecraft velocity (km per second)")
# plt.grid()
# plt.legend()
