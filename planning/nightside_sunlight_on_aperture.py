# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 20:01:12 2023

@author: iant

SOLAR INCIDENCE ANGLES ON VENSPEC-H APERTURE
"""



import os
import spiceypy as sp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from datetime import datetime, timedelta



from tools.file.paths import paths


plot = []
# plot += ["3d"]
plot += ["angles"]

SPICE_DATETIME_FMT = "%Y %b %d %H:%M:%S.%f"
SHORT_DATETIME_FMT = "%d %b %Y"

LONG_DATETIME_FMT = "%Y %b %d %H:%M:%S"

SPICE_METHOD = "INTERCEPT/ELLIPSOID"
SPICE_SHAPE_MODEL_METHOD = "Ellipsoid"

SP_DPR = sp.dpr()


SPICE_PLANET_REFERENCE_FRAME = "IAU_VENUS"
SPICE_ABERRATION_CORRECTION = "None"
SPICE_PLANET_ID = 299
# et2lst: form of longitude supplied by the variable lon
SPICE_LONGITUDE_FORM = "PLANETOCENTRIC"
# spkpos: reference frame relative to which the output position vector
# should be expressed
SPICE_REFERENCE_FRAME = "J2000"
#et2utc: string format flag describing the output time string. 'C' Calendar format, UTC
SPICE_STRING_FORMAT = "C"
# et2utc: number of decimal places of precision to which fractional seconds
# (for Calendar and Day-of-Year formats) or days (for Julian Date format) are to be computed
SPICE_TIME_PRECISION = 3

SPICE_TARGET = "VENUS"


SPICE_SHAPE_MODEL_METHOD = "Ellipsoid"
SPICE_INTERCEPT_METHOD = "INTERCEPT/ELLIPSOID"

# DREF = "TGO_NOMAD_LNO_OPS_NAD"







kernel_root_dir = paths["KERNEL_ROOT_DIRECTORY"]

#load each kernel manually

kernels = [
    "fk/envision_v00.tf",
    "envision_venspec_v01.ti",
    "lsk/naif0012.tls",
    "pck/pck00010.tpc",
    "pck/de-403-masses.tpc",
    "spk/de432s.bsp",
    ]

orbit_d ={
    "EnVision_ET1_2031_NorthVOI.bsp":{
        "title":"EnVision_ET1_2031_NorthVOI (2022/01/18)\nNominal orbit with Venus Orbit Insertion\nover the Northern hemisphere",
        "science_start":datetime(2035, 3, 18, 0, 0, 0),
        "science_end":datetime(2038, 1, 1),
        "spice_observer":"-668", #error in kernel production => should be -668
        },
}

key = "EnVision_ET1_2031_NorthVOI.bsp"
for kernel in kernels:
    sp.furnsh(os.path.join(kernel_root_dir, os.path.normcase(kernel)))
        
sp.furnsh(os.path.join(kernel_root_dir, "spk", key))


SPICE_OBSERVER = orbit_d[key]["spice_observer"]
    
    

# start time close to anti-solar point
# offset_days = 70
# seconds  = range(700, 6100, 10)

# start of mission
# offset_days = 0
# seconds  = range(0, 50000000, 100)

#end of mission
offset_days = 365*2
seconds  = range(0, 50000000, 100)

    


dt_starts = [orbit_d[key]["science_start"] + timedelta(seconds=i) + timedelta(days=offset_days) for i in seconds]
et = [sp.utc2et(datetime.strftime(dt_start, SPICE_DATETIME_FMT)) for dt_start in dt_starts]



dref = "ENVISION_SPACECRAFT"
bodyAxes = sp.bodvrd("VENUS", "RADII", 3)[1]




if "3d" in plot:
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect(aspect = (1,1,1))

    ax_scaler = 10000
    ax.axes.set_xlim3d(left=-ax_scaler, right=ax_scaler) 
    ax.axes.set_ylim3d(bottom=-ax_scaler, top=ax_scaler) 
    ax.axes.set_zlim3d(bottom=-ax_scaler, top=ax_scaler) 

# target, time, ref, abs corr, observer
ven_pos=np.asfarray([sp.spkpos("VENUS", et, "J2000", "NONE", "-668")[0] for et in et])
sun_pos=np.asfarray([sp.spkpos("SUN", et, "J2000", "NONE", "-668")[0] for et in et])# / 1.0e4

ven_norm = np.asfarray([ven_pos[i, :] / np.linalg.norm(ven_pos[i, :]) for i in range(len(et))])
sun_norm = np.asfarray([sun_pos[i, :] / np.linalg.norm(sun_pos[i, :]) for i in range(len(et))])

dots = np.asfarray([np.dot(ven_norm[i, :], sun_norm[i, :]) for i in range(len(et))])
angles = np.arccos(dots) * SP_DPR


# #%%
# fig = plt.figure(figsize=(12,12))
# ax = fig.add_subplot(111, projection='3d')
# ax.set_box_aspect(aspect = (1,1,1))
# ax_scaler = 1.1
# ax.axes.set_xlim3d(left=-ax_scaler, right=ax_scaler) 
# ax.axes.set_ylim3d(bottom=-ax_scaler, top=ax_scaler) 
# ax.axes.set_zlim3d(bottom=-ax_scaler, top=ax_scaler) 

# ax.scatter(sun_norm[:, 0], sun_norm[:, 1], sun_norm[:, 2])
# ax.scatter(ven_norm[:, 0], ven_norm[:, 1], ven_norm[:, 2])
# ax.scatter(0., 0., 0.)
# #%%
# stop()

i = 0

if "3d" in plot:
    plt.plot([0, ven_pos[i,0]], [0, ven_pos[i,1]], [0, ven_pos[i,2]])
    plt.plot([0, sun_pos[i,0]], [0, sun_pos[i,1]], [0, sun_pos[i,2]])
    
    u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
    x = np.cos(u)*np.sin(v) * bodyAxes[0] + ven_pos[i,0]
    y = np.sin(u)*np.sin(v) * bodyAxes[0] + ven_pos[i,1]
    z = np.cos(v) * bodyAxes[0] + ven_pos[i,2]
    # ax.plot_wireframe(x, y, z, color="r")
    ax.plot_surface(x, y, z, color="r", alpha=0.3)

    print(sun_pos[i, :])

angle_d = {"et":[], "angle":[]}

for i in range(len(et)):
    #distance between line [p1 -> p2] and point [p0]
    p1 = np.array([0., 0., 0.])
    p2 = sun_pos[i, :]
    p0 = ven_pos[i, :]
    
    tangent_alt = np.linalg.norm(np.cross((p0-p1), (p0-p2))) / np.linalg.norm((p2-p1)) - bodyAxes[0]
    
    
    if tangent_alt > 0.0 and angles[i] < 90.0:
        # angle_et.append(et[i])
        
        angle_d["et"].append(dt_starts[i])
        angle_d["angle"].append(angles[i])



if "angles" in plot:
    plt.figure(figsize=(15, 6))
    plt.scatter(angle_d["et"], angle_d["angle"], s=2)
    plt.xlabel("Datetime")
    plt.ylabel("Angle between EnVision->Venus vector and EnVision->Sun vector (degrees)")
    plt.title("%s\nAngle between EnVision->Venus vector and EnVision->Sun vector where Sun vector does not intersect planet" %key)
    plt.grid()
    