# -*- coding: utf-8 -*-
"""
Created on Thu May 12 09:15:54 2022

@author: iant

FOV CALCS
"""

import numpy as np
import spiceypy as sp
from datetime import datetime, timedelta

from tools.file.write_log import write_log
from tools.spice.load_spice_kernels import load_spice_kernels, SPICE_DATETIME_FMT, \
    SPICE_METHOD, SPICE_ABCORR, SPICE_SHAPE_MODEL_METHOD

from tools.general.progress import progress

#load spice kernels
orbit_name = "EnVision_ESC_T2_2032_SouthVOI"
# orbit_name = "EnVision_HEO_T2_2032_NorthVOI"
orbit_dict = load_spice_kernels("%s.bsp" %orbit_name)


#get integration time depending on observation type
obs_settings = {
    "dayside":{"integration_time":1.0},
    "nightside":{"integration_time":13.47},
    }



fov_dimensions = np.array([6.7, 0.1333]) * np.pi / 180.0 #degrees to radians


angular_error = 25 / 1000. #86 arcminutes in mrad, converted to rad. Value for 14.4 seconds




    
print("SPICE calculations")

#get science phase start/end ephemeris times
dt_start = orbit_dict["science_start"]
dt_end = dt_start + timedelta(days=112)
et_start = sp.utc2et(datetime.strftime(dt_start, SPICE_DATETIME_FMT))
et_end = sp.utc2et(datetime.strftime(dt_end, SPICE_DATETIME_FMT))

#make list of ephemeris times        
ets = np.arange(et_start, et_end, 300.0) #one value per X seconds

#make strings
dts = [sp.et2utc(et, "C", 0) for et in ets]

n_points = len(ets)

#spice calculations
venus_radius = sp.bodvrd("VENUS", "RADII", 3)[1][0] #no ellipsoid for Venus -> circular so only one radius value

#sub-observer point data
observer = orbit_dict["spice_observer"]
subpnts = [sp.subpnt(SPICE_METHOD, "VENUS", et, "IAU_VENUS", SPICE_ABCORR, observer) for et in ets]
subpnts_xyz = [subpnt[0] for subpnt in subpnts]


#convert to lat/lons in degrees
reclats = [sp.reclat(subpnt_xyz) for subpnt_xyz in subpnts_xyz]
lons_rad = [reclat[1] for reclat in reclats]
lats_rad = [reclat[2] for reclat in reclats]
lons = np.asfarray(lons_rad) * sp.dpr()
lats = np.asfarray(lats_rad) * sp.dpr()

#get orbit altitude
# find obs position/velocity rel to mars in J2000
obs2venus_spkezrs = [sp.spkezr("VENUS", et, "IAU_VENUS", SPICE_ABCORR, observer) for et in ets]

#height of observer above Mars centre
altitudes = np.array([sp.vnorm(spkezr[0][0:3]) for spkezr in obs2venus_spkezrs] - venus_radius)

#spacecraft velocities
speeds = np.array([sp.vnorm(spkezr[0][3:6]) for spkezr in obs2venus_spkezrs])

#surface speeds
circum_venus = 2.0 * np.pi * (venus_radius)
circum_orbit = 2.0 * np.pi * (venus_radius + altitudes)

surf_speeds = speeds * circum_venus / circum_orbit

#incidence angles
surf_ilumin = [sp.ilumin(SPICE_SHAPE_MODEL_METHOD, "VENUS", et, "IAU_VENUS", SPICE_ABCORR, observer, subpnt_xyz) for et, subpnt_xyz in zip(ets, subpnts_xyz)]
incidence_angles = [ilumin[3] * sp.dpr() for ilumin in surf_ilumin]


#write to file
write_log("Time\tObservation Type\tLongitude\tLatitude\tAltitude\tSpeed\tSurface Speed\tSurface FOV across-track\tSurface FOV along-track\tSurface FOV along-track error\tSurface FOV along-track with error", "%s_fovs.txt" %orbit_name)

max_values = {"fov along":0, "fov+error along":0, "error along":0}
lines = []
for i in progress(range(n_points)):
    
    
    if incidence_angles[i] < 90.0:
        obs_type = "dayside"
    else:
        obs_type = "nightside"
    
    obs_params = obs_settings[obs_type]
    distance = surf_speeds[i] * obs_params["integration_time"]

    angular_error_over_it = angular_error * obs_params["integration_time"] / 14.4 #ratio of error to 14.4s as defined as the RPE
    
    ifov = 2.0 * np.tan(fov_dimensions / 2.0) * altitudes[i]
    fov_error = np.tan(angular_error_over_it) * altitudes[i]

    #calculate field of view including s/c velocity
    fov = ifov + distance
   
    #write to file
    line = "%s\t%s\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f" \
        %(dts[i], obs_type, lons[i], lats[i], altitudes[i], speeds[i], surf_speeds[i], ifov[0], fov[1], fov_error, fov[1]+fov_error)
        
    
    lines.append(line)

    if fov[1] > max_values["fov along"]:
        max_values["fov along"] = fov[1]
    if fov[1]+fov_error > max_values["fov+error along"]:
        max_values["fov+error along"] = fov[1]+fov_error
    if fov_error > max_values["error along"]:
        max_values["error along"] = fov_error


print("Writing to file")
for line in lines:
    write_log(line, "%s_fovs.txt" %orbit_name)
        
      
print(orbit_name)
print("Max values:")
for key in max_values:
    print(key, ":", max_values[key])
