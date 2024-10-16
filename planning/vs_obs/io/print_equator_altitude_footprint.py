# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 14:11:05 2022

@author: iant
"""

import numpy as np


def print_equator_altitude_footprint(orbit_plan):
    """for each semiorbit, find footprint closest to the equator. Print orbit altitude and footprint lat/lon"""
    
    orbit_ixs = []
    for orbit in orbit_plan:
        if "footprint_s" not in orbit[1].keys():
            continue
        
        if orbit[1]["semiorbit"] in orbit_ixs:
            continue
        
        
        footprint = orbit[1]["footprint_s"]
        lat = np.mean(footprint[:, 1])
        if lat > -1 and lat < 1:
            print(orbit[1]["daynight"], orbit[1]["semiorbit"], orbit[1]["alt"])

            footprint_s = orbit[1]["footprint_s"]
            footprint_e = orbit[1]["footprint_e"]
            

            # plt.plot(footprint[:, 0], footprint[:, 1])
            print(footprint)
            

    
            orbit_ixs.append(orbit[1]["semiorbit"])
