# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 21:31:08 2022

@author: iant
"""

import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


from vs_obs.filter_selection import filter_parameters
from vs_obs.config.constants import DATETIME_FMT_SECONDS
from vs_obs.observation_cycle import observation_cycle    

from vs_obs.io.get_venus_magellan_data import get_venus_magellan_data


def plot_footprint(orbit_plan):

    
    cmap = LinearSegmentedColormap.from_list("", colors=["black", "darkorange", "lemonchiffon"], N=32)


    #plot venus map
    venus_topo = get_venus_magellan_data()
    plt.figure(figsize=(10, 6), constrained_layout=True)
    plt.imshow(venus_topo, extent=(-180, 180, -90, 90), interpolation="bilinear", cmap=cmap, aspect="auto")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    
    dt_str = datetime.strftime(orbit_plan[0][0], DATETIME_FMT_SECONDS)
    plt.title("VenSpec-H Groundtracks for %i orbits, starting %s" %(len(observation_cycle)/2, dt_str))
    
    colours = []
    for orbit in orbit_plan:
        if orbit[1]["status"] == "science":
    
            filter_ = orbit[1]["filter"]
            daynight = orbit[1]["daynight"]
            colour = filter_parameters[(filter_, daynight)]["colour"]
            
            footprint_s = orbit[1]["footprint_s"]
            footprint_e = orbit[1]["footprint_e"]
            
            # plt.plot(footprint_s[:, 0], footprint_s[:, 1])
            # plt.plot(footprint_e[:, 0], footprint_e[:, 1])
    
            #make footprint 
            footprint = np.zeros((5, 2))
            footprint[0, :] = footprint_s[1, :]
            footprint[1, :] = footprint_e[2, :]
            footprint[2, :] = footprint_e[3, :]
            footprint[3, :] = footprint_s[0, :]
            footprint[4, :] = footprint_s[1, :]
            
            
            
            if np.max(footprint[:, 0]) - np.min(footprint[:, 0]) < 100:
                if colour not in colours:
                    label = "Filter %s" %filter_
                    colours.append(colour)
                else:
                    label = ""
                plt.plot(footprint[:, 0], footprint[:, 1], color=colour, label=label)
                
    plt.legend()
    plt.grid()
    plt.savefig("vsh_groundtracks_%i_orbits_%s.png" %(len(observation_cycle)/2, dt_str.replace(" ","_").replace(":","_")))
    
