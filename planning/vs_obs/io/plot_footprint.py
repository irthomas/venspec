# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 21:31:08 2022

@author: iant
"""

import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


from planning.vs_obs.filter_selection import filter_parameters
from planning.vs_obs.config.constants import DATETIME_FMT_SECONDS
from planning.vs_obs.observation_cycle import observation_cycle, observation_name

from planning.vs_obs.io.get_venus_magellan_data import get_venus_magellan_data


def plot_footprint(orbit_plan):

    cmap = LinearSegmentedColormap.from_list("", colors=["black", "darkorange", "lemonchiffon"], N=32)

    # plot venus map
    venus_topo = get_venus_magellan_data()
    plt.figure(figsize=(10, 6), constrained_layout=False)
    plt.imshow(venus_topo, extent=(-180, 180, -90, 90), interpolation="bilinear", cmap=cmap, aspect=1, alpha=0.5)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")

    dt_str = datetime.strftime(orbit_plan[0]["dt_start"], DATETIME_FMT_SECONDS)
    plt.title("VenSpec-H groundtrack footprint for %i orbits, starting %s\n%s" % (np.ceil(len(observation_cycle)/2), dt_str, observation_name))

    colours = []
    for orbit in orbit_plan:

        # plot if science
        if orbit["status"] == "science":

            filter_ = orbit["filter"]
            daynight = orbit["daynight"]
            colour = filter_parameters[(filter_, daynight)]["colour"]

            footprint = orbit["footprint"]

            # don't plot if longitude wraps from one side of plot to the other
            if np.max(footprint[:, 0]) - np.min(footprint[:, 0]) < 100:

                # code to show legend label just once
                if colour not in colours:
                    label = "Filter %s" % filter_
                    colours.append(colour)
                else:
                    label = ""
                plt.plot(footprint[:, 0], footprint[:, 1], color=colour, label=label)

        # plot if dark (optional)
        if orbit["status"] == "dark":

            filter_ = "dk"
            daynight = orbit["daynight"]
            colour = filter_parameters[(filter_, daynight)]["colour"]

            footprint = orbit["footprint"]

            # don't plot if longitude wraps from one side of plot to the other
            if np.max(footprint[:, 0]) - np.min(footprint[:, 0]) < 100:

                # code to show legend label just once
                if colour not in colours:
                    label = "Filter %s" % filter_
                    colours.append(colour)
                else:
                    label = ""
                plt.plot(footprint[:, 0], footprint[:, 1], color=colour, label=label)

    plt.legend(loc="lower left")
    plt.grid()
    plt.savefig("vsh_groundtracks_%i_orbits_%s.png" % (len(observation_cycle)/2, dt_str.replace(" ", "_").replace(":", "_")))
