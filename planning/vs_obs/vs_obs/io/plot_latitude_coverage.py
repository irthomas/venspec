# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 13:57:39 2022

@author: iant
"""

import matplotlib.pyplot as plt
import numpy as np

from vs_obs.observation_cycle import observation_cycle


def plot_latitude_coverage(coverage_grids, lats_hr, lons_hr, title):
    
    #find number of filters with data    
    filters = []
    for filter_ in coverage_grids.keys():
        if not np.all(coverage_grids[filter_] == 0.0):
            filters.append(filter_)
    filters = sorted(filters)

    
    fig, axes = plt.subplots(figsize=(8, 5), ncols=len(filters), sharey=True)
    fig.suptitle("VenSpec-H longitudinal coverage for %i orbits\n%s" %(np.ceil(len(observation_cycle)/2), title))

    # n_lats = len(coverage_grids["1"][:, 0])
    n_lons = len(coverage_grids["1"][0, :])
    
    for i, filter_ in enumerate(filters):
    
        coverage_lons = np.sum(coverage_grids[filter_], axis=1) / n_lons * 100
        axes[i].barh(lats_hr, coverage_lons)
        axes[i].set_title("Filter %s" %filter_)
        axes[i].set_xlabel("Coverage percentage %")
        axes[i].grid()
    axes[0].set_ylabel("Latitude")

    return fig