# -*- coding: utf-8 -*-
"""
Created on Tue May 24 11:34:33 2022

@author: iant
"""


import matplotlib.pyplot as plt
import numpy as np
import matplotlib.collections as mcoll

def plot_coloured_line(ax, x, y, z=None, cmap=plt.get_cmap('viridis_r'), norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0):
    """make line plot with variable colours based on z"""

    segments = make_segments(np.asarray(x), np.asarray(y))
    lc = mcoll.LineCollection(segments, array=np.asarray(z), cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)

    ax.add_collection(lc)

    ax.set_ylim([-90, 90])
    ax.set_xlim([-180, 180])

    return lc


def make_segments(x, y):

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

