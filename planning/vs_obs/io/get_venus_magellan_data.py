# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 16:06:05 2022

@author: iant

"""

import os
import numpy as np

from planning.vs_obs.config.paths import paths


def get_venus_magellan_data():

    path = os.path.join(paths["REFERENCE_DIRECTORY"], "venus_magellan_map.txt")

    return np.loadtxt(path)
