# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:02:16 2024

@author: iant

PLOT ENVISION TRAJECTORY - NOTE THAT KERNEL envision_study_et1_north_voi_ml008_v050.tm
DOES NOT YET CONTAIN TRAJECTORY INFO
"""

import numpy as np
import spiceypy as sp
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

from planning.vs_obs.other.load_spice_kernels import load_spice_kernels
load_spice_kernels()


dt_start = datetime(2033, 6, 1)
t_delta = 1.0

ets = [sp.utc2et(str(dt_start + timedelta(seconds=i))) for i in range(1000)]

pos = np.zeros((len(ets), 3))

for i, et in enumerate(ets):
    pos[i, :] = sp.spkpos("-668", et, "J2000", "None", "Sun")[0]


ax = plt.figure().add_subplot(projection='3d')

ax.plot(pos[:, 0], pos[:, 1], pos[:, 2])
