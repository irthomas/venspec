# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:34:09 2022

@author: iant
"""

import numpy as np

from vs_obs.config.constants import FILTER_WHEEL_LAYOUT, ANGLE_BETWEEN_FILTERS, WHEEL_ROTATION_SPEED


def calculate_angle(pos1, pos2):
    """calculate angle from start and end filter position (degrees)"""
    pos1_i = FILTER_WHEEL_LAYOUT.index(pos1)
    pos2_i = FILTER_WHEEL_LAYOUT.index(pos2)
    if np.abs(pos2_i - pos1_i) < (len(FILTER_WHEEL_LAYOUT)/2. + 1.):
        step = np.abs(pos2_i - pos1_i)
    else:
        step = len(FILTER_WHEEL_LAYOUT) - np.abs(pos2_i - pos1_i)
    
    return step * ANGLE_BETWEEN_FILTERS



def calculate_wheel_rotation_time(angle):
    """calculate movement time for a single movement (seconds)"""
    movement_time = angle / WHEEL_ROTATION_SPEED
    return movement_time


