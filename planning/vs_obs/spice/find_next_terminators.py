# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 16:26:24 2024

@author: iant
"""

from planning.vs_obs.spice.spice_functions import et2dt, dt2et, get_centre_sza


def next_d2n_terminator(dt):
    """find first day-to-night terminator"""
    et = dt2et(dt)
    incidence_angle = get_centre_sza(et)

    # if already on nightside
    if incidence_angle > 90:
        while incidence_angle > 90:
            et += 1. * 60.
            incidence_angle = get_centre_sza(et)
            # print(incidence_angle)
    # print("**")
    while incidence_angle < 90:
        et += 1
        incidence_angle = get_centre_sza(et)
    # print(incidence_angle)

    return et2dt(et)


def next_n2d_terminator(dt):
    """find first day-to-night terminator"""
    et = dt2et(dt)
    incidence_angle = get_centre_sza(et)

    # if already on dayside
    if incidence_angle < 90:
        while incidence_angle < 90:
            et += 1. * 60.
            incidence_angle = get_centre_sza(et)
            # print(incidence_angle)
    # print("**")
    while incidence_angle > 90:
        et += 1
        incidence_angle = get_centre_sza(et)
    # print(incidence_angle)

    return et2dt(et)
