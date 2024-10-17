# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:34:09 2022

@author: iant
"""

import numpy as np

from planning.vs_obs.config.constants import FILTER_WHEEL_LAYOUT, ANGLE_BETWEEN_FILTERS, WHEEL_ROTATION_SPEED, DETECTOR_READOUT_TIME
from planning.vs_obs.instrument.filter_selection import filter_sequences, filter_parameters


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


def make_sequence(filter_sequence, daynight, duration, silent=False):
    """find time prior to sza crossing for dark frame and rotation if required"""
    et = 0.0  # relative to sza crossing time. At this time, the first filter is already in position ready to run
    seq = []

    # go backwards in time to calculate first dark
    first_dark_name = "dk"
    first_filter_name = filter_sequence[0]

    filter_wheel_angle = calculate_angle(first_dark_name, first_filter_name)
    rotation_time = calculate_wheel_rotation_time(filter_wheel_angle)
    int_time = filter_parameters[("dk", daynight)]["integration_time"]

    if not silent:
        print("First dark frame start", et - rotation_time - int_time)
    seq.append([et - rotation_time - int_time, first_dark_name, "int start", "Start frame acquisition filter %s" % first_dark_name])
    if not silent:
        print("First dark frame end", et - rotation_time)
    seq.append([et - rotation_time, first_dark_name, "int end", "End frame acquisition filter %s" % first_dark_name])
    if not silent:
        print("Start filter wheel rotation %s->%s" % (first_dark_name, first_filter_name), et - rotation_time)
    seq.append([et - rotation_time, first_dark_name, "rot start", "Start filter wheel rotation %s->%s" %
                (first_dark_name, first_filter_name)])
    if not silent:
        print("Filter %s now in position" % (first_filter_name), et)
    seq.append([et, first_filter_name, "rot end", "End filter wheel rotation %s->%s" %
                (first_dark_name, first_filter_name)])

    """find complete filter sequences until observation time limit is reached"""
    et = 0.0  # relative to sza crossing time. At this time, the first filter is already in position ready to run
    # seq = []

    i = 0  # loop number (filter sequence)
    while et < duration:

        # get the next filter in the sequence. Return to first filter at the end of the sequence
        current_filter_name = filter_sequence[i % len(filter_sequence)]
        next_filter_name = filter_sequence[(i + 1) % len(filter_sequence)]

        filter_wheel_angle = calculate_angle(current_filter_name, next_filter_name)
        rotation_time = calculate_wheel_rotation_time(filter_wheel_angle)
        int_time = filter_parameters[(current_filter_name, daynight)]["integration_time"]

        if not silent:
            print("Filter %s frame start" % current_filter_name, et)
        seq.append([et, current_filter_name, "int start", "Start frame acquisition filter %s" % current_filter_name])
        if not silent:
            print("Filter %s frame end" % current_filter_name, et + int_time)
        seq.append([et + int_time, current_filter_name, "int end", "End frame acquisition filter %s" % current_filter_name])

        if rotation_time > 0.0:
            if not silent:
                print("Start filter wheel rotation %s->%s" % (current_filter_name, next_filter_name), et + int_time)
            seq.append([et + int_time, current_filter_name, "rot start", "Start filter wheel rotation %s->%s" %
                        (current_filter_name, next_filter_name)])
            if not silent:
                print("Filter %s now in position" % (next_filter_name), et + int_time + rotation_time)
            seq.append([et + int_time + rotation_time, next_filter_name, "rot end", "End filter wheel rotation %s->%s" %
                        (current_filter_name, next_filter_name)])
            et += int_time + rotation_time

        else:  # if no rotation
            if not silent:
                print("No filter wheel rotation, filter %s readout start" % current_filter_name, et + int_time)
            seq.append([et + int_time, current_filter_name, "readout start", "No filter wheel rotation, filter %s readout start" % current_filter_name])
            if not silent:
                print("No filter wheel rotation, filter %s readout end" % current_filter_name, et + int_time + DETECTOR_READOUT_TIME)
            seq.append([et + int_time + DETECTOR_READOUT_TIME, current_filter_name, "readout start",
                       "No filter wheel rotation, filter %s readout start" % current_filter_name])
            et += int_time + DETECTOR_READOUT_TIME

        i += 1

    # once end time is reached, check what filter is running
    # the time must have been reached during detector acquisition or filter wheel rotation/readout
    # if during rotation/readout, then keep the measurement and do a final dark (even if dark)
    # if during detector acq, then stop immediately, move to dark, measure dark

    # get sequences where duration is exceeded
    seqs_over = [each_seq for each_seq in seq if each_seq[0] > duration]

    if seqs_over[0][2] == "int end":
        # end time reached during detector acquisition
        # stop observation immediately, then move to dark (if not there already) and measure dark frame

        last_filter_name = seqs_over[0][1]
        last_dark_name = "dk"
        et = duration

        # remove final rotation / readout from the sequence (never reaches this point)
        seq.pop(-1)
        seq.pop(-1)
        seq.pop(-1)
        seq.pop(-1)

        if not silent:
            print("Observation ends during filter %s integration" % last_filter_name, duration)
        seq.append([duration, last_filter_name, "filter discard",
                   "Observation ends during filter %s integration" % last_filter_name])

    else:
        # end time reached during filter wheel rotation or detector readout
        last_filter_name = seqs_over[0][1]
        last_dark_name = "dk"
        et = duration

        if not silent:
            print("Observation ends between measurements", duration)
        seq.append([duration, last_filter_name, "obs end",
                   "Observation ends between measurements"])

    filter_wheel_angle = calculate_angle(last_filter_name, last_dark_name)
    rotation_time = calculate_wheel_rotation_time(filter_wheel_angle)
    int_time = filter_parameters[(last_dark_name, daynight)]["integration_time"]

    # move to final dark
    if rotation_time > 0.0:
        if not silent:
            print("Start filter wheel rotation %s->%s" % (last_filter_name, last_dark_name), et)
        seq.append([et, last_filter_name, "rot start", "Start filter wheel rotation %s->%s" %
                    (last_filter_name, last_dark_name)])
        if not silent:
            print("Filter %s now in position" % (last_dark_name), et + rotation_time)
        seq.append([et + rotation_time, last_dark_name, "rot end", "End filter wheel rotation %s->%s" %
                    (last_filter_name, last_dark_name)])

    if not silent:
        print("Filter %s frame start" % last_dark_name, et + rotation_time)
    seq.append([et, last_dark_name, "int start", "Start frame acquisition filter %s" % last_dark_name])
    if not silent:
        print("Filter %s frame end" % last_dark_name, et + rotation_time + int_time)
    seq.append([et + rotation_time + int_time, last_dark_name, "int end", "Observation end"])

    if not silent:
        print("Observation end", et + rotation_time + int_time)
    seq.append([et + rotation_time + int_time, last_dark_name, "obs end", "Observation end"])

    return seq
