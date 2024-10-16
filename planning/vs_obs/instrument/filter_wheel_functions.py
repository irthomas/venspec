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


# def calculate_sequences_duration(filter_sequence, daynight, silent=True):
#     """calculate the duration of the filter sequence plus start from dark and return to start position time"""

#     sequence = {"first dark": [], "filter": [], "last dark": []}

#     first_dark_name = "dk"
#     first_filter_name = filter_sequence[0]
#     if "dk" in first_filter_name:
#         first_filter_name = "dk"

#     # dark is always in position at observation start
#     current_filter = "dk"

#     # if first filter is not the dark though, find the integration time and the rotation time
#     if first_dark_name != first_filter_name:
#         first_dark_time = 0.0
#         sequence["first dark"].append([first_dark_time, current_filter, "int start", "Start frame acquisition filter %s" % current_filter])
#         if not silent:
#             print("First filter is not dark - dark frame start", first_dark_time)
#         first_dark_time += filter_parameters[("dk", daynight)]["integration_time"]
#         sequence["first dark"].append([first_dark_time, current_filter, "int end", "End frame acquisition filter %s" % current_filter])
#         if not silent:
#             print("First dark frame end", first_dark_time)

#         filter_wheel_angle = calculate_angle(first_dark_name, first_filter_name)
#         sequence["first dark"].append([first_dark_time, current_filter, "rot start", "Start filter wheel rotation %s->%s" %
#                                       (first_dark_name, first_filter_name)])
#         first_dark_time += calculate_wheel_rotation_time(filter_wheel_angle)
#         current_filter = first_filter_name
#         sequence["first dark"].append([first_dark_time, current_filter, "rot end", "End filter wheel rotation %s->%s" % (first_dark_name, first_filter_name)])
#         if not silent:
#             print("Rotated to first filter %s for nominal sequence start" % current_filter, first_dark_time)

#     else:  # if dark is already the first measurement in the sequence, set this to 0 and include the time in the filter sequence
#         first_dark_time = 0.0
#         sequence["first dark"].append([first_dark_time, current_filter, "", ""])
#         if not silent:
#             print("First filter in sequence is dark, time is included in filter sequence")

#     # first filter in sequence is now in position
#     seconds = 0.0
#     for i in range(len(filter_sequence)):
#         filter_name = filter_sequence[i]

#         if i < len(filter_sequence) - 1:
#             next_filter_name = filter_sequence[i + 1]
#         else:
#             next_filter_name = filter_sequence[0]  # return to first filter at the end

#         if not silent:
#             print("Filter position %s, frame start" % current_filter, seconds)
#         sequence["filter"].append([seconds, current_filter, "int start", "Start frame acquisition filter %s" % current_filter])
#         seconds += filter_parameters[(filter_name, daynight)]["integration_time"]
#         sequence["filter"].append([seconds, current_filter, "int end", "End frame acquisition filter %s" % current_filter])

#         # filter wheel rotation
#         # first change dk2 and dk4 to dk
#         if "dk" in filter_name:
#             filter_name = "dk"
#         if "dk" in next_filter_name:
#             next_filter_name = "dk"

#         if filter_name != next_filter_name:  # readout detector while filter is rotated. Assume readout << wheel movement time

#             if not silent:
#                 print("Filter position %s frame end, start rotation" % current_filter, seconds)
#             filter_wheel_angle = calculate_angle(filter_name, next_filter_name)
#             sequence["filter"].append([seconds, current_filter, "rot start", "Start filter wheel rotation %s->%s" % (filter_name, next_filter_name)])
#             seconds += calculate_wheel_rotation_time(filter_wheel_angle)
#             current_filter = next_filter_name
#             sequence["filter"].append([seconds, current_filter, "rot end", "End filter wheel rotation %s->%s" % (filter_name, next_filter_name)])

#         else:  # just read out detector if not filter change needed
#             if not silent:
#                 print("Frame end, filter already in position %s, start detector readout" % current_filter, seconds)
#             sequence["filter"].append([seconds, current_filter, "readout start", "Start detector readout filter %s" % current_filter])
#             seconds += DETECTOR_READOUT_TIME
#             sequence["filter"].append([seconds, current_filter, "readout end", "End detector readout filter %s" % current_filter])

#     if not silent:
#         print("Filter sequence end in filter position %s" % current_filter, seconds)

#     last_filter_name = filter_sequence[-1]
#     last_dark_name = "dk"
#     if "dk" in last_filter_name:
#         last_filter_name = "dk"

#     if last_filter_name != last_dark_name:  # if last filter is not the dark frame, need to include an extra rotation and integration at the end
#         last_dark_time = 0.0
#         if not silent:
#             print("Start rotation from filter %s" % last_filter_name, last_dark_time)
#         sequence["last dark"].append([last_dark_time, current_filter, "rot start", "Start filter wheel rotation %s->%s" % (last_filter_name, last_dark_name)])
#         filter_wheel_angle = calculate_angle(last_filter_name, last_dark_name)
#         last_dark_time += calculate_wheel_rotation_time(filter_wheel_angle)
#         current_filter = last_filter_name
#         sequence["last dark"].append([last_dark_time, current_filter, "rot end", "End filter wheel rotation %s->%s" % (last_filter_name, last_dark_name)])
#         if not silent:
#             print("Rotated to dark frame", last_dark_time)

#         sequence["last dark"].append([last_dark_time, current_filter, "int start", "Start frame acquisition filter %s" % current_filter])
#         if not silent:
#             print("Last filter is not dark - dark frame start", last_dark_time)
#         last_dark_time += filter_parameters[("dk", daynight)]["integration_time"]
#         sequence["last dark"].append([last_dark_time, current_filter, "int end", "End frame acquisition filter %s" % current_filter])
#         if not silent:
#             print("Last dark frame end", last_dark_time)

#     else:  # if already in position
#         last_dark_time = 0.0
#         sequence["last dark"].append([last_dark_time, current_filter, "", ""])
#         if not silent:
#             print("Last filter in sequence is dark, time is included in filter sequence")

#     if not silent:
#         print("Filter sequence:", filter_sequence)
#         print("Time taken to measure first dark and rotate to first filter position (if first filter in sequence is not dark):", first_dark_time)
#         print("Time taken to run through the filter sequence:", seconds)
#         print("Time taken to rotate to last dark position and measure the final dark (if last filter in sequence is not dark):", last_dark_time)
#     return sequence


# filter_sequence = filter_sequences["night_123d"]
# daynight = "n"

# first_dark_time, seconds, last_dark_time = calculate_sequences_duration(filter_sequence, daynight)


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
    seq.append([et, "rot end", first_filter_name, "End filter wheel rotation %s->%s" %
                (first_dark_name, first_filter_name)])

    """find complete filter sequences until observation time limit is reached"""
    et = 0.0  # relative to sza crossing time. At this time, the first filter is already in position ready to run
    # seq = []

    i = 0  # loop number (filter sequence)
    while et < duration:

        # get the next filter in the sequence. Return to first filter at the end of the sequence
        current_filter_name = filter_sequence[i % 4]
        next_filter_name = filter_sequence[(i + 1) % 4]

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
