# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 21:09:28 2022

@author: iant
"""

import numpy as np
from datetime import timedelta


from planning.vs_obs.config.constants import PRECOOLING_DURATION, STANDBY_DURATION, DETECTOR_READOUT_TIME, SAVE_HTML, CHANNEL_ID
from planning.vs_obs.spice.spice_functions import dt2et, et2dt, fov_lonlat, sc_alt, get_fov_vectors
from planning.vs_obs.spice.find_next_terminators import next_d2n_terminator, next_n2d_terminator
from planning.vs_obs.spice.find_sza_crossings import get_et_crossing_sza

from planning.vs_obs.instrument.filter_selection import obs_types, filter_parameters
from planning.vs_obs.instrument.filter_wheel_functions import calculate_angle, calculate_wheel_rotation_time
from planning.vs_obs.html.plan_table_row import plan_table_row
from planning.vs_obs.observation_cycle import observation_cycle
from planning.vs_obs.html.save_plan_table import save_plan_table


def make_orbit_plan(utc_start_time):

    # get the FOV corner vectors
    _, fov_vectors = get_fov_vectors(CHANNEL_ID)

    orbit_plan = []
    h = ""

    # start at night to day terminator
    dt_start = next_n2d_terminator(utc_start_time)

    orbit_plan.append({"dt_start": dt_start, "daynight": "d", "status": "off"})

    # find next sunset
    dt_next_sunset = next_d2n_terminator(dt_start)

    # find precooling start
    dt_standby_start = dt_next_sunset - timedelta(seconds=(PRECOOLING_DURATION + STANDBY_DURATION))
    orbit_plan.append({"dt_standby_start": dt_standby_start, "daynight": "d", "status": "standby"})

    dt_precooling_start = dt_next_sunset - timedelta(seconds=PRECOOLING_DURATION)
    orbit_plan.append({"dt_precooling_start": dt_precooling_start, "daynight": "d", "status": "precooling"})

    # start time
    et = dt2et(dt_next_sunset)

    day_ix = 0
    night_ix = 0
    for obs_ix, observation in enumerate(observation_cycle):

        daynight = observation["daynight"]

        if observation["status"] == "off":
            # no observations - increment et
            if daynight == "n":
                # find next sunrise
                dt_next_sunrise = next_n2d_terminator(et2dt(et + 60.0))
                et = dt2et(dt_next_sunrise)
                night_ix += 1

            if daynight == "d":
                # find next sunrise
                dt_next_sunset = next_d2n_terminator(et2dt(et + 60.0))
                et = dt2et(dt_next_sunset)
                day_ix += 1
            continue

        filter_cycle_name = observation["filter_cycle"]
        filter_cycle = obs_types[filter_cycle_name] * 1000

        if daynight == "n":
            # find next sunrise, add 60 seconds to ensure the et is on the nightside
            dt_next_sunrise = next_n2d_terminator(et2dt(et + 60.0))

            # loop through nightside
            i = 0
            alts = []

            while et < dt2et(dt_next_sunrise):

                et_start = float(et)

                filter_name = filter_cycle[i + night_ix]  # start with a different filter each time
                next_filter_name = filter_cycle[i + 1 + night_ix]  # start with a different filter each time
                filter_dict = filter_parameters[(filter_name, daynight)]
                integration_time = filter_dict["integration_time"]

                # start of integration time
                lonlats_s = fov_lonlat(et_start, fov_vectors)

                et += integration_time
                et_end = float(et)

                # end of integration time
                lonlats_e = fov_lonlat(et_end, fov_vectors)

                # get altitude
                alt = sc_alt(et_end)
                alts.append(alt)

                # make footprint polygon including integration time
                footprint = np.zeros((5, 2))
                footprint[0, :] = lonlats_s[1, :]
                footprint[1, :] = lonlats_e[2, :]
                footprint[2, :] = lonlats_e[3, :]
                footprint[3, :] = lonlats_s[0, :]
                footprint[4, :] = lonlats_s[1, :]

                # filter wheel rotation
                # first change dk2 and dk4 to dk
                if "dk" in filter_name:
                    filter_name = "dk"
                if "dk" in next_filter_name:
                    next_filter_name = "dk"

                # save footprint even if dark
                orbit_d = {"dt_acq_s": et2dt(et_start), "dt_acq_e": et2dt(et_end),
                           "daynight": daynight, "status": "science",
                           "semiorbit": obs_ix, "alt": alt, "filter": filter_name,
                           "footprint_s": lonlats_s, "footprint_e": lonlats_e,
                           "footprint": footprint, "day_ix": day_ix, "night_ix": night_ix, "spec_ix": i}

                if filter_name != next_filter_name:  # readout detector while filter is rotated

                    filter_wheel_angle = calculate_angle(filter_name, next_filter_name)

                    orbit_d["dt_rot_s"] = et2dt(et)
                    rotation_time = calculate_wheel_rotation_time(filter_wheel_angle)
                    et += rotation_time
                    orbit_d["filter_transition"] = {"start": filter_name, "end": next_filter_name,
                                                    "angle": filter_wheel_angle, "time": rotation_time}
                    orbit_d["dt_rot_e"] = et2dt(et)

                else:  # just read out detector

                    orbit_d["dt_read_s"] = et2dt(et)
                    et += DETECTOR_READOUT_TIME
                    orbit_d["filter_transition"] = {"angle": 0.0}
                    orbit_d["dt_read_e"] = et2dt(et)

                i += 1
                orbit_plan.append(orbit_d)
            # print(day_ix, night_ix, daynight, np.min(alts), np.mean(alts), np.max(alts))

            # pop last entry if over the sunrise
            if et > dt2et(dt_next_sunrise):
                orbit_plan.pop(-1)

            night_ix += 1

        elif daynight == "d":

            """repeat for dayside"""
            # find next sunset, add 60 seconds to ensure the et is on the dayside
            dt_next_sunset = next_d2n_terminator(et2dt(et + 60.))

            # loop through dayside
            i = 0
            while et < dt2et(dt_next_sunset):

                et_start = float(et)

                filter_name = filter_cycle[i + day_ix]  # start with a different filter each time
                next_filter_name = filter_cycle[i + 1 + day_ix]  # start with a different filter each time
                filter_dict = filter_parameters[(filter_name, daynight)]
                integration_time = filter_dict["integration_time"]

                # start of integration time
                lonlats_s = fov_lonlat(et_start, fov_vectors)

                et += integration_time
                et_end = float(et)

                # end of integration time
                lonlats_e = fov_lonlat(et_end, fov_vectors)

                # get altitude
                alt = sc_alt(et_end)
                alts.append(alt)

                # make footprint polygon including integration time
                footprint = np.zeros((5, 2))
                footprint[0, :] = lonlats_s[1, :]
                footprint[1, :] = lonlats_e[2, :]
                footprint[2, :] = lonlats_e[3, :]
                footprint[3, :] = lonlats_s[0, :]
                footprint[4, :] = lonlats_s[1, :]

                # filter wheel rotation
                # first change dk2 and dk4 to dk
                if "dk" in filter_name:
                    filter_name = "dk"
                if "dk" in next_filter_name:
                    next_filter_name = "dk"

                # save footprint even if dark
                orbit_d = {"et_acq_s": et_start, "et_acq_e": et_end,
                           "daynight": daynight, "status": "science",
                           "semiorbit": obs_ix, "alt": alt, "filter": filter_name,
                           "footprint_s": lonlats_s, "footprint_e": lonlats_e,
                           "footprint": footprint, "day_ix": day_ix, "night_ix": night_ix, "spec_ix": i}

                if filter_name != next_filter_name:

                    filter_wheel_angle = calculate_angle(filter_name, next_filter_name)

                    orbit_d["dt_rot_s"] = et2dt(et)
                    rotation_time = calculate_wheel_rotation_time(filter_wheel_angle)
                    et += rotation_time
                    orbit_d["filter_transition"] = {"start": filter_name, "end": next_filter_name,
                                                    "angle": filter_wheel_angle, "time": rotation_time}
                    orbit_d["dt_rot_e"] = et2dt(et)

                else:

                    orbit_d["dt_read_s"] = et2dt(et)
                    et += DETECTOR_READOUT_TIME
                    orbit_d["filter_transition"] = {"angle": 0.0}
                    orbit_d["dt_read_e"] = et2dt(et)

                i += 1
                orbit_plan.append(orbit_d)
            # print(day_ix, night_ix, daynight, np.min(alts), np.mean(alts), np.max(alts))

            # pop last entry if over the sunrise
            if et > dt2et(dt_next_sunset):
                orbit_plan.pop(-1)

            day_ix += 1

        # add 20 seconds of standby
        et += 20.

        orbit_d = {"daynight": "", "status": "standby"}

        orbit_plan.append(orbit_d)

    return orbit_plan
