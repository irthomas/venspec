# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 09:45:05 2022

@author: iant

VENSPEC-H OPERATIONS PLAN

HIGH DATA RATE:
    DAYSIDE:
        MINIMAL ACROSS-TRACK BINNING
        SHORT INTEGRATION TIME - MAY NEED MULTIPLE ACCUMULATIONS?
        DARK FRAMES EVERY N MEASUREMENTS
    NIGHTSIDE:
        BINNED ACROSS-TRACK
        LONG INTEGRATION TIME - SINGLE ACCUMULATION
        DARK FRAMES EVERY N MEASUREMENTS

"""
from scipy.spatial import ConvexHull
import os
# import spiceypy as sp
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from planning.vs_obs.spice.find_lat_crossings import get_et_crossing_lats
from planning.vs_obs.spice.find_sza_crossings import get_et_crossing_szas, get_corners_sza
# from tools.general.progress_bar import progress
from planning.vs_obs.spice.spice_functions import et2dt, get_centre_lonlat, get_corners_lonlats
from planning.vs_obs.spice.get_orbit_number import get_orbit_number, read_orbit_numbers
from planning.vs_obs.instrument.filter_selection import filter_sequences
from planning.vs_obs.instrument.filter_wheel_functions import make_sequence

from planning.vs_obs.config.constants import PRECOOLING_DURATION, STANDBY_DURATION, CCU_RECONFIGURATION_TIME

from planning.vs_obs.io.get_venus_magellan_data import get_venus_magellan_data


# save to a log file?
# OUTPUT_LOG = True
OUTPUT_LOG = False


orbit_number_dts = read_orbit_numbers()

if os.path.exists("event_log.txt"):
    print("Deleting existing log file")
    os.remove("event_log.txt")


def write_event(etdt, comment):
    with open("event_log.txt", "a") as f:
        if isinstance(etdt, datetime):
            line = "%s\t%s" % (datetime.strftime(etdt, "%Y-%m-%dT%H:%M:%S.%f")[:-4], comment)
        elif isinstance(etdt, float):
            line = "%s\t%s" % (datetime.strftime(et2dt(etdt), "%Y-%m-%dT%H:%M:%S.%f")[:-4], comment)
        elif etdt == "":
            line = "%s" % comment
        f.write(line + "\n")


# must start with a nightside in this version of the code
utc_start_time = datetime(2035, 6, 1)

nightside_sza_limit = 100.0  # for all points in FOV
dayside_sza_limit = 95.0  # for centre of FOV


# plots = ["sza night", "sza day", "lat night", "lat day", "groundtracks"]  # plot everything
# plots = ["lat night", "lat day", "groundtracks"]  # plot latitudes and groundtracks (recommended)
# plots = ["lat night", "lat day"]  # plot latitudes
plots = ["groundtracks"]  # plot just groundtracks
# plots = []  # plot nothing

n_blocks = 5  # number of 4-orbit blocks to be simulated (more = longer calculation)
# n_blocks = 1  # number of 4-orbit blocks to be simulated

n_orbits = 4  # number of orbits within each venspec block (normal operations)
# n_orbits = 1  # number of orbits within each venspec block (for testing only)

nightside_filter_sequence_name = "night_123d"
dayside_filter_sequence_name = "day_22h44hd"

"""start the calculations"""
# get filter sequences, 1 for each of the 4 orbits
nightside_filter_sequences = filter_sequences[nightside_filter_sequence_name]
dayside_filter_sequences = filter_sequences[dayside_filter_sequence_name]


# calculate the orbit block start and end times
start_orbit_number = get_orbit_number(orbit_number_dts, utc_start_time)
start_orbit_numbers = [start_orbit_number + i for i in np.arange(n_blocks) * 15]
end_orbit_numbers = [i + n_orbits for i in start_orbit_numbers]


print("Planning VenSpec orbits %i to %i" % (start_orbit_numbers[0], end_orbit_numbers[-1]))

v_block_dt_starts = [k[0] - timedelta(seconds=10) for i, k in enumerate(orbit_number_dts.keys()) if i in start_orbit_numbers]
v_block_dt_ends = [k[0] + timedelta(minutes=10) for i, k in enumerate(orbit_number_dts.keys()) if i in end_orbit_numbers]

# v_block_dt_starts = [utc_start_time + timedelta(seconds=int(i)) for i in np.arange(n_blocks) * 5620 * 15]
# v_block_dt_ends = [v_block_dt_start + timedelta(minutes=390) for v_block_dt_start in v_block_dt_starts]

utc_start_str = str(v_block_dt_starts[0])
utc_end_str = str(v_block_dt_ends[-1])

orbit_plan = []
for v_block_dt_start, v_block_dt_end in zip(v_block_dt_starts, v_block_dt_ends):

    # calculate the theoretical nightside start and end times for the 4 orbits, considering only sza/lat constraints
    print("Calculating nightside timings")
    plot = "sza night" in plots
    sza_crossings_night = get_et_crossing_szas(v_block_dt_start, v_block_dt_end, nightside_sza_limit, plot=plot)
    if sza_crossings_night["type"] != "Nightside first":
        print("Error: at present, code is only tested for nightside first")

    # check szas
    # for et in sza_crossings["ascending"].keys():
    #     print(sza_crossings["ascending"][et]["corner szas"])

    # calculate the start and end latitudes where the SZA constraint is met for all 4 corners on all 4 orbits
    lat_nightside_start = np.max([sza_crossings_night["ascending"][et]["centre lat"] for et in sza_crossings_night["ascending"].keys()])
    lat_nightside_end = np.min([sza_crossings_night["descending"][et]["centre lat"] for et in sza_crossings_night["descending"].keys()])

    # once the real nightside start and end latitudes are known, we can calculate all start/end times for the 4 orbits
    plot = "lat night" in plots
    lat_nightside_crossings = get_et_crossing_lats(v_block_dt_start, v_block_dt_end, [lat_nightside_start, lat_nightside_end], plot=plot)
    if lat_nightside_crossings["type"] != "Latitudes descending":
        print("Error: at present, code is only tested for descending nightside latitudes")

    # check szas
    # for et in lat_crossings["descending"].keys():
    #     print(lat_crossings["descending"][et]["corner szas"])

    # get descending start times. These are the start times for the first Venus obs (not the dark)
    ets_nightside_start = [et for et in list(lat_nightside_crossings["descending"].keys())[::2]]
    # sza crossing times near end of observation. The real obs end time will depend on the final dark
    # this is the time the TC is sent to end the observation
    ets_nightside_end = [et for et in list(lat_nightside_crossings["descending"].keys())[1::2]]

    dts_nightside_start = [et2dt(et) for et in ets_nightside_start]
    dts_nightside_end = [et2dt(et) for et in ets_nightside_end]

    if "lat night" in plots:
        for i, (et, dt) in enumerate(zip(ets_nightside_start, dts_nightside_start)):
            plt.text(et, lat_nightside_start, "Start %i %s" % (i + 1, str(dt)[:-7]))
        for i, (et, dt) in enumerate(zip(ets_nightside_end, dts_nightside_end)):
            plt.text(et, lat_nightside_end, "End %i %s" % (i + 1, str(dt)[:-7]))

    # calculate the theoretical dayside start and end times for the 4 orbits
    print("Calculating dayside timings")
    plot = "sza day" in plots
    sza_crossings_day = get_et_crossing_szas(v_block_dt_start, v_block_dt_end, dayside_sza_limit, limit_type="mean", plot=plot)
    if sza_crossings_day["type"] != "Nightside first":
        print("Error: at present, code is only tested for nightside first")

    # check szas
    # for et in sza_crossings_day["descending"].keys():
    #     print(sza_crossings_day["descending"][et]["corner szas"])

    # this time we use descending sza for observation start and ascending sza for observation end
    # calculate the start and end latitudes where the SZA constraint is met on all 4 orbits
    lat_dayside_start = np.max([sza_crossings_day["descending"][et]["centre lat"] for et in sza_crossings_day["descending"].keys()])
    lat_dayside_end = np.min([sza_crossings_day["ascending"][et]["centre lat"] for et in sza_crossings_day["ascending"].keys()])

    # once the real nightside start and end latitudes are known, we can calculate all start/end times for the 4 orbits
    # the first dayside goes from the 2nd descending latitude to the 3rd descending latitude
    plot = "lat day" in plots
    lat_dayside_crossings = get_et_crossing_lats(v_block_dt_start, v_block_dt_end, [lat_dayside_start, lat_dayside_end], plot=plot)
    if lat_dayside_crossings["type"] != "Latitudes descending":
        print("Error: at present, code is only tested for descending nightside latitudes")

    # check szas
    # for et in lat_dayside_crossings["descending"].keys():
    #     print(lat_dayside_crossings["descending"][et]["corner szas"])

    # get descending start times. These are the start times for the first dayside Venus obs (not the dark)
    ets_dayside_start = [et for et in list(lat_dayside_crossings["descending"].keys())[1:-1:2]]
    # sza crossing times near end of observation. The real obs end time will depend on filter measurements
    ets_dayside_end = [et for et in list(lat_dayside_crossings["descending"].keys())[2::2]]

    dts_dayside_start = [et2dt(et) for et in ets_dayside_start]
    dts_dayside_end = [et2dt(et) for et in ets_dayside_end]

    if "lat day" in plots:
        for i, (et, dt) in enumerate(zip(ets_dayside_start, dts_dayside_start)):
            plt.text(et, lat_dayside_start, "Start %i %s" % (i + 1, str(dt)[:-7]))
        for i, (et, dt) in enumerate(zip(ets_dayside_end, dts_dayside_end)):
            plt.text(et, lat_dayside_end, "End %i %s" % (i + 1, str(dt)[:-7]))

    # now figure out how many filters to run in this time
    for orbit_ix in range(n_orbits):

        # orbit_number = 0
        nightside_filter_sequence = nightside_filter_sequences[orbit_ix]
        dayside_filter_sequence = dayside_filter_sequences[orbit_ix]
        print("Running orbit %i in block, filter sequences:" % (orbit_ix + 1), nightside_filter_sequence, dayside_filter_sequence)

        # nightside
        daynight = "n"

        et_nightside_start = ets_nightside_start[orbit_ix]
        et_nightside_end = ets_nightside_end[orbit_ix]
        dt_nightside_start = et2dt(et_nightside_start)
        dt_nightside_end = et2dt(et_nightside_end)
        duration = et_nightside_end - et_nightside_start

        lon_nightside_start, lat_nightside_start = get_centre_lonlat(et_nightside_start)
        lon_nightside_end, lat_nightside_end = get_centre_lonlat(et_nightside_end)

        sequence_nightside_relative = make_sequence(nightside_filter_sequence, daynight, duration, silent=True)

        sequence_nightside_et = [{"et": seq[0] + et_nightside_start, "filter": seq[1],
                                  "operation": seq[2], "comment": seq[3]} for seq in sequence_nightside_relative]

        # calculate geometry for each sequence, add to each sequence
        ets = [seq["et"] for seq in sequence_nightside_et]
        lonlats = get_corners_lonlats(ets)
        szas = get_corners_sza(ets)
        for i, sequence in enumerate(sequence_nightside_et):
            if sequence["operation"] in ["int start", "int end"]:
                sequence["corner lonlats"] = lonlats[i, :, :]
                sequence["corner szas"] = szas[i, :]
            sequence["dt"] = et2dt(sequence["et"])

        # dayside
        daynight = "d"

        et_dayside_start = ets_dayside_start[orbit_ix]
        et_dayside_end = ets_dayside_end[orbit_ix]
        dt_dayside_start = et2dt(et_dayside_start)
        dt_dayside_end = et2dt(et_dayside_end)
        duration = et_dayside_end - et_dayside_start

        lon_dayside_start, lat_dayside_start = get_centre_lonlat(et_dayside_start)
        lon_dayside_end, lat_dayside_end = get_centre_lonlat(et_dayside_end)

        sequence_dayside_relative = make_sequence(dayside_filter_sequence, daynight, duration, silent=True)

        sequence_dayside_et = [{"et": seq[0] + et_dayside_start, "filter": seq[1], "operation": seq[2], "comment": seq[3]} for seq in sequence_dayside_relative]

        # timings are now calculated
        # calculate geometry for each sequence, add to each sequence
        ets = [seq["et"] for seq in sequence_dayside_et]
        lonlats = get_corners_lonlats(ets)
        szas = get_corners_sza(ets)
        for i, sequence in enumerate(sequence_dayside_et):
            if sequence["operation"] in ["int start", "int end"]:
                sequence["corner lonlats"] = lonlats[i, :, :]
                sequence["corner szas"] = szas[i, :]
            sequence["dt"] = et2dt(sequence["et"])

        # add precooling and standby time before start of first half orbit
        if orbit_ix == 0:
            # get standby start time from first nightside start time

            et_block_start = np.min((sequence_nightside_et[0]["et"], sequence_dayside_et[0]["et"]))
            dt_block_start = np.min((sequence_nightside_et[0]["dt"], sequence_dayside_et[0]["dt"]))

            et_standby_start = et_block_start - PRECOOLING_DURATION - STANDBY_DURATION
            dt_standby_start = dt_block_start - timedelta(seconds=(PRECOOLING_DURATION + STANDBY_DURATION))

            # get precooling start
            et_precooling_start = et_block_start - PRECOOLING_DURATION
            dt_precooling_start = dt_block_start - timedelta(seconds=(PRECOOLING_DURATION))

            orbit_number = get_orbit_number(orbit_number_dts, dt_standby_start)
            orbit = {"orbit number": orbit_number, "event": "Standby start", "daynight": "", "start time": {"et": et_standby_start, "dt": dt_standby_start}}
            orbit_plan.append(orbit)

            orbit_number = get_orbit_number(orbit_number_dts, dt_precooling_start)
            orbit = {"orbit number": orbit_number, "event": "Standby end", "daynight": "", "start time": {"et": et_precooling_start, "dt": dt_precooling_start}}
            orbit_plan.append(orbit)

            orbit = {"orbit number": orbit_number, "event": "Precooling start", "daynight": "",
                     "start time": {"et": et_precooling_start, "dt": dt_precooling_start}
                     }
            orbit_plan.append(orbit)

            orbit_number = get_orbit_number(orbit_number_dts, dt_block_start)
            orbit = {"orbit number": orbit_number, "event": "Precooling end", "daynight": "",
                     "start time": {"et": et_block_start, "dt": dt_block_start}
                     }
            orbit_plan.append(orbit)

            # record telecommands
            orbit_number = get_orbit_number(orbit_number_dts, dt_standby_start)
            orbit = {"orbit number": orbit_number, "event": "Send TC start", "daynight": "", "start time": {"et": et_standby_start, "dt": dt_standby_start}}
            orbit_plan.append(orbit)
        else:
            # record telecommands
            et_nightside_tc_start = sequence_nightside_et[0]["et"] - CCU_RECONFIGURATION_TIME
            dt_nightside_tc_start = sequence_nightside_et[0]["dt"] - timedelta(seconds=CCU_RECONFIGURATION_TIME)

            orbit_number = get_orbit_number(orbit_number_dts, dt_nightside_tc_start)
            orbit = {"orbit number": orbit_number, "event": "Send TC start", "daynight": "n",
                     "start time": {"et": et_nightside_tc_start, "dt": dt_nightside_tc_start}}
            orbit_plan.append(orbit)

            et_dayside_tc_start = sequence_dayside_et[0]["et"] - CCU_RECONFIGURATION_TIME
            dt_dayside_tc_start = sequence_dayside_et[0]["dt"] - timedelta(seconds=CCU_RECONFIGURATION_TIME)

            orbit_number = get_orbit_number(orbit_number_dts, dt_dayside_tc_start)
            orbit = {"orbit number": orbit_number, "event": "Send TC start", "daynight": "d",
                     "start time": {"et": et_dayside_tc_start, "dt": dt_dayside_tc_start}}
            orbit_plan.append(orbit)

        orbit_number = get_orbit_number(orbit_number_dts, dt_nightside_end)
        orbit = {"orbit number": orbit_number, "event": "Send TC end", "daynight": "n", "start time": {"et": et_nightside_end, "dt": dt_nightside_end}}
        orbit_plan.append(orbit)

        orbit_number = get_orbit_number(orbit_number_dts, dt_dayside_end)
        orbit = {"orbit number": orbit_number, "event": "Send TC end", "daynight": "d", "start time": {"et": et_dayside_end, "dt": dt_dayside_end}}
        orbit_plan.append(orbit)

        # collect nightside info
        orbit_number = get_orbit_number(orbit_number_dts, sequence_nightside_et[0]["dt"])
        orbit_nightside = {"orbit number": orbit_number, "obs type": "nadir", "daynight": "n",
                           "sza start time": {"et": et_nightside_start, "dt": dt_nightside_start},
                           "sza end time": {"et": et_nightside_end, "dt": dt_nightside_end},
                           "sza start geom": {"lon": lon_nightside_start, "lat": lat_nightside_start, "sza": nightside_sza_limit},
                           "sza end geom": {"lon": lon_nightside_end, "lat": lat_nightside_end, "sza": nightside_sza_limit},

                           "obs start time": {"et": sequence_nightside_et[0]["et"], "dt": sequence_nightside_et[0]["dt"]},
                           "obs end time": {"et": sequence_nightside_et[-1]["et"], "dt": sequence_nightside_et[-1]["dt"]},
                           # "obs start geom": {"lon": lon_nightside_obs_start, "lat": lat_nightside_obs_start},
                           # "obs end geom": {"lon": lon_nightside_obs_end, "lat": lat_nightside_obs_end},

                           "filter_sequence": sequence_nightside_et,
                           }

        # collect dayside info
        orbit_number = get_orbit_number(orbit_number_dts, sequence_dayside_et[0]["dt"])
        orbit_dayside = {"orbit number": orbit_number, "obs type": "nadir", "daynight": "d",
                         "sza start time": {"et": et_dayside_start, "dt": dt_dayside_start},
                         "sza end time": {"et": et_dayside_end, "dt": dt_dayside_end},
                         "sza start geom": {"lon": lon_dayside_start, "lat": lat_dayside_start, "sza": dayside_sza_limit},
                         "sza end geom": {"lon": lon_dayside_end, "lat": lat_dayside_end, "sza": dayside_sza_limit},

                         "obs start time": {"et": sequence_dayside_et[0]["et"], "dt": sequence_dayside_et[0]["dt"]},
                         "obs end time": {"et": sequence_dayside_et[-1]["et"], "dt": sequence_dayside_et[-1]["dt"]},
                         # "obs start geom": {"lon": lon_dayside_obs_start, "lat": lat_dayside_obs_start},
                         # "obs end geom": {"lon": lon_dayside_obs_end, "lat": lat_dayside_obs_end},

                         "filter_sequence": sequence_dayside_et,
                         }

        orbit_number = get_orbit_number(orbit_number_dts, dt_nightside_start)
        orbit_nightside_start = {"orbit number": orbit_number, "daynight": "n", "event": "Nightside start SZA >= %i degrees" %
                                 nightside_sza_limit, "start time": {"et": et_nightside_start, "dt": dt_nightside_start}}

        orbit_number = get_orbit_number(orbit_number_dts, dt_nightside_end)
        orbit_nightside_end = {"orbit number": orbit_number, "daynight": "n", "event": "Nightside end SZA <= %i degrees" %
                               nightside_sza_limit, "start time": {"et": et_nightside_end, "dt": dt_nightside_end}}

        orbit_number = get_orbit_number(orbit_number_dts, dt_dayside_start)
        orbit_dayside_start = {"orbit number": orbit_number, "daynight": "n", "event": "Dayside start SZA <= %i degrees" %
                               dayside_sza_limit, "start time": {"et": et_dayside_start, "dt": dt_dayside_start}}

        orbit_number = get_orbit_number(orbit_number_dts, dt_dayside_end)
        orbit_dayside_end = {"orbit number": orbit_number, "daynight": "n", "event": "Dayside end SZA >= %i degrees" %
                             dayside_sza_limit, "start time": {"et": et_dayside_end, "dt": dt_dayside_end}}

        orbit_plan.append(orbit_nightside_start)
        orbit_plan.append(orbit_nightside_end)
        orbit_plan.append(orbit_dayside_start)
        orbit_plan.append(orbit_dayside_end)

        if sequence_nightside_et[0]["et"] < sequence_dayside_et[0]["et"]:
            # if nightside comes before dayside

            # save nightside to orbit plan
            orbit_plan.append(orbit_nightside)
            # save dayside to orbit plan
            orbit_plan.append(orbit_dayside)
        else:
            # if nightside comes after dayside

            # save dayside to orbit plan
            orbit_plan.append(orbit_dayside)
            # save nightside to orbit plan

            orbit_plan.append(orbit_nightside)


# sort by orbit plan timing
ets_all = np.zeros(len(orbit_plan)) - 999
for i, orbit in enumerate(orbit_plan):
    if "start time" in orbit.keys():
        ets_all[i] = orbit["start time"]["et"]
    if "obs start time" in orbit.keys():
        ets_all[i] = orbit["obs start time"]["et"]


# sort orbits into chronological order
sort_ixs = np.argsort(ets_all)
orbit_plan = [orbit_plan[sort_ix] for sort_ix in sort_ixs]

# write orbit plan to file
if OUTPUT_LOG:
    print("Writing log file")
    comments = []
    headers = "\t".join(["Datetime", "Orbit Number", "Day/night", "Action", "Lon1", "Lat1", "SZA1",
                        "Lon2", "Lat2", "SZA2", "Lon3", "Lat3", "SZA3", "Lon4", "Lat4", "SZA4"])
    write_event("", headers)
    for orbit in orbit_plan:
        daynight = {"d": "Day", "n": "Night", "": "-"}[orbit["daynight"]]
        if "event" in orbit.keys():
            comment1 = "\t".join(["%04i" % orbit["orbit number"], daynight, orbit["event"]])
        else:
            comment1 = "\t".join(["%04i" % orbit["orbit number"], daynight])

        if "filter_sequence" in orbit.keys():
            for i, sequence in enumerate(orbit["filter_sequence"]):
                if "corner lonlats" in sequence.keys():
                    comment2 = sequence["comment"] + "\t" + "\t".join(["%0.3f\t%0.3f\t%0.3f" % v for v in zip(sequence["corner lonlats"]
                                                                      [:, 0], sequence["corner lonlats"][:, 1], sequence["corner szas"][:])])
                else:
                    comment2 = sequence["comment"]
                comments.append([sequence["dt"], comment1 + "\t" + comment2])
        else:
            comments.append([orbit["start time"]["dt"], comment1])

    # sort comments into chronological order
    sort_ixs = np.argsort([comment[0] for comment in comments])
    comments = [comments[sort_ix] for sort_ix in sort_ixs]

    # write to event file
    for comment in comments:
        write_event(comment[0], comment[1])


if "groundtracks" in plots:
    print("Plotting results")
    filter_colours = {"1": "blueviolet", "2": "firebrick", "2h": "tomato", "3": "limegreen", "4": "blue", "4h": "lightblue", "dk": "black"}

    cmap = LinearSegmentedColormap.from_list("", colors=["black", "darkorange", "lemonchiffon"], N=32)

    # plot venus map
    venus_topo = get_venus_magellan_data()
    plt.figure(figsize=(10, 6), constrained_layout=False)
    plt.imshow(venus_topo, extent=(-180, 180, -90, 90), interpolation="bilinear", cmap=cmap, aspect=1, alpha=0.5)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title("VenSpec-H Observation Plan for Orbits %i-%i\n%s to %s" % (start_orbit_numbers[0], end_orbit_numbers[-1], utc_start_str, utc_end_str))

    for orbit in orbit_plan:
        tracks = {}

        if "filter_sequence" in orbit:
            daynight = orbit["daynight"]
            for i in range(len(orbit["filter_sequence"])):
                if orbit["filter_sequence"][i]["operation"] == "int start":
                    points1 = orbit["filter_sequence"][i]["corner lonlats"]
                    points2 = orbit["filter_sequence"][i+1]["corner lonlats"]

                    points = np.concatenate((points1, points2))

                    # make footprint polygon including integration time
                    hull = ConvexHull(points)

                    vertices = np.concatenate((hull.vertices, hull.vertices[0:1]))

                    filter_name = orbit["filter_sequence"][i]["filter"]
                    if filter_name not in tracks.keys():
                        tracks[filter_name] = []
                    tracks[filter_name].append(points[vertices, :])

        for filter_name in tracks.keys():
            colour = filter_colours[filter_name]
            for polygon in tracks[filter_name]:
                if np.std(polygon[:, 0]) < 20.0:
                    plt.plot(polygon[:, 0], polygon[:, 1], color=colour)

        if "sza start geom" in orbit:
            plt.scatter(orbit["sza start geom"]["lon"], orbit["sza start geom"]["lat"], color="k")
        if "sza end geom" in orbit:
            plt.scatter(orbit["sza end geom"]["lon"], orbit["sza end geom"]["lat"], color="k")


# print("Ops start time", utc_start_time)

# print("Making orbit plan")
# start = time.time()
# orbit_plan = make_orbit_plan(utc_start_time, fov_vectors)
# end = time.time()
# print("Orbit plan made in %0.1fs" % (end - start))

# print("Plotting footprints")
# start = time.time()
# plot_footprint(orbit_plan)
# end = time.time()
# print("Footprints plotted in %0.1fs" % (end - start))
