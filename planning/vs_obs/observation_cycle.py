# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 15:37:39 2022

@author: iant
"""


"""typical orbit: 4 nightside and 4 daysides"""
observation_name = "Typical orbit: nightside filters: 1, 2, 3, dark; dayside filters: 2, 4, dark"
observation_cycle = []
# #add 4 orbits on
observation_cycle.extend([
    {"daynight": "n", "status": "on", "filter_cycle": "night_123d"},
    {"daynight": "d", "status": "on", "filter_cycle": "day_22h4d"},
] * 4)
# add 11 orbits off
observation_cycle.extend([
    {"daynight": "n", "status": "off"},
    {"daynight": "d", "status": "off"},
] * 11)

# repeat cycle N times
observation_cycle *= 2  # N cycles


"""scenario 1: 4 nightsides only; filter 1, 2, 3 dark"""
# observation_name = "Nightside only: nightside filters: 1, 2, 3, dark"
# observation_cycle = []
# # add 4 orbits on
# observation_cycle.extend([
#     {"daynight": "n", "status": "on", "filter_cycle": "night_123d"},
#     {"daynight": "d", "status": "off"},
# ] * 4)
# # add 11 orbits off
# observation_cycle.extend([
#     {"daynight": "n", "status": "off"},
#     {"daynight": "d", "status": "off"},
# ] * 11)

# # repeat cycle N times
# observation_cycle *= 5  # N cycles


"""scenario 2: 4 nightsides only; filter 1, 1, 1, dark on orbit 2; rest are 1, 2, 3 dark"""
# observation_name = "Nightside only: filter 1, 1, 1, dark on orbit 2; filters 1, 2, 3, dark on orbits 1, 3, 4"
# observation_cycle = []
# #add 4 orbits on
# observation_cycle.extend([
#     {"daynight":"n", "status":"on", "filter_cycle":"night_123d"},
#     {"daynight":"d", "status":"off"},
#     ] * 1)
# observation_cycle.extend([
#     {"daynight":"n", "status":"on", "filter_cycle":"night_111d"},
#     {"daynight":"d", "status":"off"},
#     ] * 1)
# observation_cycle.extend([
#     {"daynight":"n", "status":"on", "filter_cycle":"night_123d"},
#     {"daynight":"d", "status":"off"},
#     ] * 2)
# #add 11 orbits off
# observation_cycle.extend([
#     {"daynight":"n", "status":"off"},
#     {"daynight":"d", "status":"off"},
#     ] * 11)

# #repeat cycle N times
# observation_cycle *= 28 #N cycles


"""scenario 3: 4 nightsides only; filter 1, 1, 1, dark on orbits 2 and 4; rest are 1, 2, 3 dark"""
# observation_name = "Nightside only: filter 1, 1, 1, dark on orbits 2 and 4; filters 1, 2, 3, dark on orbits 1 and 3"
# observation_cycle = []
# #add 4 orbits on
# observation_cycle.extend([
#     {"daynight":"n", "status":"on", "filter_cycle":"night_123d"},
#     {"daynight":"d", "status":"off"},
#     ] * 1)
# observation_cycle.extend([
#     {"daynight":"n", "status":"on", "filter_cycle":"night_111d"},
#     {"daynight":"d", "status":"off"},
#     ] * 1)
# observation_cycle.extend([
#     {"daynight":"n", "status":"on", "filter_cycle":"night_123d"},
#     {"daynight":"d", "status":"off"},
#     ] * 1)
# observation_cycle.extend([
#     {"daynight":"n", "status":"on", "filter_cycle":"night_111d"},
#     {"daynight":"d", "status":"off"},
#     ] * 1)
# #add 11 orbits off
# observation_cycle.extend([
#     {"daynight":"n", "status":"off"},
#     {"daynight":"d", "status":"off"},
#     ] * 11)

# #repeat cycle N times
# observation_cycle *= 28 #N cycles
