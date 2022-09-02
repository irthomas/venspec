# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 14:56:03 2022

@author: iant

CALCULATE FILTER WHEEL MOVEMENTS


"""

import numpy as np
import itertools

POLARISATION = True
# POLARISATION = False

if POLARISATION:
    FILTER_NAMES = ["D", "F1", "F2", "F3", "F4", "F2P"] #DARK, NIGHT, DAY+NIGHT, NIGHT, DAY, DAY-POLARISED
else:
    FILTER_NAMES = ["D", "F1", "F2", "F3", "F4", "O"] #DARK, NIGHT, DAY+NIGHT, NIGHT, DAY, OPEN


SAVE_OUTPUT = True
# SAVE_OUTPUT = False


# VERBOSE = True #warning, slow!
VERBOSE = False

ANGLE_BETWEEN_FILTERS = 60.0
WHEEL_ROTATION_SPEED = 80.0 #degrees per second

OBSERVATION_TIME = 2223.3 #1800.0 #seconds for one day or night observation
MISSION_DURATION = 4. * 365. #days of mission
# OBSERVING_ORBITS_PER_DAY = 4. #orbits on which observations are made per day i.e. 4 = 4 daysides and 4 nightsides


filter_cycles = {
    "Nightside H2O low altitude":{"cycle":["F1"], "readout_time":14.4, "orbits_per_day":5.},
    "Nightside SO2 + CO":{"cycle":["F2"], "readout_time":14.4, "orbits_per_day":5.},
    "Nightside Filter 3 only":{"cycle":["F3"], "readout_time":14.4, "orbits_per_day":5.},
    "Nightside H2O vertical profile":{"cycle":["F1", "F3"], "readout_time":14.4, "orbits_per_day":5.},
    "Nightside H2O and SO2 anticorrelation":{"cycle":["F3", "F2"], "readout_time":14.4, "orbits_per_day":5.},
    "Nightside Complete set":{"cycle":["F1", "F3", "F2"], "readout_time":14.4, "orbits_per_day":5.},

    "Dayside SO2 + CO":{"cycle":["F2"], "readout_time":1.9, "orbits_per_day":4.},
    "Dayside H2O and SO2 anticorrelation":{"cycle":["F2", "F4"], "readout_time":1.9, "orbits_per_day":4.},
    "Dayside Filter 4 only":{"cycle":["F4"], "readout_time":1.9, "orbits_per_day":4.},

    "Nightside H2O low altitude alternate darks":{"cycle":["D", "F1"], "readout_time":14.4, "orbits_per_day":5.},
    "Nightside SO2 + CO alternate darks":{"cycle":["D", "F2"], "readout_time":14.4, "orbits_per_day":5.},
    "Nightside Filter 3 alternate darks":{"cycle":["D", "F3"], "readout_time":14.4, "orbits_per_day":5.},
    "Nightside H2O vertical profile alternate darks":{"cycle":["D", "F1", "D", "F3"], "readout_time":14.4, "orbits_per_day":5.},
    "Nightside H2O and SO2 anticorrelation alternate darks":{"cycle":["D", "F3", "D", "F2"], "readout_time":14.4, "orbits_per_day":5.},
    "Nightside Complete set alternate darks":{"cycle":["D", "F1", "D", "F3", "D", "F2"], "readout_time":14.4, "orbits_per_day":5.},

    "Dayside SO2 + CO alternate darks":{"cycle":["D", "F2"], "readout_time":1.9, "orbits_per_day":4.},
    "Dayside H2O and SO2 anticorrelation alternate darks":{"cycle":["D", "F2", "D", "F4"], "readout_time":1.9, "orbits_per_day":4.},
    "Dayside Filter 4 alternate darks":{"cycle":["D", "F4"], "readout_time":1.9, "orbits_per_day":4.},

}

if POLARISATION:
    filter_cycles["Dayside SO2 + CO and polarisation"] = {"cycle":["F2", "F2P"], "readout_time":1.9, "orbits_per_day":4.}
    filter_cycles["Dayside H2O + SO2 anticorrelation and polarisation"] = {"cycle":["F2", "F2P", "F4"], "readout_time":1.9, "orbits_per_day":4.}
    filter_cycles["Dayside SO2 + CO alternate darks and polarisation"] = {"cycle":["D", "F2", "F2P"], "readout_time":1.9, "orbits_per_day":4.}
    filter_cycles["Dayside H2O + SO2 anticorrelation alternate darks and polarisation"] = {"cycle":["D", "F2", "F2P", "D", "F4"], "readout_time":1.9, "orbits_per_day":4.}




def wheel_filter_permutations(filter_positions):
    
    #no need to permutate the first filter position (Dark)
    permutations = list(itertools.permutations(filter_positions[1:], 5))
    
    #add dark back into permutations
    positions = [["D"] + list(p) for p in permutations]
    return positions






def calculate_angle(filter_positions, pos1, pos2):
    """calculate angle from start and end filter position (degrees)"""
    pos1_i = filter_positions.index(pos1)
    pos2_i = filter_positions.index(pos2)
    if np.abs(pos2_i - pos1_i) < 4:
        step = np.abs(pos2_i - pos1_i)
    else:
        step = 6 - np.abs(pos2_i - pos1_i)
    
    return step * ANGLE_BETWEEN_FILTERS




def calculate_wheel_rotation_time(angle):
    """calculate movement time for a single movement (seconds)"""
    movement_time = angle / WHEEL_ROTATION_SPEED
    return movement_time





def make_movement_dict(filter_positions):
    """make a dictionary of angles and times between all filters for a given filter wheel"""
    movements = {}
    
    for pos1 in filter_positions:
        for pos2 in filter_positions:
            if pos1 != pos2:
                angle = calculate_angle(filter_positions, pos1, pos2)
                movements[(pos1, pos2)] = {
                    "angle":angle,
                    "time":calculate_wheel_rotation_time(angle)
                }

    return movements







def mission_n_movements(n_movements, orbits_per_day):
    """calculate number of movements over the full 4 years"""

    n_orbits_mission = MISSION_DURATION * orbits_per_day
    n_movements_mission = n_orbits_mission * n_movements
    
    return n_movements_mission











#save headers to text file
filter_cycle_names = ["%s (%s)" %(k, ",".join(v["cycle"])) for k, v in filter_cycles.items()]
line = "\t".join(["Filter wheel positions"] + filter_cycle_names)

if SAVE_OUTPUT:
    with open("wheel_cycles.tsv", "w") as f:
        f.write("%s\n" %line)


#get all possible filter wheel position permutations, where dark is always filter 1
filter_position_permutations = wheel_filter_permutations(list(FILTER_NAMES))


#loop through permutations
for filter_positions in filter_position_permutations:
    
    line = ",".join(filter_positions)
    
    #loop through desired observation filter cycles
    for filter_cycle_name, filter_cycle_info in filter_cycles.items():
        
        
        filter_cycle = filter_cycle_info["cycle"]
        readout_time = filter_cycle_info["readout_time"]
        orbits_per_day = filter_cycle_info["orbits_per_day"]


        movement_dict = make_movement_dict(filter_positions)
        
        #start with a dark
        #iterate through observation filter cycle until observation ends
        #end with a dark
        
        loop = 0
        obs_time = 0
        obs_angle = 0
        pos1 = ""
        pos2 = ""
        
        #loop until observation time is filled
        while obs_time < OBSERVATION_TIME:
            
            if obs_time == 0:
                pos1 = "D" #dark first
                
                #if first observation is also a dark, skip it
                if filter_cycle[0] == "D":
                    loop = 1
                
            else:
                pos1 = pos2[:] #get previous filter position, set as starting position

            #reset filter cycle loop to zero if end of cycle
            if loop == len(filter_cycle):
                loop = 0
                
            pos2 = filter_cycle[loop] #set end position from next position in filter cycle
            loop += 1
            
            obs_time += readout_time #add integration + readout time
            
            if pos1 != pos2:
                #if starting and ending positions are not the same, add angle and rotation time
                obs_angle += movement_dict[(pos1, pos2)]["angle"]
                obs_time += movement_dict[(pos1, pos2)]["time"]
            
            if VERBOSE:
                print(pos1, pos2, obs_time, obs_angle)
            
        #add final dark
        pos1 = pos2[:]
        pos2 = "D"

        obs_time += readout_time
        
        if pos1 != pos2:
            obs_angle += movement_dict[(pos1, pos2)]["angle"]
            obs_time += movement_dict[(pos1, pos2)]["time"]

        
        obs_wheel_cycles = obs_angle / 360.0
        
        if VERBOSE:
            print(pos1, pos2, obs_time, obs_angle)
        
        if VERBOSE:
            print("Total 360 wheel cycles in observation", obs_wheel_cycles)
        
        mission_wheel_cycles = mission_n_movements(obs_wheel_cycles, orbits_per_day)
        
        if VERBOSE:
            print(line, "Total 360 wheel cycles in mission", mission_wheel_cycles)
        
        line += "\t%.0f" %mission_wheel_cycles
    
    if SAVE_OUTPUT:
        with open("wheel_cycles.tsv", "a") as f:
            f.write("%s\n" %line)