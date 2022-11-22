# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 21:09:28 2022

@author: iant
"""

import numpy as np
from datetime import timedelta


from vs_obs.config.constants import PRECOOLING_DURATION, STANDBY_DURATION, DETECTOR_READOUT_TIME, SAVE_HTML
from vs_obs.planning.spice_functions import next_d2n_terminator, next_n2d_terminator, dt2et, et2dt, fov_lonlat, sc_alt
from vs_obs.filter_selection import obs_types, filter_parameters
from vs_obs.planning.filter_wheel_functions import calculate_angle, calculate_wheel_rotation_time
from vs_obs.html.plan_table_row import plan_table_row
from vs_obs.observation_cycle import observation_cycle
from vs_obs.html.save_plan_table import save_plan_table


def make_orbit_plan(utc_start_time, fov_vectors):

    


    orbit_plan = []
    h = ""
    
    #start at night to day terminator
    dt_start = next_n2d_terminator(utc_start_time)
    
    if SAVE_HTML:
        h += plan_table_row(dt_start, "Day to night terminator")
    
    orbit_plan.append([dt_start, {"daynight":"d", "status":"off"}])
    
    #find next sunset
    dt_next_sunset = next_d2n_terminator(dt_start)
    
    #find precooling start
    dt_standby_start = dt_next_sunset - timedelta(seconds=(PRECOOLING_DURATION + STANDBY_DURATION))
    orbit_plan.append([dt_standby_start, {"daynight":"d", "status":"standby"}])
    
    if SAVE_HTML:
        h += plan_table_row(dt_standby_start, "Standby start")
    
    
    dt_precooling_start = dt_next_sunset - timedelta(seconds=PRECOOLING_DURATION)
    orbit_plan.append([dt_precooling_start, {"daynight":"d", "status":"precooling"}])
    if SAVE_HTML:
        h += plan_table_row(dt_precooling_start, "Precooling start")
    
    
    
    #start time
    et = dt2et(dt_next_sunset)

    day_ix = 0
    night_ix = 0
    for obs_ix, observation in enumerate(observation_cycle):
        
        daynight = observation["daynight"]

        if observation["status"] == "off":
            #no observations - increment et
            if daynight == "n":
                #find next sunrise
                dt_next_sunrise = next_n2d_terminator(et2dt(et + 60.0))
                et = dt2et(dt_next_sunrise)
                night_ix += 1

            if daynight == "d":
                #find next sunrise
                dt_next_sunset = next_d2n_terminator(et2dt(et + 60.0))
                et = dt2et(dt_next_sunset)
                day_ix += 1
            continue
        

        filter_cycle_name = observation["filter_cycle"]
        filter_cycle = obs_types[filter_cycle_name] * 1000


        if daynight == "n":
            #find next sunrise, add 60 seconds to ensure the et is on the nightside
            dt_next_sunrise = next_n2d_terminator(et2dt(et + 60.0))


            
            #loop through nightside
            i = 0
            alts = []
            
            while et < dt2et(dt_next_sunrise):
                
                filter_name = filter_cycle[i + night_ix] #start with a different filter each time
                next_filter_name = filter_cycle[i + 1 + night_ix] #start with a different filter each time
                filter_dict = filter_parameters[(filter_name, daynight)]
                integration_time = filter_dict["integration_time"]
                
                #start of integration time
                lonlats_s = fov_lonlat(et, fov_vectors)
                if SAVE_HTML:
                    h += plan_table_row(et2dt(et), "Filter %s start (%0.2f, %0.2f)" %(filter_name, *lonlats_s[0, :]))

                et += integration_time

                #end of integration time
                lonlats_e = fov_lonlat(et, fov_vectors)
                if SAVE_HTML:
                    h += plan_table_row(et2dt(et), "Filter %s end (%0.2f, %0.2f)" %(filter_name, *lonlats_e[0, :]))
    
                
                #get altitude
                alt = sc_alt(et)
                alts.append(alt)
                
        
                #make footprint polygon including integration time
                footprint = np.zeros((5, 2))
                footprint[0, :] = lonlats_s[1, :]
                footprint[1, :] = lonlats_e[2, :]
                footprint[2, :] = lonlats_e[3, :]
                footprint[3, :] = lonlats_s[0, :]
                footprint[4, :] = lonlats_s[1, :]
        
                
                if "dk" not in filter_name:
                    orbit_plan.append([et2dt(et), {"daynight":daynight, "status":"science", \
                                                   "semiorbit":obs_ix, "alt":alt, "filter":filter_name, \
                                                   "footprint_s":lonlats_s, "footprint_e":lonlats_e, \
                                                   "footprint":footprint}])
                else:
                    orbit_plan.append([et2dt(et), {"daynight":daynight, "status":"dark"}])
            
            
            
                #filter wheel rotation
                #first change dk2 and dk4 to dk
                if "dk" in filter_name:
                    filter_name = "dk"
                if "dk" in next_filter_name:
                    next_filter_name = "dk"
                
                
                if filter_name != next_filter_name:
                    
                    filter_wheel_angle = calculate_angle(filter_name, next_filter_name)
                    if SAVE_HTML:
                        h += plan_table_row(et2dt(et), "Filter movement start %s" %filter_name)
                    et += calculate_wheel_rotation_time(filter_wheel_angle)
                    if SAVE_HTML:
                        h += plan_table_row(et2dt(et), "Filter movement end %s" %next_filter_name)
                else:
                    if SAVE_HTML:
                        h += plan_table_row(et2dt(et), "Detector readout start %s" %filter_name)
                    et += DETECTOR_READOUT_TIME
                    if SAVE_HTML:
                        h += plan_table_row(et2dt(et), "Detector readout end %s" %filter_name)
                
                i += 1
            # print(day_ix, night_ix, daynight, np.min(alts), np.mean(alts), np.max(alts))
                
                
            #pop last entry if over the sunrise
            if et > dt2et(dt_next_sunrise):
                orbit_plan.pop(-1)
                if SAVE_HTML:
                    h += plan_table_row(et2dt(et), "Last observation removed")
                
                
            night_ix += 1
        
        
        elif daynight == "d":
        
            """repeat for dayside"""
            #find next sunset, add 60 seconds to ensure the et is on the dayside
            dt_next_sunset = next_d2n_terminator(et2dt(et + 60.))
            
    
            #loop through dayside
            i = 0
            while et < dt2et(dt_next_sunset):
                
                filter_name = filter_cycle[i + day_ix] #start with a different filter each time
                next_filter_name = filter_cycle[i + 1 + day_ix] #start with a different filter each time
                filter_dict = filter_parameters[(filter_name, daynight)]
                integration_time = filter_dict["integration_time"]
                
                #start of integration time
                lonlats_s = fov_lonlat(et, fov_vectors)
    
    
                if SAVE_HTML:
                    h += plan_table_row(et2dt(et), "Filter %s start" %filter_name)
                et += integration_time
                if SAVE_HTML:
                    h += plan_table_row(et2dt(et), "Filter %s end" %filter_name)
    
                
                #end of integration time
                lonlats_e = fov_lonlat(et, fov_vectors)


                #get altitude
                alt = sc_alt(et)
                alts.append(alt)
                
                #make footprint polygon including integration time
                footprint = np.zeros((5, 2))
                footprint[0, :] = lonlats_s[1, :]
                footprint[1, :] = lonlats_e[2, :]
                footprint[2, :] = lonlats_e[3, :]
                footprint[3, :] = lonlats_s[0, :]
                footprint[4, :] = lonlats_s[1, :]


                if "dk" not in filter_name:
                    orbit_plan.append([et2dt(et), {"daynight":daynight, "status":"science", \
                                                   "semiorbit":obs_ix, "alt":alt, "filter":filter_name, \
                                                   "footprint_s":lonlats_s, "footprint_e":lonlats_e, \
                                                   "footprint":footprint}])
                else:
                    orbit_plan.append([et2dt(et), {"daynight":daynight, "status":"dark"}])
            
        
        
                #filter wheel rotation
                #first change dk2 and dk4 to dk
                if "dk" in filter_name:
                    filter_name = "dk"
                if "dk" in next_filter_name:
                    next_filter_name = "dk"
                
                
                if filter_name != next_filter_name:
                    
                    filter_wheel_angle = calculate_angle(filter_name, next_filter_name)
                    if SAVE_HTML:
                        h += plan_table_row(et2dt(et), "Filter movement start %s" %filter_name)
                    et += calculate_wheel_rotation_time(filter_wheel_angle)
                    if SAVE_HTML:
                        h += plan_table_row(et2dt(et), "Filter movement end %s" %next_filter_name)
                else:
                    if SAVE_HTML:
                        h += plan_table_row(et2dt(et), "Detector readout start %s" %filter_name)
                    et += DETECTOR_READOUT_TIME
                    if SAVE_HTML:
                        h += plan_table_row(et2dt(et), "Detector readout end %s" %filter_name)
                
                i += 1
            # print(day_ix, night_ix, daynight, np.min(alts), np.mean(alts), np.max(alts))
                
                
            #pop last entry if over the sunrise
            if et > dt2et(dt_next_sunset):
                orbit_plan.pop(-1)
                if SAVE_HTML:
                    h += plan_table_row(et2dt(et), "Last observation removed")

            day_ix += 1


        #add 20 seconds of standby
        if SAVE_HTML:
            h += plan_table_row(et2dt(et), "Standby start")
        et += 20.
        if SAVE_HTML:
           h += plan_table_row(et2dt(et), "Standby end")

        orbit_plan.append([et2dt(et), {"daynight":"", "status":"standby"}])
        
        if SAVE_HTML:
            save_plan_table(h, "obs_%s_%02i" %(daynight, obs_ix))
            h = "" #clear output
        
    
    return orbit_plan, h