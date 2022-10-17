# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:37:54 2022

@author: iant
"""




def make_orbit_plan(utc_start_time):
    
    
    
    orbit_plan = []
    
    #start at night to day terminator
    dt_start = next_n2d_terminator(utc_start_time)
    
    
    orbit_plan.append([dt_start, {"daynight":"d", "status":"off"}])
    
    #find next sunset
    dt_next_sunset = next_d2n_terminator(dt_start)
    
    #find precooling start
    dt_standby_start = dt_next_sunset - timedelta(seconds=(PRECOOLING_DURATION + STANDBY_DURATION))
    orbit_plan.append([dt_standby_start, {"daynight":"d", "status":"standby"}])
    
    
    dt_precooling_start = dt_next_sunset - timedelta(seconds=PRECOOLING_DURATION)
    orbit_plan.append([dt_precooling_start, {"daynight":"d", "status":"precooling"}])
    
    
    
    #start time
    et = dt2et(dt_next_sunset)
    
    
    for orbit_loop in range(3):
    
    
        #find next sunrise
        dt_next_sunrise = next_n2d_terminator(et2dt(et))
        
        
        #loop through nightside
        daynight = "n"
        filter_cycle = obs_type_dict["night_123d"] * 1000
        
        i = 0
        while et < dt2et(dt_next_sunrise):
            
            filter_name = filter_cycle[i]
            next_filter_name = filter_cycle[i+1]
            filter_dict = filters_dict[(filter_name, daynight)]
            integration_time = filter_dict["integration_time"]
            
            #start of integration time
            lonlats_s = fov_lonlat(et, fov_vectors)
            et += integration_time
            
            #get altitude
            # alt = sc_alt(et)
            # print(alt)
            
            #end of integration time
            lonlats_e = fov_lonlat(et, fov_vectors)
    
    
    
            
            if "dk" not in filter_name:
                orbit_plan.append([et2dt(et), {"daynight":daynight, "status":"science", "filter":filter_name, "footprint_s":lonlats_s, "footprint_e":lonlats_e}])
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
                et += calculate_wheel_rotation_time(filter_wheel_angle)
            else:
                et += DETECTOR_READOUT_TIME
            
            i += 1
            
            
        #pop last entry if over the sunrise
        if et > dt2et(dt_next_sunrise):
            orbit_plan.pop(-1)
        
        
        et += 20.
        #add 20 seconds of standby
        orbit_plan.append([et2dt(et), {"daynight":"d", "status":"standby"}])
        
        
        
        
        """repeat for dayside"""
        #find next sunset, add 100 seconds to ensure the et is on the dayside
        dt_next_sunset = next_d2n_terminator(et2dt(et + 100.))
        
        
        #loop through dayside
        daynight = "d"
        filter_cycle = obs_type_dict["day_2d4d"] * 1000
        
        #start time
        i = 0
        while et < dt2et(dt_next_sunset):
            
            filter_name = filter_cycle[i]
            next_filter_name = filter_cycle[i+1]
            filter_dict = filters_dict[(filter_name, daynight)]
            integration_time = filter_dict["integration_time"]
            
            #start of integration time
            lonlats_s = fov_lonlat(et, fov_vectors)
            et += integration_time
            
            #end of integration time
            lonlats_e = fov_lonlat(et, fov_vectors)
            
            if "dk" not in filter_name:
                orbit_plan.append([et2dt(et), {"daynight":daynight, "status":"science", "filter":filter_name, "footprint_s":lonlats_s, "footprint_e":lonlats_e}])
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
                et += calculate_wheel_rotation_time(filter_wheel_angle)
            else:
                et += DETECTOR_READOUT_TIME
            
            i += 1
            
            
        #pop last entry if over the sunrise
        if et > dt2et(dt_next_sunset):
            orbit_plan.pop(-1)
