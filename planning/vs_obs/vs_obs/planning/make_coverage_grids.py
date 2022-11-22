# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:36:42 2022

@author: iant
"""

import numpy as np
from itertools import product

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from vs_obs.config.constants import FILTER_WHEEL_LAYOUT


def make_coverage_grids(orbit_plan, daynight, semiorbits=None):
    #footprint statistics
    #make HR lon lat grid
    lat_minmax = [-60., 60.]
    lon_minmax = [-180., 180.]
    
    # delta_lat = 0.1
    # delta_lon = 0.1
    delta_lat = 0.25
    delta_lon = 0.25
    
    lats_hr = np.arange(lat_minmax[0], lat_minmax[1], delta_lat)
    lons_hr = np.arange(lon_minmax[0], lon_minmax[1], delta_lon)
    coverage_grids = {filter_:np.zeros((len(lats_hr), len(lons_hr))) for filter_ in FILTER_WHEEL_LAYOUT}
    
    
    
    
    coverage_coords = np.asfarray(list(product(lons_hr, lats_hr)))
    
    matching_indices = [i for i in range(len(orbit_plan)) if orbit_plan[i][1]["status"] == "science"]
    matching_indices = [i for i in matching_indices if orbit_plan[i][1]["daynight"] == daynight]
    
    if semiorbits:
        matching_indices = [i for i in matching_indices if orbit_plan[i][1]["semiorbit"] in semiorbits]
        
        
        
    
    polygons = [orbit_plan[i][1]["footprint"] for i in matching_indices]
    filters = [orbit_plan[i][1]["filter"] for i in matching_indices]
    
    
    #get all lats from all footprints
    all_polygon_lats = [p[:, 1] for p in polygons]
    
    #find indices of all polygons at mid latitudes
    indices_in_lat_zone = [i for i, lats in enumerate(all_polygon_lats) if (np.min(lats) < 60.0) & (np.max(lats) > -60.0)]
    
    #find all polygoins in the latitude zone
    polygons_in_lat_zone = [polygons[i] for i in indices_in_lat_zone]
    filters_in_lat_zone = [filters[i] for i in indices_in_lat_zone]
    
    
    
    for polygon, filter_ in zip(polygons_in_lat_zone, filters_in_lat_zone):
    
    
        #get min and max lons of polygon
        min_polygon_lon = np.min(polygon[:, 0])
        max_polygon_lon = np.max(polygon[:, 0])
        min_polygon_lat = np.min(polygon[:, 1])
        max_polygon_lat = np.max(polygon[:, 1])
    
    
        #get indices of coordinates in the required longitude zone
        indices_in_lonlat_zone = np.where((coverage_coords[:, 0] > (min_polygon_lon - delta_lon/2)) & (coverage_coords[:, 0] < (max_polygon_lon + delta_lon/2)) & \
                                          (coverage_coords[:, 1] > (min_polygon_lat - delta_lat/2)) & (coverage_coords[:, 1] < (max_polygon_lat + delta_lat/2)))[0]
    
        for i in indices_in_lonlat_zone:
        
            coverage_coord = coverage_coords[i, :] + [delta_lon/2, -delta_lat/2]
            #check if centre of each grid point is within each footprint
            if Polygon(polygon).contains(Point(coverage_coord)):
                #calculate indices in coverage grid
                lat_ix = len(lats_hr) - np.mod(i, len(lats_hr))
                lon_ix = i // len(lats_hr)
                # print(coverage_coord, lat_ix, lon_ix)
                #record each filter with a different number
                
                if lat_ix < len(lats_hr) and lon_ix < len(lons_hr):
                    
                    coverage_grids[filter_][lat_ix, lon_ix] = 1
    

    return coverage_grids, lats_hr, lons_hr