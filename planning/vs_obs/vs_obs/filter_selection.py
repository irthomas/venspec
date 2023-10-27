# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:33:05 2022

@author: iant
"""



filter_parameters = {
    ("1", "n"):{"integration_time":14.29, "colour":"red"},
    # ("1", "n"):{"integration_time":27., "colour":"red"},
    ("2", "d"):{"integration_time":1.2, "colour":"blue"},
    ("2h", "d"):{"integration_time":1.2, "colour":"deepskyblue"},
    ("2v", "d"):{"integration_time":1.2, "colour":"lightskyblue"},
    ("2", "n"):{"integration_time":14.29, "colour":"blue"},
    ("3", "n"):{"integration_time":14.29, "colour":"yellow"},
    ("4", "d"):{"integration_time":0.44, "colour":"green"},
    ("4h", "d"):{"integration_time":0.44, "colour":"lime"},

    ("dk", "n"):{"integration_time":14.29, "colour":"black"},
    ("dk2", "d"):{"integration_time":1.2, "colour":"black"},
    ("dk4", "d"):{"integration_time":0.44, "colour":"black"},
}



obs_types = {
    "day_24_polarisation_h":["2", "2h", "4", "4h"],
    "day_24":["2", "4"],
    "day_2_polarisation":["2", "2h"],
    "day_4_polarisation":["4", "4h"],

    "day_2d4d_polarisation_h":["2", "2h", "dk2", "4", "4h", "dk4"],
    "day_2d4d":["2", "dk2", "4", "dk4"],
    "day_2d_polarisation":["2", "2h", "dk2"],
    "day_4d_polarisation":["4", "4h", "dk4"],

    "night_123":["1", "2", "3"],
    "night_123d":["1", "2", "3", "dk"],
    "night_111d":["1", "1", "1", "dk"],
    "night_1d2d3d":["1", "dk", "2", "dk", "3", "dk"],
}

