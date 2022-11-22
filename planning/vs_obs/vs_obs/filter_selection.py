# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:33:05 2022

@author: iant
"""



filter_parameters = {
    # ("1", "n"):{"integration_time":13.5, "colour":"red"},
    ("1", "n"):{"integration_time":27., "colour":"red"},
    ("2", "d"):{"integration_time":0.982, "colour":"blue"},
    ("2h", "d"):{"integration_time":0.982, "colour":"deepskyblue"},
    ("2v", "d"):{"integration_time":0.982, "colour":"lightskyblue"},
    ("2", "n"):{"integration_time":13.5, "colour":"blue"},
    ("3", "n"):{"integration_time":13.5, "colour":"yellow"},
    ("4", "d"):{"integration_time":0.370, "colour":"green"},
    ("4h", "d"):{"integration_time":0.370, "colour":"lime"},

    ("dk", "n"):{"integration_time":13.5},
    ("dk2", "d"):{"integration_time":0.982},
    ("dk4", "d"):{"integration_time":0.370},
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

