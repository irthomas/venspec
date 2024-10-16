# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:33:05 2022

@author: iant
"""


filter_parameters = {
    ("1", "n"): {"integration_time": 14.29, "colour": "red"},
    # ("1", "n"):{"integration_time":27., "colour":"red"},
    ("2", "d"): {"integration_time": 1.0, "colour": "blue"},
    ("2h", "d"): {"integration_time": 1.0, "colour": "deepskyblue"},
    ("2v", "d"): {"integration_time": 1.0, "colour": "lightskyblue"},
    ("2", "n"): {"integration_time": 14.29, "colour": "blue"},
    ("3", "n"): {"integration_time": 14.29, "colour": "yellow"},
    ("4", "d"): {"integration_time": 1.0, "colour": "green"},
    ("4h", "d"): {"integration_time": 1.0, "colour": "lime"},

    ("dk", "n"): {"integration_time": 14.29, "colour": "black"},

    # ("dk2", "d"): {"integration_time": 1.0, "colour": "black"},
    # ("dk4", "d"): {"integration_time": 1.0, "colour": "black"},
    ("dk", "d"): {"integration_time": 1.0, "colour": "black"},  # generic dark, use the shortest int time
}


filter_sequences = {
    "night_123d": [["1", "2", "3", "dk"], ["2", "3", "dk", "1"], ["3", "dk", "1", "2"], ["dk", "1", "2", "3"]],
    "day_22h44hd": [["2", "2h", "4", "4h", "dk"], ["4", "4h", "dk", "2", "2h"], ["dk", "2", "2h", "4", "4h"], ["2", "2h", "4", "4h", "dk"]],


    # "day_22h44h": ["2", "2h", "4", "4h"],
    # "day_24": ["2", "4"],
    # "day_22h": ["2", "2h"],
    # "day_44h": ["4", "4h"],

    # "day_22hd44hd": ["2", "2h", "dk2", "4", "4h", "dk4"],
    # "day_2d4d": ["2", "dk2", "4", "dk4"],
    # "day_24d": ["2", "4", "dk4"],
    # "day_22hd": ["2", "2h", "dk2"],
    # "day_44hd": ["4", "4h", "dk4"],
    # "day_22h4d": ["2", "2h", "4", "dk4"],

    # "night_123": ["1", "2", "3"],
    # "night_111d": ["1", "1", "1", "dk"],
    # "night_1d2d3d": ["1", "dk", "2", "dk", "3", "dk"],
}
