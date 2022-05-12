# -*- coding: utf-8 -*-
"""
Created on Thu May 12 09:45:43 2022

@author: iant

LOAD SPICE KERNELS
"""


import os
import spiceypy as sp
from datetime import datetime

from tools.file.paths import paths
# from tools.plotting.get_colours import get_colours

SPICE_DATETIME_FMT = "%Y %b %d %H:%M:%S.%f"
SHORT_DATETIME_FMT = "%d %b %Y"

SPICE_SHAPE_MODEL_METHOD = "Ellipsoid"
SPICE_METHOD = "INTERCEPT/ELLIPSOID"
SPICE_ABCORR = "NONE"
# SPICE_OBSRVR = "-668"

SPICE_TARGET = "VENUS"
SPICE_LONGITUDE_FORM = "PLANETOCENTRIC"
SPICE_PLANET_ID = 299 # venus



def load_spice_kernels(scenario):
    # SPICE_FORMATSTR = "C"
    # SPICE_PRECISION = 0
    
    kernel_root_dir = paths["KERNEL_ROOT_DIRECTORY"]
    
    #load each kernel manually
    
    kernels = [
        "fk/envision_v00.tf",
        "fk/envision_sc_venus_npo_v00.tf",
        "lsk/naif0012.tls",
        "pck/pck00010.tpc",
        "pck/de-403-masses.tpc",
        "spk/de432s.bsp",
        ]
    
    orbit_d ={
        "EnVision_ALT_T1_2032_NorthVOI.bsp":{
            "title":"ALT_T1_2032_NorthVOI (2022/01/18)\nAlternative backup scenario with Venus Orbit Insertion over the Northern hemisphere",
            "science_start":datetime(2034, 12, 20),
            "science_end":datetime(2038, 12, 18),
            "spice_observer":"-999", #error in kernel production => should be -668
            },
        "EnVision_ALT_T4_2032_SouthVOI.bsp":{
            "title":"ALT_T4_2032_SouthVOI (2022/01/18)\nAlternative baseline scenario with Venus Orbit Insertion over the Southern hemisphere",
            "science_start":datetime(2035, 3, 27),
            "science_end":datetime(2039, 3, 25),
            "spice_observer":"-999", #error in kernel production => should be -668
            },
        "EnVision_ESC_T4_2032_NorthVOI.bsp":{
            "title":"ESC_T4_2032_NorthVOI (2020/06/10)\nLaunch into direct escape at -5.0 deg declination and Venus Orbit Insertion over the Northern hemisphere.",
            "science_start":datetime(2035, 6, 15),
            # "science_end":datetime(2035, 9, 14),
            "science_end":datetime(2039, 6, 14),
            "spice_observer":"-668",
            },
        "EnVision_ESC_T2_2032_SouthVOI.bsp":{
            "title":"ESC_T2_2032_SouthVOI (2019/07/26)\nLaunch into direct escape at -3.0 deg declination and Venus Orbit Insertion over the Southern hemisphere",
            "science_start":datetime(2034, 11, 26),
            "science_end":datetime(2040, 3, 24),
            "spice_observer":"-668",
            },
        "EnVision_HEO_T2_2032_NorthVOI.bsp":{
            "title":"HEO_T2_2032_NorthVOI (2019/07/26)\nLaunch into HEO (310x100,000 km) and Venus Orbit Insertion over the Northern hemisphere",
            "science_start":datetime(2034, 11, 26),
            "science_end":datetime(2040, 3, 24),
            "spice_observer":"-668",
            },
        }
    

    if scenario in orbit_d.keys():
        for kernel in kernels:
            sp.furnsh(os.path.join(kernel_root_dir, os.path.normcase(kernel)))
            
        sp.furnsh(os.path.join(kernel_root_dir, "spk", scenario))
        
        return orbit_d[scenario]
    
    else:
        print("Error: scenario must be one of")
        for key in orbit_d.keys():
            print(key)

