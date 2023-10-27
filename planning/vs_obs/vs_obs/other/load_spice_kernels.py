# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:30:37 2022

@author: iant
"""
import os
# from datetime import datetime
import spiceypy as sp


from vs_obs.config.paths import paths



def load_spice_kernels():
    kernel_root_dir = paths["KERNEL_ROOT_DIRECTORY"]
    
    
    
    
    #load each kernel manually
    
    kernels = [
        "fk/envision_v00.tf",
        "fk/envision_sc_venus_npo_v00.tf",
        "envision_venspec_v01.ti", #nadir along track orientation
        "lsk/naif0012.tls",
        "pck/pck00010.tpc",
        "pck/de-403-masses.tpc",
        "sclk/envision_200715_fict.tsc",
        "spk/de432s.bsp",
        ]
    
    # orbit_d ={
    #     "EnVision_ESC_T2_2032_SouthVOI.bsp":{
    #         "title":"ESC_T2_2032_SouthVOI (2019/07/26)\nLaunch into direct escape at -3.0 deg declination and Venus Orbit Insertion over the Southern hemisphere",
    #         "science_start":datetime(2034, 11, 26),
    #         "science_end":datetime(2040, 3, 24),
    #         "spice_observer":"-668",
    #         },
    #     }
    
    
    kernel_name = "EnVision_ET1_2031_NorthVOI.bsp" #baseline launch
    # spice_observer = orbit_d[kernel_name]["spice_observer"]
    
    print("KERNEL_DIRECTORY=%s, METAKERNEL_NAME=%s" %(kernel_root_dir, kernel_name))
    
    for kernel in kernels:
        sp.furnsh(os.path.join(kernel_root_dir, os.path.normcase(kernel)))
    sp.furnsh(os.path.join(kernel_root_dir, "spk", kernel_name))
            
