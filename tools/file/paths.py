# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 09:14:10 2020

@author: iant

PATHS

"""
import os
import platform

if platform.system() == "Windows":
    SYSTEM = "Windows"
else:
    SYSTEM = "Linux"

paths = {}

if SYSTEM == "Linux":  # linux system
    paths["BASE_DIRECTORY"] = "/bira-iasb/projects/EnVision/OPERATIONS/PLANNING/Python"
    paths["KERNEL_ROOT_DIRECTORY"] = "/bira-iasb/projects/EnVision/OPERATIONS/PLANNING/envision_kernels/envision/kernels"
    paths["KERNEL_DIRECTORY"] = "/bira-iasb/projects/EnVision/OPERATIONS/PLANNING/envision_kernels/envision/kernels/mk"
    paths["REFERENCE_DIRECTORY"] = os.path.join(paths["BASE_DIRECTORY"], "reference_files")


elif os.path.exists(os.path.normcase(r"C:\Users\iant\Dropbox\VenSpec\Python")):  # outside BIRA
    paths["BASE_DIRECTORY"] = os.path.normcase(r"C:\Users\iant\Dropbox\VenSpec\Python")
    paths["KERNEL_ROOT_DIRECTORY"] = os.path.normcase(r"C:\Users\iant\Documents\DATA\envision_kernels\envision\kernels")
    paths["KERNEL_DIRECTORY"] = os.path.join(paths["KERNEL_ROOT_DIRECTORY"], "mk")
    paths["REFERENCE_DIRECTORY"] = os.path.join(paths["BASE_DIRECTORY"], "reference_files")
