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
    paths["BASE_DIRECTORY"] = "/bira-iasb/projects/EnVision/OPERATIONS/PLANNING/Python/venspec/"
    paths["PLANNING_DIRECTORY"] = os.path.join(paths["BASE_DIRECTORY"], "Planning", "vs_obs")
    paths["REFERENCE_DIRECTORY"] = os.path.join(paths["BASE_DIRECTORY"], "reference_files")

    paths["KERNEL_ROOT_DIRECTORY"] = "/bira-iasb/projects/EnVision/OPERATIONS/PLANNING/envision_kernels/envision/kernels/"
    paths["KERNEL_DIRECTORY"] = os.path.join(paths["KERNEL_ROOT_DIRECTORY"], "mk")


elif os.path.exists(os.path.normcase(r"C:\Users\iant\Dropbox\VenSpec\Python")):  # outside BIRA

    paths["BASE_DIRECTORY"] = os.path.normcase(r"C:\Users\iant\Dropbox\VenSpec\Python")
    paths["PLANNING_DIRECTORY"] = os.path.normcase(r"C:\Users\iant\Dropbox\VenSpec\Python\Planning\vs_obs")
    paths["REFERENCE_DIRECTORY"] = os.path.normcase(r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files")

    paths["KERNEL_ROOT_DIRECTORY"] = os.path.normcase(r"C:\Users\iant\Documents\DATA\envision_kernels\envision")
    paths["KERNEL_DIRECTORY"] = os.path.normcase(r"C:\Users\iant\Documents\DATA\envision_kernels\envision\mk")
