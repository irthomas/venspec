# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 14:46:30 2022

@author: iant

MAKE VERY HIGH SPECTRAL RESOLUTION VENUS NADIR SPECTRA

"""
import os
import numpy as np
import matplotlib.pyplot as plt

from tools.file.paths import paths
from tools.spectra.fft_zerofilling import fft_hr_nu_spectrum


def hr_nadir_spectra(zerofilling):
    
    # day_filename = "Simulations_Haus_LIDORT_dayside_scattering_molecule_aer_Rayleigh_sp_nonoise.dat"
    # night_filename = "Simulations_Haus_LIDORT_nightside_scattering_molecule_aer_Rayleigh_sp_nonoise.dat"

    day_filename = "Simulations_LIDORT_dayside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"
    night_filename = "Combined_spectrum_modified.dat"
    # night_filename = "Simulations_LIDORT_nightside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"
    
    day_nm, day_Wm2um = np.loadtxt(
        os.path.join(paths["BASE_DIRECTORY"], "reference_files", day_filename), skiprows=4, unpack=True)
    night_nm, night_Wm2um = np.loadtxt(
        os.path.join(paths["BASE_DIRECTORY"], "reference_files", night_filename), skiprows=4, unpack=True)
    
    day_um = np.flip(day_nm) / 1000.
    day_Wm2um = np.flip(day_Wm2um)
    
    night_um = np.flip(night_nm) / 1000.
    night_Wm2um = np.flip(night_Wm2um)
    
    if zerofilling > 1:
        day_um_hr, day_wm2um_hr = fft_hr_nu_spectrum(day_um, day_Wm2um, zerofilling=zerofilling)
        night_um_hr, night_wm2um_hr = fft_hr_nu_spectrum(night_um, night_Wm2um, zerofilling=zerofilling)
        
        return {"day_um":day_um_hr, "day_Wm2um":day_wm2um_hr, "night_um":night_um_hr, "night_Wm2um":night_wm2um_hr}

    else:
        return {"day_um":day_um, "day_Wm2um":day_Wm2um, "night_um":night_um, "night_Wm2um":night_Wm2um}


def extend_to_2p6um(nadir):
    """extend to 2.6um"""
    um = nadir["day_um"]
    um_delta = um[1] - um[0]
    um_extra = np.arange(um[-1]+um_delta, 2.6000000001, um_delta)
    new_um = np.concatenate((um, um_extra))
    nadir["day_um"] = new_um

    rad = nadir["day_Wm2um"]
    rad_extra = np.zeros_like(um_extra) + rad[-1]
    new_rad = np.concatenate((rad, rad_extra))
    nadir["day_Wm2um"] = new_rad

    plt.plot(new_um, new_rad)
    plt.plot(um, rad, alpha=0.5)

    um = nadir["night_um"]
    um_delta = um[1] - um[0]
    um_extra = np.arange(um[-1]+um_delta, 2.6000000001, um_delta)
    new_um = np.concatenate((um, um_extra))
    nadir["night_um"] = new_um

    rad = nadir["night_Wm2um"]
    rad_extra = np.zeros_like(um_extra) + rad[-1]
    new_rad = np.concatenate((rad, rad_extra))
    nadir["night_Wm2um"] = new_rad
    
    plt.plot(new_um, new_rad)
    plt.plot(um, rad, alpha=0.5)
    plt.yscale("log")
    
    return nadir




# nadir = hr_nadir_spectra(1)
# nadir = extend_to_2p6um(nadir)