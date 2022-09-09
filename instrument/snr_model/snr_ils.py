# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 12:26:06 2022

@author: iant

VENSPEC INSTRUMENT LINE SHAPE
"""

import numpy as np
import matplotlib.pyplot as plt

from input_spectra import load_asimut_spectra


day_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\Simulations_LIDORT_dayside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"
night_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\Simulations_LIDORT_nightside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"


class test(object):
    
    def __init__(self):

        
        
        self.day_input_filepath = day_input_filepath
        self.night_input_filepath = night_input_filepath
        
        load_asimut_spectra(self)


snr = test()

um_grid = snr.night_um
venus_signal = snr.night_Wm2um

grid_width_scalar = 3.0 #calculate ILS for range (-delta lambda * scalar) to (+delta lambda * scalar)


def gaussian(x, a, b, c, d):
    return a * np.exp(-((x - b)/c)**2.0) + d


def convolve_ils(lambda_um, delta_lambda_um):

    #convert spectral resolution to gauss FWHM
    gauss_width = delta_lambda_um / 2.355
    
    #calculate wavelength range of gaussian and make grid
    um_grid_start = -1.0 * gauss_width * grid_width_scalar
    um_grid_end = gauss_width * grid_width_scalar
    
    #find indices in hr grid covering +- wavelength range of gaussian
    nu_hr_start_ix = np.searchsorted(um_grid, lambda_um + um_grid_start) #start index
    nu_hr_end_ix = np.searchsorted(um_grid, lambda_um + um_grid_end) #end index
    
    px_grid = um_grid[nu_hr_start_ix:nu_hr_end_ix] - lambda_um
    px_venus_signal = venus_signal[nu_hr_start_ix:nu_hr_end_ix]
    
    ils_gauss = gaussian(px_grid, 1.0, 0.0, gauss_width, 0.0)
    ils_venus_signal = px_venus_signal * ils_gauss
    
    px_venus_signal[255]
    px_venus_convolved = np.sum(ils_venus_signal) / np.sum(ils_gauss)

    fig, ax = plt.subplots()
    ax.plot(px_grid + lambda_um, ils_gauss)
    ax2 = ax.twinx()
    ax2.plot(px_grid + lambda_um, px_venus_signal)
    ax2.plot(px_grid + lambda_um, ils_venus_signal)
    ax2.plot([px_grid[0] + lambda_um, px_grid[-1] + lambda_um], [px_venus_convolved, px_venus_convolved])    
    
    return px_venus_convolved


delta_lambda_um = 0.002
lambda_um = 2.4


convolve_ils(lambda_um, delta_lambda_um)