# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 12:26:06 2022

@author: iant

VENSPEC INSTRUMENT LINE SHAPE
"""

import numpy as np
import matplotlib.pyplot as plt




def gaussian(x, a, b, c, d):
    #amplitude is a variable (a)
    return a * np.exp(-(x - b)**2.0/(2.0 * c**2.0)) + d


def convolve_ils(lambda_um, delta_lambda_um, ils_gaussian_width, venus_um_grid, venus_radiance, ax=None):

    #convert spectral resolution to gauss FWHM
    gauss_width = delta_lambda_um / 2.355
    
    #calculate wavelength range of gaussian and make grid
    um_grid_start = -1.0 * gauss_width * ils_gaussian_width
    um_grid_end = gauss_width * ils_gaussian_width
    
    #find indices in hr grid covering +- wavelength range of gaussian
    nu_hr_start_ix = np.searchsorted(venus_um_grid, lambda_um + um_grid_start) #start index
    nu_hr_end_ix = np.searchsorted(venus_um_grid, lambda_um + um_grid_end) #end index
    
    if nu_hr_start_ix == nu_hr_end_ix:
        #pixel ILS is beyond the range of the input spectrum -> using nearest value instead
        px_venus_convolved = venus_radiance[-1]

        warning = True
        
    else:
        
        px_grid = venus_um_grid[nu_hr_start_ix:nu_hr_end_ix] - lambda_um
        px_venus_signal = venus_radiance[nu_hr_start_ix:nu_hr_end_ix]
        
        ils_gauss = gaussian(px_grid, 1.0, 0.0, gauss_width, 0.0)
        ils_venus_signal = px_venus_signal * ils_gauss
        
        px_venus_convolved = np.sum(ils_venus_signal) / np.sum(ils_gauss)
        
        warning = False

        # test ILS
        if ax:
            ax[0].plot(px_grid + lambda_um, ils_gauss, "k")
            ax[1].plot(px_grid + lambda_um, px_venus_signal, "orange")
            # ax[1].plot(px_grid + lambda_um, ils_venus_signal)
            # ax[1].plot([px_grid[0] + lambda_um, px_grid[-1] + lambda_um], [px_venus_convolved, px_venus_convolved])  
        # stop()
    
    return px_venus_convolved, warning


def detector_ils(lambda_ums, delta_lambda_um, ils_gaussian_width, venus_um_grid, venus_radiance, plot_ils=False):
    
    venus_signal_convolved = np.zeros_like(lambda_ums)
    
    warning_shown = False
    if plot_ils:
        fig_ils, ax_ils = plt.subplots(figsize=(17, 5), constrained_layout=True)
        ax_ils2 = ax_ils.twinx()
        ax_ils.set_title("ILS vs radiance")
        ax_ils.set_xlabel("Wavelength (um)")
        ax_ils.set_ylabel("Normalised ILS")
        ax_ils2.set_yscale("log")
        ax_ils2.set_ylabel("Venus radiance (log scale)")
        ax = [ax_ils, ax_ils2]
    for px_ix, lambda_um in enumerate(lambda_ums):
        if plot_ils:
            venus_signal_convolved[px_ix], warning = convolve_ils(lambda_um, delta_lambda_um, ils_gaussian_width, venus_um_grid, venus_radiance, ax=ax)
        else:
            venus_signal_convolved[px_ix], warning = convolve_ils(lambda_um, delta_lambda_um, ils_gaussian_width, venus_um_grid, venus_radiance)

        #display error just once
        if warning and not warning_shown:
            print("Warning: pixel ILS is beyond the range of the input spectrum -> using nearest value instead")
            warning_shown = True

    
    return venus_signal_convolved



# test ILS plotting - whole detector
# first run snr_model.py to get the snr object, then run this script to plot the ILS
# if __name__ == "__main__":
#     lambda_ums = snr.px_um
#     delta_lambda_um = snr.px_delta_lambda_um
#     # delta_lambda_um = snr.real_px_delta_lambda_um
#     ils_gaussian_width = snr.ils_gaussian_width

#     venus_um_grid = snr.night_um
#     venus_radiance = snr.night_Wm2um
    
    
    
#     fig, ax = plt.subplots()
#     ax.set_ylabel("Gaussian ILS")
#     ax.set_xlabel("Wavelength (um)")
#     ax.set_title("Detector pixel ILS")
    
#     for px_ix, lambda_um in enumerate(lambda_ums):
    
#         #convert spectral resolution to gauss FWHM
#         gauss_width = delta_lambda_um / 2.355
        
#         #calculate wavelength range of gaussian and make grid
#         um_grid_start = -1.0 * gauss_width * ils_gaussian_width
#         um_grid_end = gauss_width * ils_gaussian_width
        
#         #find indices in hr grid covering +- wavelength range of gaussian
#         nu_hr_start_ix = np.searchsorted(venus_um_grid, lambda_um + um_grid_start) #start index
#         nu_hr_end_ix = np.searchsorted(venus_um_grid, lambda_um + um_grid_end) #end index
    
#         px_grid = venus_um_grid[nu_hr_start_ix:nu_hr_end_ix] - lambda_um
#         px_venus_signal = venus_radiance[nu_hr_start_ix:nu_hr_end_ix]
        
#         ils_gauss = gaussian(px_grid, 1.0, 0.0, gauss_width, 0.0)
#         ils_venus_signal = px_venus_signal * ils_gauss
        
#         px_venus_convolved = np.sum(ils_venus_signal) / np.sum(ils_gauss)
    
#         ax.plot(px_grid + lambda_um, ils_gauss)
    
