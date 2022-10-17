# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 11:14:55 2022

@author: iant
"""

import numpy as np


def spectral_calibration(band):

    n_columns = 384 #number of detector spectral pixels
    # fno = 3.0 #f/# number
    fsr = 144.88 #cm-1
    theoretical_resolving_power = 11000.
    real_resolving_power = 7000.
    slit_width = 1. #pixels

        
    #desired spectral range for each band
    um_range_aim = {
        "1":[1.16, 1.18],
        "2a":[2.34, 2.422],
        "2b":[2.422, 2.51],
        "3":[1.704, 1.747],
        "4":[1.367, 1.394],
        }[band]
    
    #find corresponding diffraction order
    order_calc = (10000. / np.mean(um_range_aim)) / fsr
    order = int(np.round(order_calc))
    
    #wavelength of centre of diffraction order
    band_centre_um = (10000. / order) / fsr
    
    cm_centre = 10000. / band_centre_um
    
    #diffration order wavenumber range: centre +- half the FSR
    order_range_cm = np.array([cm_centre + fsr/2, cm_centre - fsr/2])
    order_range_um = 10000. / order_range_cm
    
    
    
    #spectral resolution: width of gaussian encompassing wavelengths received by a pixel * 2.355
    #resolving power = lambda/delta_lambda
    px_delta_lambda_um = band_centre_um / theoretical_resolving_power
    real_px_delta_lambda_um = band_centre_um / real_resolving_power #for ILS calculation, use reduced RP 
    
    #spectral sampling: wavelength between centre of adjacent pixels
    px_sampling_um = px_delta_lambda_um / slit_width
    
    #total wavelength range allowed on detector
    detector_um_width = px_sampling_um * n_columns
    detector_um_range = np.array([
        band_centre_um - detector_um_width/2., 
        band_centre_um + detector_um_width/2.
        ])
    
    #if FSR is bigger than detector, light cannot be detected => reduce band spectral range
    #if detector width is bigger than FSR, extra pixels will not be illuminated -> reduce band spectral range
    band_um_range = np.array([
        np.max([detector_um_range[0], order_range_um[0]]),
        np.min([detector_um_range[1], order_range_um[1]]),
        ])

    #centre wavelengths on each pixel even if outside FSR
    detector_px_centres = np.arange(detector_um_range[0], detector_um_range[-1]-0.000001, px_sampling_um)

    #centre wavelengths on each pixel (real detector)
    px_illuminated_ixs = np.where((detector_px_centres >= band_um_range[0]) & (detector_px_centres <= band_um_range[1])) #indices of illuminated pixels
    px_um = detector_px_centres[px_illuminated_ixs] #real detector wavelengths
    # print("Band", band, "pixels used:", len(bands[band]["detector_ums"]))
        
    return px_um



spectral_calibration("1")