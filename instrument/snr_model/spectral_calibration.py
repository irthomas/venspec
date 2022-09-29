# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 10:32:35 2022

@author: iant
"""

import numpy as np



def spectral_calibration(self, band):

        
    #desired spectral range for each band
    um_range_aim = {
        "1":[1.16, 1.18],
        "2a":[2.34, 2.422],
        "2b":[2.422, 2.51],
        "3":[1.704, 1.747],
        "4":[1.367, 1.394],
        }[band]
    
    #find corresponding diffraction order
    order_calc = (10000. / np.mean(um_range_aim)) / self.fsr
    self.order = int(np.round(order_calc))
    
    #wavelength of centre of diffraction order
    self.band_centre_um = (10000. / self.order) / self.fsr
    
    self.cm_centre = 10000. / self.band_centre_um
    
    #diffration order wavenumber range: centre +- half the FSR
    order_range_cm = np.array([self.cm_centre + self.fsr/2, self.cm_centre - self.fsr/2])
    order_range_um = 10000. / order_range_cm
    
    
    
    #spectral resolution: width of gaussian encompassing wavelengths received by a pixel * 2.355
    #resolving power = lambda/delta_lambda
    self.px_delta_lambda_um = self.band_centre_um / self.theoretical_resolving_power
    self.real_px_delta_lambda_um = self.band_centre_um / self.real_resolving_power #for ILS calculation, use reduced RP 
    
    #spectral sampling: wavelength between centre of adjacent pixels
    px_sampling_um = self.px_delta_lambda_um / self.slit_width
    
    #total wavelength range allowed on detector
    detector_um_width = px_sampling_um * self.n_columns
    detector_um_range = np.array([
        self.band_centre_um - detector_um_width/2., 
        self.band_centre_um + detector_um_width/2.
        ])
    
    #if FSR is bigger than detector, light cannot be detected => reduce band spectral range
    #if detector width is bigger than FSR, extra pixels will not be illuminated -> reduce band spectral range
    band_um_range = np.array([
        np.max([detector_um_range[0], order_range_um[0]]),
        np.min([detector_um_range[1], order_range_um[1]]),
        ])

    #centre wavelengths on each pixel even if outside FSR
    self.detector_px_centres = np.arange(detector_um_range[0], detector_um_range[-1]-0.000001, px_sampling_um)

    #centre wavelengths on each pixel (real detector)
    self.px_illuminated_ixs = np.where((self.detector_px_centres >= band_um_range[0]) & (self.detector_px_centres <= band_um_range[1])) #indices of illuminated pixels
    self.px_um = self.detector_px_centres[self.px_illuminated_ixs] #real detector wavelengths
    # print("Band", band, "pixels used:", len(bands[band]["detector_ums"]))
        
    return self