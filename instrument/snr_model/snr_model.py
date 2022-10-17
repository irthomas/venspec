# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 10:12:56 2022

@author: iant


BIRA SNR MODEL

AIMS:
1) APPROXIMATELY REPRODUCE THE PHASE A OIP SNR STUDY RESULTS. RESULTS WITHIN A FEW % ARE ACCEPTABLE FOR COMPARISON AT THIS STAGE
2) WHEN THE PHASE B SNR MODEL IS FINALISED, THE RESULTS MUST MATCH EXACTLY
2) EXPAND THE MODEL TO TAKE REAL HIGH RESOLUTION INPUT SPECTRA, CONVOLVE TO DETECTOR PIXELS, AND CALCULATE ERROR ON EACH PIXEL RATHER THAN MEAN SNR
3) USE THE MODEL TO INVESTIGATE THE BEST OBSERVATION PARAMETERS FOR VENSPEC-H OPERATIONS


NOT YET IMPLEMENTED

CONTRIBUTION OF WARM SECTION (I.E. RADIATION PASSING THROUGH THE SLIT) TO THE NOISE CALCULATION (LOW IF WARM SECTION NOT TOO HOT)
RUNNING SIMULATION FOR BOTH DETECTOR GAIN SETTINGS AND CHOOSING THE BEST ONE
SIMULATE REAL SPECTRAL ACQUISITIONS: CALCULATE VENUS FRAME AND DARK FRAME SEPARATELY, DIGITISE OUTPUT, AND THEN SUBTRACT


    
"""

import numpy as np


from spectral_calibration import spectral_calibration
from snr_detector import detector_settings, detector, detector_qe_interp
from snr_signal import venus, signal
from snr_thermal_background import cold_shield, detector_window, cold_section
from save_output import prepare_plot, save_plot#, save_generic_output
from save_output_ASIMUT import save_cm1_clean_radiance, save_cm1_radiance_error#, save_cm1_radiance_with_noise
from input_spectra import load_asimut_spectra



#insert paths to high resolution spectra - assume Asimut output with 2 columns (nm and W m-2 sr-1 µm-1), 4 header rows and nm in decreasing order
#if left blank, the OIP values will be used instead for comparison with their SNR report

#use OIP max/mean/min radiance values from the report if filepaths are left blank
day_input_filepath = ""
night_input_filepath = ""




#output filepath (band number will be appended to each filename). Output files contain:
#Wavelength (um), Radiance (W.m-2.sr-1.um-2), SNR (1 pixel), SNR (binned), Radiance error 1 pixel (W.m-2.sr-1.um-2), Radiance error binned (W.m-2.sr-1.um-2)
output_filepath = "Simulations_LIDORT.txt"

#save SNR png? Output filepath will be used
# plot_output = False
plot_output = True

#save output data files?
save_output = True
# save_output = False

#convolve high resolution spectrum to pixel ILS? If False, interpolate input spectrum to pixel central wavelength
#if using OIP max/mean/min radiance values then set this to False
# ils_convolution = False
ils_convolution = True

#band 2b SNRs are consistently 5-7% too small, probably due to the quantum efficiency figure that was approximated from OIP's study
#If you set the scale_2b to True, it increases the Venus signal by 6% to get consistent results with the OIP model.
#increase band 2b radiance by 6% to match OIP calculations?
scale_2b = True


#the max integration time is calculated by the model to fill the detector to 80% of the full well, up to a maximum of 14.4 seconds
#the best detector gain setting is not yet calculated
#auto mode not yet implemented => use OIP's chosen gain
#if choosing a temperature not defined by OIP, the settings for cold temperatures will be used
# auto = False


#scale the Venus input radiance by a factor. 1.0 = no change.
input_radiance_scalar = 1.0


#simulate blaze function: shape can be 'flat' or 'sinc2'
blaze_shape = "flat" #OIP SNR study used a flat grating efficiency
# blaze_shape = "sinc2"

#grating efficiency (at the peak of the blaze if not flat)
blaze_peak = 0.9 #OIP SNR study assumed 90%


#transmittance of spectrometer sections
#not including grating - see blaze_shape and blaze_peak above to set these values independently
transmittance_cold_section = 0.516 #total transmittance of the spectrometer section without grating
transmittance_complete = 0.41 #total transmittance of the instrument without grating




#run Venus and dark frames separately for subtraction on ground?
#not yet implemented
# separate_bg = False


#specify bands to analyse, dayside (d) and/or nightside (n)
bands = [
    ["1", "n"],
    ["2a", "d"],
    ["2a", "n"],
    ["2b", "d"],
    ["2b", "n"],
    ["3", "n"],
    ["4", "d"],
]


#select instrument section temperatures
t_cold_section = 228. #spectrometer section temperature (K)
t_detector_window = 253. #detector window temperature (K)


#warm section contribution through slit
#not yet implemented
# t_warm_section = 253


#dark current variations with temperature not implemented, but thermal background from detector is included
t_detector = 130.  #detector and cold shield temperature (K)


#define resolving powers
theoretical_resolving_power = 11000. #this is the value for the order calculations and spectral calibration of the detector pixels

#in reality the instrument line shape (ILS) is reduced due to aberations etc. Use this value to calculate the pixel ILS
# real_resolving_power = 8500. #8500 = general minimum from OIP optical design report
real_resolving_power = 7000. #7000 = minimum for science


# maximum allowed integration time based on the maximum allowed nightside footprint
max_integration_time = 14.401 #seconds. 14.401 = value from OIP SNR report; results in groundtrack approx 105km


#ADC number of bits
number_of_bits = 14.


#SNR model assumes that the instrument line shape (ILS) of each pixel is a Gaussian of width sigma (like NOMAD). TBD the best approximation
#choose the width of the ILS (in values of sigma) i.e. calculate ILS for range (-ils_gaussian width * sigma) to (+ils_gaussian width * sigma)
ils_gaussian_width = 3.0 #3-sigma = 99.7% of the radiance included in the ILS. Lower number = faster simulation






class snr_model_calc(object):
    
    def __init__(self, bands):
        
        self.bands = bands
        
        self.day_input_filepath = day_input_filepath
        self.night_input_filepath = night_input_filepath
        self.output_filepath = output_filepath
        
        self.n_bits = number_of_bits

        self.n_columns = 384 #number of detector spectral pixels
        self.fno = 3.0 #f/# number
        self.fsr = 144.88 #cm-1
        self.theoretical_resolving_power = theoretical_resolving_power
        self.real_resolving_power = real_resolving_power
        self.slit_width = 1. #pixels
        
        self.blaze_shape = blaze_shape
        self.blaze_peak = blaze_peak


        self.t_cold_section = t_cold_section
        self.t_detector_window = t_detector_window
        self.t_detector = t_detector
        
        self.ils_convolution = ils_convolution
        self.ils_gaussian_width = ils_gaussian_width
        self.scale_2b = scale_2b
        self.input_radiance_scalar = input_radiance_scalar
        
        self.max_integration_time = max_integration_time
        
        self.transmittance_cold_section = transmittance_cold_section
        self.transmittance_complete = transmittance_complete
        
        # self.det_window_tb_scalar = det_window_tb_scalar
        # self.cold_section_tb_scalar = cold_section_tb_scalar

        load_asimut_spectra(self)
        
        if plot_output:
            prepare_plot(self)


    
        #print OIP comparison results
        if self.day_input_filepath == "" and self.night_input_filepath == "":
            print("SNR 1 pixel values for comparison to the OIP SNR study report")
            print("E.g. If t_cold_section=228 then values should be equal to VENSPEC-H-REP-OIP-003 v3_0 (SNR Analysis) table 4-6 page 52")


        for band_ix, (band, daynight) in enumerate(self.bands):
        
            self.band = band
            self.daynight = daynight
    
            spectral_calibration(self, band)
            
            venus(self, band, daynight)

            detector_settings(self, band, daynight)
            detector(self)
            detector_qe_interp(self)


            cold_shield(self)
            detector_window(self)
            cold_section(self)
    
            signal(self, band)
            
            #print OIP comparison results
            if len(self.signal["snr"]) == 3:
                print("Band", band, {"d":"day", "n":"night"}[daynight], " - SNRs for 1 pixel: min R={:.2f}, mean R={:.2f}, max R=={:.2f}".format(*self.signal["snr"]))
            else: #else print min/mean/max SNR of real spectrum for 1 pixel
                print(band, daynight, "SNR 1 px:", np.min(self.signal["snr"]), np.mean(self.signal["snr"]), np.max(self.signal["snr"]))
            
        
        
            if self.output_filepath != "":
                if save_output:
                    # save_generic_output(self, band, daynight)
                    save_cm1_clean_radiance(self, band, daynight)
                    # save_cm1_radiance_with_noise(self, band, daynight)
                    save_cm1_radiance_error(self, band, daynight)
            
                

            if plot_output:
                save_plot(self, band, daynight, band_ix)              




# run calculations for temperature given above
snr = snr_model_calc(bands)
