# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 10:12:56 2022

@author: iant
"""

import numpy as np

from spectral_calibration import spectral_calibration
from snr_detector import detector_settings, detector, detector_qe_interp
from snr_signal import venus, signal
from snr_thermal_background import cold_shield, detector_window, cold_section
from save_output import save_file, prepare_plot, save_plot
from input_spectra import load_asimut_spectra

#paths to high resolution spectra - assume Asimut output with 2 columns (nm and W m-2 sr-1 µm-1), 4 header rows and nm in decreasing order
#if left blank, the OIP values will be used instead for comparison with their SNR report
# day_input_filepath = "Simulations_LIDORT_dayside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"
# night_input_filepath = "Simulations_LIDORT_nightside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"

#windows 
day_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\Simulations_LIDORT_dayside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"
night_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\Simulations_LIDORT_nightside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"

#use OIP max/mean/min radiance values from the report
# day_input_filepath = ""
# night_input_filepath = ""


#output filepath (band number will be appended to each filename). Output files contain:
#Wavelength (um), Radiance (W.m-2.sr-1.um-2), SNR (1 pixel), SNR (binned), Radiance error 1 pixel (W.m-2.sr-1.um-2), Radiance error binned (W.m-2.sr-1.um-2)
output_filepath = "Simulations_LIDORT.txt"

#save SNR png? Output filepath will be used
plot_output = True

#convolve high resolution spectrum to pixel ILS? If False, take interpolate input spectrum to pixel central wavelength
#ILS convolution not yet implemented
# ils_convolution = False
# ils_convolution = True

#increase band 2b radiance by 6% to match OIP calculations?
scale_2b = True


#automatically calculate best integration time and gain settings based on temperature inputs?
#auto mode not yet implemented => use OIP gain and integration times
#if choosing a temperature not defined by OIP, the settings for cold temperatures will be used
# auto = False


#simulate blaze function?
#not yet implemented
# blaze = False

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
# t_warm_section = 273

#dark current variations with temperature not implemented, but thermal background from detector is included
t_detector = 130.  #detector and cold shield temperature (K)


#ADC number of bits
number_of_bits = 14.



#scalars to make the thermal background work
det_window_tb_scalar = 0.32
cold_section_tb_scalar = 0.76


transmittance_cold_section = 0.46 #total transmittance of the spectrometer section
transmittance_complete = 0.37 #total transmittance of the instrument



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
        self.resolving_power = 11000.
        self.slit_width = 1. #pixels


        self.t_cold_section = t_cold_section
        self.t_detector_window = t_detector_window
        self.t_detector = t_detector
        self.scale_2b = scale_2b
        
        self.transmittance_cold_section = transmittance_cold_section
        self.transmittance_complete = transmittance_complete
        
        self.det_window_tb_scalar = det_window_tb_scalar
        self.cold_section_tb_scalar = cold_section_tb_scalar

        load_asimut_spectra(self)
        
        if plot_output:
            prepare_plot(self)


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
                print(band, daynight, self.signal["snr"])
            else: #else print min/mean/max SNR of real spectrum
                print(band, daynight, np.min(self.signal["snr"]), np.mean(self.signal["snr"]), np.max(self.signal["snr"]))
            
        
        
            if self.output_filepath != "":
                save_file(self, band, daynight)
            
                

            if plot_output:
                save_plot(self, band, daynight, band_ix)              


#run calculations
snr = snr_model_calc(bands)

#run for different cold section temperatures
# for t_cold_section in [228., 238., 248., 258., 268.]:
#     snr = snr_model_calc(bands)
