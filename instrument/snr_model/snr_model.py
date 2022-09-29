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
from save_output import prepare_plot, save_plot#, save_generic_output
from save_output_ASIMUT import save_cm1_clean_radiance, save_cm1_radiance_error#, save_cm1_radiance_with_noise
from input_spectra import load_asimut_spectra

#paths to high resolution spectra - assume Asimut output with 2 columns (nm and W m-2 sr-1 µm-1), 4 header rows and nm in decreasing order
#if left blank, the OIP values will be used instead for comparison with their SNR report

### LINUX ###
# day_input_filepath = "Simulations_LIDORT_dayside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"
# night_input_filepath = "Simulations_LIDORT_nightside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"

# bands 2a and 2b only
# day_input_filepath = "nominal_test2_LidortG_AER_Haus_fact050_contXS_3E8_fact2_EMlinelist_specR_040_dayside_SP1_FEN1_rad_forw.dat"
# night_input_filepath = "nominal_test2_LidortG_AER_Haus_fact050_contXS_3E8_fact2_EMlinelist_specR_030_nightside_SP1_FEN1_rad_forw.dat"


### WINDOWS ###

# day_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\Simulations_LIDORT_dayside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"
# night_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\Simulations_LIDORT_nightside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"

# all bands high res - old
# day_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\nominal_test2_LidortG_AER_Haus_fact050_contXS_3E8_fact2_HITRAN2020_FBord_125_final_specR_040_dayside_SP1_FEN1_rad_forw.dat"
# night_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\nominal_test2_LidortG_AER_Haus_fact050_contXS_3E8_fact2_HITRAN2020_FBord_125_final_specR_030_nightside_SP1_FEN1_rad_forw.dat"

# bands 2a and 2b only
# day_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\nominal_test2_LidortG_AER_Haus_fact050_contXS_3E8_fact2_EMlinelist_specR_040_dayside_SP1_FEN1_rad_forw.dat"
# night_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\nominal_test2_LidortG_AER_Haus_fact050_contXS_3E8_fact2_EMlinelist_specR_030_nightside_SP1_FEN1_rad_forw.dat"

# all bands high res - 27 Sep 2022
day_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\nominal_test2_LidortG_AER_Haus_fact050_contXS_3E8_fact2_HITRAN2020_FBord_125_dayside_specR_025_SP1_FEN1_radNC.dat"
night_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\nominal_test2_LidortG_AER_Haus_fact050_contXS_3E8_fact2_HITRAN2020_FBord_125_nightside_specR_025_SP1_FEN1_radNC.dat"


#use OIP max/mean/min radiance values from the report
# day_input_filepath = ""
# night_input_filepath = ""


#output filepath (band number will be appended to each filename). Output files contain:
#Wavelength (um), Radiance (W.m-2.sr-1.um-2), SNR (1 pixel), SNR (binned), Radiance error 1 pixel (W.m-2.sr-1.um-2), Radiance error binned (W.m-2.sr-1.um-2)
output_filepath = "Simulations_LIDORT.txt"

#save SNR png? Output filepath will be used
plot_output = False
# plot_output = True

#save output data files?
# save_output = True
save_output = False

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


#scale the Venus input radiance by a factor i.e. 0.5 means that all input radiance values are divided by 2. 1.0 = no change.
input_radiance_scalar = 0.1


#simulate blaze function?
# blaze_shape = "flat"
blaze_shape = "sinc2"
blaze_peak = 0.9


#run Venus and dark frames separately for subtraction on ground?
#not yet implemented
# separate_bg = False

#specify bands to analyse, dayside (d) and/or nightside (n)
bands = [
    # ["1", "n"],
    # ["2a", "d"],
    # ["2a", "n"],
    # ["2b", "d"],
    # ["2b", "n"],
    ["3", "n"],
    # ["4", "d"],
]


#select instrument section temperatures
t_cold_section = 228. #spectrometer section temperature (K)
t_detector_window = 253. #detector window temperature (K)

#warm section contribution through slit
#not yet implemented
# t_warm_section = 273

#dark current variations with temperature not implemented, but thermal background from detector is included
t_detector = 130.  #detector and cold shield temperature (K)


#define resolving powers
theoretical_resolving_power = 11000. #this is the value for the spectral calibration
#reduced ils due to aberations
real_resolving_power = 8500. #this is the value for the ILS, taking into account aberations etc. from OIP optical design report


#ADC number of bits
number_of_bits = 14.



#scalars to make the thermal background work
det_window_tb_scalar = 0.32
cold_section_tb_scalar = 0.76


# including grating
# transmittance_cold_section = 0.46 #total transmittance of the spectrometer section
# transmittance_complete = 0.37 #total transmittance of the instrument

#not including grating
transmittance_cold_section = 0.516 #total transmittance of the spectrometer section
transmittance_complete = 0.41 #total transmittance of the instrument


#choose the width of the ILS (in values of sigma)
ils_gaussian_width = 3.0 #calculate ILS for range (-ils_gaussian width * scalar) to (+ils_gaussian width * scalar)


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
                print(band, daynight, "SNR 1 px:", self.signal["snr"])
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


#run calculations for temperature given above
snr = snr_model_calc(bands)


#run for different cold section temperatures OIP model
# for t_cold_section in [228., 238., 248., 258., 268.]:
#     snr = snr_model_calc(bands)


#run for different cold section temperatures
# for t_cold_section in [273., 263., 253., 243., 233., 228.]:
#     snr = snr_model_calc(bands)



#run for different ILS widths
# plt.figure()
# for ils_gaussian_width in [3., 5., 10.]:
#     snr = snr_model_calc(bands)
#     plt.plot(snr.px_um, snr.venus_wm2um, label=ils_gaussian_width)
# plt.legend()



#run for different cold section temperatures and different input radiance scalar
# for t_cold_section in [273., 263., 253., 243., 233., 228.]:
#     for input_radiance_scalar in [1.0, 1.0/3.0, 0.1]:
#         snr = snr_model_calc(bands)



#plot HR input spectrum and convolved output for various RPs
#choose 1 band only
# from input_spectra import load_asimut_spectra_cm1
# from save_output_ASIMUT import convert_units_add_noise
# import matplotlib.pyplot as plt

# input_radiance_scalar = 0.1
# plt.figure(figsize=(10, 5), constrained_layout=True)

# for i, real_resolving_power in enumerate([5000., 7000., 8000., 9000., 10000.0, 11000., 15000.]):

#     snr = snr_model_calc(bands)
#     if i == 0:
#         load_asimut_spectra_cm1(snr)
#         plt.plot(snr.night_cm1, snr.night_Wcm2cm1, label="High resolution input spectrum")

#     venus_cm1, venus_wcm2cm1, _, _ = convert_units_add_noise(snr.px_um, snr.venus_wm2um, snr.signal["snr"])
#     plt.plot(venus_cm1, venus_wcm2cm1, label="Convolved spectrum on detector RP=%0.0f" %real_resolving_power)
# plt.legend()            
# plt.grid()
# plt.yscale("log")
# plt.xlim((venus_cm1[0]-1, venus_cm1[-1]+1))
# # plt.ylim((np.min(venus_wcm2cm1[0])/1e3, np.max(venus_wcm2cm1[-1]*1e3)))
# plt.ylim((np.min(venus_wcm2cm1[0])/1e2, np.max(venus_wcm2cm1[-1]*1e1)))
# plt.xlabel("Wavenumber cm-1")
# plt.ylabel("Radiance W/cm2/sr/cm-1")
# plt.title("Band %s %s" %(snr.band, {"d":"day", "n":"night"}[snr.daynight]))
# plt.savefig("resolving_power_band_%s_%s.png" %(snr.band, {"d":"day", "n":"night"}[snr.daynight]))



#plot SNR for various blaze efficiencies and shapes
import matplotlib.pyplot as plt
fig1, ax1 = plt.subplots(figsize=(10, 5), constrained_layout=True)
fig2, ax2 = plt.subplots(figsize=(10, 5), constrained_layout=True)

for blaze_shape in ["flat", "sinc2"]:
    for blaze_peak in [0.8, 0.9]:
        snr = snr_model_calc(bands)
        ax1.plot(snr.px_um, snr.signal["snr"], label="Blaze: %s, peak %0.1f" %(blaze_shape, blaze_peak))
        ax2.plot(snr.px_um, snr.signal["snr_binned"], label="Blaze: %s, peak %0.1f" %(blaze_shape, blaze_peak))

ax1.legend()            
ax1.grid()
ax1.set_xlabel("Wavelength um")
ax1.set_ylabel("SNR 1 pixel")
ax1.set_title("Band %s %s SNR (1 pixel)\nBlaze shape: %s peak: %0.1f" %(snr.band, {"d":"day", "n":"night"}[snr.daynight], blaze_shape, blaze_peak))
fig1.savefig("blaze_snr1px_band_%s_%s.png" %(snr.band, {"d":"day", "n":"night"}[snr.daynight]))
    
ax2.legend()            
ax2.grid()
ax2.set_xlabel("Wavelength um")
ax2.set_ylabel("SNR Binned")
ax2.set_title("Band %s %s SNR (binned)\nBlaze shape: %s peak: %0.1f" %(snr.band, {"d":"day", "n":"night"}[snr.daynight], blaze_shape, blaze_peak))
fig2.savefig("blaze_snrbinned_band_%s_%s.png" %(snr.band, {"d":"day", "n":"night"}[snr.daynight]))
    
