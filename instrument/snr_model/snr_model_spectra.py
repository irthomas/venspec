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
day_input_filepath = r"C:/Users/iant/Nextcloud/VenSpecH.Staging/VenSpecH.BIRA.Staging/SN (Science Note)/ENVIS-VS-VEH-SN-0017-iss0rev0+AER-DaySpectrum-20221116.dat"
night_input_filepath = r"C:/Users/iant/Nextcloud/VenSpecH.Staging/VenSpecH.BIRA.Staging/SN (Science Note)/ENVIS-VS-VEH-SN-0016-iss0rev0+AER-NightSpectrum-20221116.dat"


#use OIP max/mean/min radiance values from the report if filepaths are left blank
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
# scale_2b = True
scale_2b = False


#the max integration time is calculated by the model to fill the detector to 80% of the full well, up to a maximum of 14.4 seconds
#the best detector gain setting is not yet calculated
#auto mode not yet implemented => use OIP's chosen gain
#if choosing a temperature not defined by OIP, the settings for cold temperatures will be used
# auto = False


#scale the Venus input radiance by a factor. 1.0 = no change; 0.1 = closest match of Sev output to previous studies
input_radiance_scalar = 1.0


#simulate blaze function: shape can be 'flat' or 'sinc2'
blaze_shape = "flat" #OIP SNR study used a flat grating efficiency
# blaze_shape = "sinc2"

#grating efficiency (at the peak of the blaze if not flat)
# blaze_peak = 0.9 #OIP SNR study assumed 90%


#transmittance of spectrometer sections
#not including grating - see blaze_shape and blaze_peak above to set these values independently
# transmittance_cold_section = 0.516 #total transmittance of the spectrometer section without grating
# transmittance_complete = 0.41 #total transmittance of the instrument without grating

transmittance_band = {"1":0.339, "2a":0.328, "2b":0.313, "3":0.332, "4":0.306}
transmittance_cold_section = 0.4 #now including the grating, but same for all bands




#run Venus and dark frames separately for subtraction on ground?
#not yet implemented
# separate_bg = False


#specify bands to analyse, dayside (d) and/or nightside (n)
bands = [
    ["1", "n"],
    # ["2a", "d"],
    # ["2a", "n"],
    # ["2b", "d"],
    # ["2b", "n"],
    # ["3", "n"],
    # ["4", "d"],
]


#select instrument section temperatures
t_cold_section = 228. #spectrometer section temperature (K)
t_detector_window = 253. #detector window temperature (K)


#warm section contribution through slit
#not yet implemented
# t_warm_section = 253


#dark current variations with temperature not implemented, but thermal background from detector is included
t_detector = 120.  #detector and cold shield temperature (K)


#define resolving powers
theoretical_resolving_power = 11000. #this is the value for the order calculations and spectral calibration of the detector pixels

#in reality the instrument line shape (ILS) is reduced due to aberations etc. Use this value to calculate the pixel ILS
# real_resolving_power = 8500. #8500 = general minimum from OIP optical design report
real_resolving_power = 9000. #7000 = minimum for science


# maximum allowed integration time based on the maximum allowed nightside footprint
max_integration_time = 14.3 #seconds. 14.401 = value from OIP SNR report; results in groundtrack approx 105km


#ADC number of bits
number_of_bits = 14.


#SNR model assumes that the instrument line shape (ILS) of each pixel is a Gaussian of width sigma (like NOMAD). TBD the best approximation
#choose the width of the ILS (in values of sigma) i.e. calculate ILS for range (-ils_gaussian width * sigma) to (+ils_gaussian width * sigma)
ils_gaussian_width = 5.0 #3-sigma = 99.7% of the radiance included in the ILS. Lower number = faster simulation


# when the pixel wavelengths have been calculated, shift each wavelength by this amount
spectral_shift = 0.0 #um shift from the nominal wavelength assignment





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
        # self.blaze_peak = blaze_peak
        
        self.spectral_shift = spectral_shift


        self.t_cold_section = t_cold_section
        self.t_detector_window = t_detector_window
        self.t_detector = t_detector
        
        self.ils_convolution = ils_convolution
        self.ils_gaussian_width = ils_gaussian_width
        self.scale_2b = scale_2b
        self.input_radiance_scalar = input_radiance_scalar
        
        self.max_integration_time = max_integration_time
        
        self.transmittance_cold_section = transmittance_cold_section
        # self.transmittance_complete = transmittance_complete
        self.transmittance_band = transmittance_band
        
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


#run for different cold section temperatures OIP model
# for t_cold_section in [228., 238., 248., 258., 268.]:
#     snr = snr_model_calc(bands)


#run for different cold section temperatures
# for t_cold_section in [273., 263., 253., 243., 233., 228.]:
#     snr = snr_model_calc(bands)



#for Justin
#run for different ILS widths
# plt.figure()
# for ils_gaussian_width in [3., 5., 10.]:
#     snr = snr_model_calc(bands)
#     plt.plot(snr.px_um, snr.venus_wm2um, label=ils_gaussian_width)
# plt.legend()


#for Severine
#run for different cold section temperatures and different input radiance scalar
# for t_cold_section in [273., 263., 253., 243., 233., 228.]:
#     for input_radiance_scalar in [1.0, 1.0/3.0, 0.1]:
#         snr = snr_model_calc(bands)

#for Severine
#run for different footprint sizes and different input radiance scalar
# for max_integration_time in [14.4, 28.8, 43.2, 57.6, 72.0]:
#     for input_radiance_scalar in [0.1, 1.0]:
#         snr = snr_model_calc(bands)



# plot HR input spectrum and convolve output for various RPs
# choose 1 band only
# from input_spectra import load_asimut_spectra_cm1
# from save_output_ASIMUT import convert_units_add_noise
# import matplotlib.pyplot as plt

# input_radiance_scalar = 0.1
# plt.figure(figsize=(10, 5), constrained_layout=True)

# plot_output = False
# save_output = False

# for i, real_resolving_power in enumerate([2000., 7000., 8400., 11000.]):

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
# import matplotlib.pyplot as plt
# fig1, ax1 = plt.subplots(figsize=(10, 5), constrained_layout=True)
# fig2, ax2 = plt.subplots(figsize=(10, 5), constrained_layout=True)

# for blaze_shape in ["flat", "sinc2"]:
#     for blaze_peak in [0.8, 0.9]:
#         snr = snr_model_calc(bands)
#         ax1.plot(snr.px_um, snr.signal["snr"], label="Blaze: %s, peak %0.1f" %(blaze_shape, blaze_peak))
#         ax2.plot(snr.px_um, snr.signal["snr_binned"], label="Blaze: %s, peak %0.1f" %(blaze_shape, blaze_peak))

# ax1.legend()            
# ax1.grid()
# ax1.set_xlabel("Wavelength um")
# ax1.set_ylabel("SNR 1 pixel")
# ax1.set_title("Band %s %s SNR (1 pixel)\nBlaze shape: %s peak: %0.1f" %(snr.band, {"d":"day", "n":"night"}[snr.daynight], blaze_shape, blaze_peak))
# fig1.savefig("blaze_snr1px_band_%s_%s.png" %(snr.band, {"d":"day", "n":"night"}[snr.daynight]))
    
# ax2.legend()            
# ax2.grid()
# ax2.set_xlabel("Wavelength um")
# ax2.set_ylabel("SNR Binned")
# ax2.set_title("Band %s %s SNR (binned)\nBlaze shape: %s peak: %0.1f" %(snr.band, {"d":"day", "n":"night"}[snr.daynight], blaze_shape, blaze_peak))
# fig2.savefig("blaze_snrbinned_band_%s_%s.png" %(snr.band, {"d":"day", "n":"night"}[snr.daynight]))
    


# plot HR input spectrum vs planck functions of different temperatures
# import matplotlib.pyplot as plt
# from input_spectra import load_asimut_file_cm1

# cm1, wcm2cm1 = load_asimut_file_cm1(night_input_filepath)

# um = 10000.0 / cm1
# wm2um = wcm2cm1 * 1.0e8 / um**2

# fig1, ax1 = plt.subplots(figsize=(10, 5), constrained_layout=True)
# fig2, ax2 = plt.subplots(figsize=(10, 5), constrained_layout=True)
# ax1.plot(cm1, wcm2cm1, label="High resolution input spectrum")
# ax2.plot(um, wm2um, label="High resolution input spectrum")

# c1 = 1.191e-5
# c2 = 1.43877

# for t in [500, 600, 700]:
#     lcm1 = ((c1 * cm1**3) / (np.exp(c2 * cm1 / t)-1))/1.0e7
#     lum = lcm1 * 1.0e8 / um**2
#     ax1.plot(cm1, lcm1, label="%iK Planck" %t)
#     ax2.plot(um, lum, label="%iK Planck" %t)
# # ax1.set_yscale("log")
# ax1.legend()
# # ax2.set_yscale("log")
# ax2.legend()
# ax2.set_xlim(1.16, 1.18)
# ax2.set_ylim(0, 0.5)




# plot HR input spectrum and convolve output for various RPs
# choose 1 band only
# from input_spectra import load_asimut_spectra_cm1
# from save_output_ASIMUT import convert_units_add_noise
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages

# input_radiance_scalar = 0.1

# plot_output = False
# save_output = False

# shift = 0.0001063522949549256 #band 1 spectral sampling
# real_resolving_power = 11000.

# with PdfPages("spectral_shift.pdf") as pdf: #open pdf


#     for i, spectral_shift_scalar in enumerate(np.arange(10.0) /10.0):
    
#         fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
#         spectral_shift = shift * spectral_shift_scalar
        
        
#         snr = snr_model_calc(bands)
#         load_asimut_spectra_cm1(snr)
#         ax.plot(snr.night_cm1, snr.night_Wcm2cm1, label="High resolution input spectrum")
    
#         venus_cm1, venus_wcm2cm1, _, _ = convert_units_add_noise(snr.px_um, snr.venus_wm2um, snr.signal["snr"])
#         ax.plot(venus_cm1, venus_wcm2cm1)
    
#         # plt.legend()            
#         ax.grid()
#         ax.set_yscale("log")
#         ax.set_xlim((venus_cm1[0]-1, venus_cm1[-1]+1))
#         # plt.ylim((np.min(venus_wcm2cm1[0])/1e3, np.max(venus_wcm2cm1[-1]*1e3)))
#         ax.set_ylim((np.min(venus_wcm2cm1)/2, np.max(venus_wcm2cm1)*2))
#         ax.set_xlabel("Wavenumber cm-1")
#         ax.set_ylabel("Radiance W/cm2/sr/cm-1")
#         ax.set_title("Band %s %s\nPixel centre wavelengths shifted by %0.06fum" %(snr.band, {"d":"day", "n":"night"}[snr.daynight], spectral_shift))

#         pdf.savefig(fig)
#         plt.close(fig)
