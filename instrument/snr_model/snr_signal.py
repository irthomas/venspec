# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 11:27:12 2022

@author: iant
"""

import numpy as np
import matplotlib.pyplot as plt

from snr_functions import F_omegaA, F_signal_det, F_thermal_background, F_adc
from snr_ils import detector_ils



def venus(self, band, daynight):
    #venus signal

    
    if daynight == "d":
        if self.day_input_filepath != "":
            
            
            #perform pixel ILS convolution
            if self.ils_convolution:
                venus_wm2um = detector_ils(self.px_um, self.real_px_delta_lambda_um, self.ils_gaussian_width, self.day_um, self.day_Wm2um)
            else:
                venus_wm2um = np.interp(self.px_um, self.day_um, self.day_Wm2um) #interpolate venus radiance to pixel wavelengths

            self.venus_wm2um = venus_wm2um
            
        else:
            #if no filename given, manually select OIP radiances
            #set wavelength to mean of order
            self.px_um = np.zeros(3) + np.mean(self.px_um)
            self.venus_wm2um = np.zeros(3)
            self.venus_wm2um[2] = {"2a":23, "2b":18, "4":176}[band]
            self.venus_wm2um[1] = {"2a":21, "2b":17, "4":165}[band]
            self.venus_wm2um[0] = {"2a":11, "2b":13, "4":126}[band]



    if daynight == "n":
        if self.night_input_filepath != "":
            
            if self.ils_convolution:
                venus_wm2um = detector_ils(self.px_um, self.real_px_delta_lambda_um, self.ils_gaussian_width, self.night_um, self.night_Wm2um)
            else:
                venus_wm2um = np.interp(self.px_um, self.night_um, self.night_Wm2um) #interpolate venus radiance to pixel wavelengths
                
            self.venus_wm2um = venus_wm2um
            
        else:
            #if no filename given, manually select OIP radiances
            #set wavelength to mean of order
            self.px_um = np.zeros(3) + np.mean(self.px_um)
            self.venus_wm2um = np.zeros(3)
            self.venus_wm2um[2] = {"1":0.1051, "2a":0.1507, "2b":0.0451, "3":0.2635}[band]
            self.venus_wm2um[1] = {"1":0.0314, "2a":0.0702, "2b":0.0252, "3":0.0667}[band]
            self.venus_wm2um[0] = {"1":0.0022, "2a":0.0298, "2b":0.0099, "3":0.0133}[band]


    
    
def calc_integration_time(full_well, total_signal_e_per_s):
    
    integration_time = (full_well * 0.8) / (total_signal_e_per_s)
    return integration_time



def F_sinc2(x, a, b): #inverted sinc2
    a = a + 0.0001 #fudge to stop infinity at peak
    return ((np.sin((x - a) / b * 2.783)**2.0) / (((x - a) / b * 2.783)**2.0))



def blaze_transmittance(self):
    
    if self.blaze_shape == "flat":
        F_blaze = np.zeros_like(self.px_um) + self.blaze_peak
    if self.blaze_shape == "sinc2":
        
        #fsr in wavenumbers
        cm1 = 10000. / self.px_um
        cm1_centre = 10000. / np.mean(self.px_um)
        F_blaze = F_sinc2(cm1, cm1_centre, self.fsr) * self.blaze_peak
        
    return F_blaze
    


def signal(self, band):
    signal = {}
    omegaA = F_omegaA(self.fno, self.detector_area_m2)
    
    #apply blaze
    self.transmittance_with_blaze = self.transmittance_complete * blaze_transmittance(self)
    # plt.figure()
    
    # fig, ax_new = plt.subplots()
    # # ax_new.plot(self.px_um, blaze_transmittance(self))
    # ax_new.plot(self.px_um, F_signal_det(self.venus_wm2um, self.px_um, self.px_delta_lambda_m, self.transmittance_with_blaze, omegaA, self.qe_px))
    # ax_new.plot(self.px_um, F_signal_det(self.venus_wm2um, self.px_um, self.px_delta_lambda_m, self.transmittance_complete, omegaA, self.qe_px))
    
    
    signal["venus_es"] = F_signal_det(self.venus_wm2um, self.px_um, self.px_delta_lambda_m, self.transmittance_with_blaze, omegaA, self.qe_px)

    signal["cold_shield_es"] = F_thermal_background(self.cold_shield["b_m"], self.qe_full_range, self.cold_shield["omega"], \
                     self.detector_area_m2, self.nu_full_range_m, self.cold_shield["emis"], self.cold_shield["trans"])
    signal["det_window_es"] = F_thermal_background(self.det_window["b_m"], self.qe_full_range, self.det_window["omega"], \
                    self.detector_area_m2, self.nu_full_range_m, self.det_window["emis"], self.det_window["trans"])
    signal["spectrometer_es"] = F_thermal_background(self.cold_section["b_m"], self.qe_full_range, self.cold_section["omega"], \
                    self.detector_area_m2, self.nu_full_range_m, self.cold_section["emis"], self.cold_section["trans"])
    


    #fudge to make 2b signal work
    if self.scale_2b and band == "2b":
        signal["venus_es"] *= 1.06
        
    


    #TODO: do warm section contribution
    signal["warmsection_es"] = 0.0
    
    # signal["tb_total_es"] = signal["cold_shield_es"] + signal["det_window_es"] + signal["spectrometer_es"] + signal["warmsection_es"]
    
    #TO DO: do this calculation properly rather than analytically
    signal["tb_total_es"] = signal["cold_shield_es"] + signal["det_window_es"] * self.det_window_tb_scalar +\
        signal["spectrometer_es"] * self.cold_section_tb_scalar + signal["warmsection_es"]



    #calculate optimum integration time
    total_signal_e_per_s = np.max(signal["venus_es"]) + self.detector_dc_es + signal["tb_total_es"]
    optimum_it = calc_integration_time(self.full_well, total_signal_e_per_s)
    if optimum_it < self.integration_time:
        print("%s %s %iK: integration time too large: reducing from %0.2f to %0.2f" %(self.band, self.daynight, self.t_cold_section, self.integration_time, optimum_it))
        print(self.band, self.daynight, "max signal:", total_signal_e_per_s, "from", np.max(signal["venus_es"]), self.detector_dc_es, signal["tb_total_es"])
        self.integration_time = optimum_it



    
    
    signal["tb_total_e"] = signal["tb_total_es"] * self.integration_time
    
    # """adc noise"""
    signal["adc_venus"] = F_adc(signal["venus_es"], self.detector_dc_es, signal["tb_total_es"], self.integration_time, self.n_bits)
    signal["adc_dark"] = F_adc(0.0, self.detector_dc_es, signal["tb_total_es"], self.integration_time, self.n_bits) #no signal on dark
    
    
    signal["venus_e"] = signal["venus_es"] * self.integration_time
    signal["noise_e"] = np.sqrt(signal["venus_e"] + signal["adc_venus"]**2. + signal["adc_dark"]**2. + 2. * \
              ((self.detector_dc_es + signal["tb_total_es"]) * self.integration_time + self.readout_noise**2 ))
    
    signal["snr"] = signal["venus_e"] / signal["noise_e"]
    
    signal["snr_binned"] = signal["snr"] * np.sqrt(self.n_rows)
    
    signal["S+TB+DC_e"] = signal["venus_e"] + signal["tb_total_e"] + self.detector_dc_e
    signal["TB+DC_e"] = signal["tb_total_e"] + self.detector_dc_e
    
    signal["S_LSB"] = signal["venus_e"] * 2**self.n_bits / self.full_well
    signal["TB+DC_LSB"] = signal["TB+DC_e"] * 2**self.n_bits / self.full_well
    
    self.signal = signal





# blaze_shape = "flat"
# # blaze_shape = "sinc2"
# blaze_peak = 0.9
# px_um = snr.px_um
# fsr = snr.fsr

# plt.plot(px_um, blaze_transmittance(blaze_shape, blaze_peak, 10000.0 / px_um, fsr))

