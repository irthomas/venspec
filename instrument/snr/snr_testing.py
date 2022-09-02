# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 16:26:52 2022

@author: iant

SNR MODEL COMPARING TO VALUES IN OIP SNR REPORT V3.0
"""
# import os
import numpy as np
import scipy.constants as spc
import matplotlib.pyplot as plt



from instrument.snr.spectrometer_calculations import spectral_cal
from instrument.snr.hr_nadir_spectra import hr_nadir_spectra
from instrument.snr.detector_qe import detector_qe



#scalars to make the thermal background work
det_window_tb_scalar = 0.32
cold_section_tb_scalar = 0.76





def F_planck(nu_m, T): #W/m2/sr/m

    a = 2. * spc.h * spc.c**2.
    b =  spc.h* spc.c /(nu_m * spc.k * T)
    planck = a / ( (nu_m**5.) * (np.exp(b) - 1.) )
    return planck



def F_signal_det(px_Wm2m, px_m, px_fwhm_m, trans, omegaA, qe):
    
    signal = px_Wm2m * px_m * px_fwhm_m * trans * omegaA * qe / (spc.c * spc.h)
    
    return signal



def F_thermal_background(planck, qe, omega, px_area_m, nu_m, emis, trans):
    
    tb = planck * qe * omega * px_area_m * nu_m * emis * trans / (spc.c * spc.h)
    
    # plt.figure()
    # plt.plot(nu_m, tb)
    integrated_tb = np.trapz(tb, x=nu_m) #W/m2/sr
    return integrated_tb



def F_adc(signal, dc, tb, it, n_bits):
    #set signal = 0 for dark frames
    
    adc = ((signal + dc + tb) * it) / (np.sqrt(12.) * 2.**n_bits)
    
    return adc
    
    

def F_solid_angle(fno):
    
    omega = 2. * np.pi * (1. - np.sqrt(1. - (1./(4. * fno**2.))))
    return omega
    


def F_omegaA(fno, px_m2):
    
    omegaA = F_solid_angle(fno) * px_m2
    return omegaA

            
         
            
         
def detector_dict(band_dict, daynight, t):
    
    px_um = np.zeros(3) + np.mean(band_dict["detector_ums"]) #simplified
    px_m = px_um / 1.0e6
    px_fwhm_um = band_dict["delta_lambda_um"]
    px_fwhm_m = px_fwhm_um / 1.0e6
    
    fno = 3.
    n_bits = 14.

    if t["cold_section"] == 273:
        gain = 2
        
        if daynight == "d":
            integration_time = {"2a":0.948, "2b":0.948, "4":0.365}[band]
        if daynight == "n":
            integration_time = {"1":14.401, "2a":14.401, "2b":14.401, "3":14.401}[band]

    elif t["cold_section"] == 253:
        
        if daynight == "d":
            gain = 2
            integration_time = {"2a":0.976, "2b":0.976, "4":0.369}[band]
        if daynight == "n":
            gain = 1
            integration_time = {"1":14.401, "2a":14.401, "2b":14.401, "3":14.401}[band]

    elif t["cold_section"] == 228:
        
        if daynight == "d":
            gain = 2
            integration_time = {"2a":0.982, "2b":0.982, "4":0.370}[band]
        if daynight == "n":
            gain = 1
            integration_time = {"1":14.401, "2a":14.401, "2b":14.401, "3":14.401}[band]

    elif t["cold_section"] == 213:
        
        if daynight == "d":
            gain = 2
            integration_time = {"2a":0.982, "2b":0.982, "4":0.370}[band]
        if daynight == "n":
            gain = 1
            integration_time = {"1":14.401, "2a":14.401, "2b":14.401, "3":14.401}[band]



    readout_noise = {1:125., 2:235.}[gain]
    full_well = {1:340.0e3, 2:1.0e6}[gain]
    detector = {"gain":gain, "readout_noise":readout_noise, "full_well":full_well, "it":integration_time, "fno":fno, "n_bits":n_bits} # electons
            
            
        
        

    
    
    # """signal on detector"""
    # """readout noise (only for gain=2)"""
    # detector["qe"] = 0.75 #simplified from OIP model
    detector["n_px"] = band_dict["n_px"]
    
    detector["px_m"] = px_m
    detector["px_um"] = px_um
    detector["px_fwhm_um"] = px_fwhm_um
    detector["px_fwhm_m"] = px_fwhm_m
    
    detector["size_um"] = 24 #um
    detector["size_m"] = 24 / 1.0e6 #um

    detector["area_um2"] = detector["size_um"] ** 2.
    detector["area_m2"] = detector["size_m"] ** 2.


    detector["dc_nAcm2"] = 0.023 # nA/cm2
    detector["dc_nAum2"] = detector["dc_nAcm2"] * (100.)**2. / (1.0e6)**2
    detector["dc_Aum2"] = detector["dc_nAum2"] * 1.0e-9
    detector["dc_A"] = detector["dc_Aum2"] * detector["area_um2"]
    detector["dc_es"] = 6.242e18 * detector["dc_A"]
    detector["dc_e"] = detector["dc_es"] * detector["it"]


    #define spectral range of detector for thermal background
    detector["nu_full_range_um"] = np.arange(0.8, 2.6, 0.0001)
    detector["nu_full_range_m"] = detector["nu_full_range_um"] * 1.0e-6

    #TO DO: do QE as function of wavelength
    qe_um, qe = detector_qe()
    
    detector["qe_px"] = np.interp(px_um, qe_um, qe)
    detector["qe_full_range"] = np.interp(detector["nu_full_range_um"], qe_um, qe)
    
    return detector



def venus_dict(daynight, band):
    
    #venus signal
    venus = {"b_wm2um":np.zeros(3)}
    #manually select OIP radiances
    if daynight == "d":
        venus["b_wm2um"][2] = {"2a":23, "2b":18, "4":176}[band]
        venus["b_wm2um"][1] = {"2a":21, "2b":17, "4":165}[band]
        venus["b_wm2um"][0] = {"2a":11, "2b":13, "4":126}[band]
    if daynight == "n":
        venus["b_wm2um"][2] = {"1":0.1051, "2a":0.1507, "2b":0.0451, "3":0.2635}[band]
        venus["b_wm2um"][1] = {"1":0.0314, "2a":0.0702, "2b":0.0252, "3":0.0667}[band]
        venus["b_wm2um"][0] = {"1":0.0022, "2a":0.0298, "2b":0.0099, "3":0.0133}[band]

    return venus




def cold_shield_dict(detector, t):
    cold_shield_omega = 2. * np.pi - F_solid_angle(detector["fno"])
    cold_shield = {"trans":1.0, "emis":1.0, "fno":detector["fno"], "omega":cold_shield_omega}
    
    cold_shield["b_m"] = F_planck(detector["nu_full_range_m"], t["cold_shield"])
    
    return cold_shield
    


        
def det_window_dict(detector, t):
    #TODO: replace with table values. At 228K thermal background changes by just 2 e-/s"""
    reflectance = 0.01
    thickness = 0.1 #cm
    abs_coeff = 0.002 #cm-1 #see page 182 for full table https://link.springer.com/content/pdf/10.1007/978-0-387-85695-7.pdf
    
    #note: typo in OIP formula
    #see equation 1 here: https://avs.scitation.org/doi/am-pdf/10.1116/1.4954211#:~:text=For%20light%20impinging%20normally%20on,1%20%E2%80%93%20R2exp(%2D2%CE%B1d)%5D.
    #for thermal purposes, transmission from window to detector = 1
    det_window_trans = ((1. - reflectance)**2 * np.exp(-abs_coeff * thickness))/(1. - np.exp(-abs_coeff * thickness) * reflectance**2)
    det_window_omega = F_solid_angle(detector["fno"])
    det_window = {"trans":1., "emis":1. - det_window_trans, "omega":det_window_omega}
    det_window["b_m"] = F_planck(detector["nu_full_range_m"], t["det_window"])
    
    return det_window
    


        
def spectrometer_dict(detector, t):
    spectrometer_omega = F_solid_angle(detector["fno"])
    # spectrometer = {"trans":1., "emis":1., "omega":spectrometer_omega}
    spectrometer = {"trans":transmittance["cold"], "emis":1., "omega":spectrometer_omega}
    
    spectrometer["b_m"] = F_planck(detector["nu_full_range_m"], t["cold_section"])
    return spectrometer
        



def signal_dict():
    signal = {}
    signal["omegaA"] = F_omegaA(detector["fno"], detector["area_m2"])
    #TODO: replace with ILS calc and blaze
    signal["venus_es"] = F_signal_det(venus["b_wm2um"], detector["px_um"], detector["px_fwhm_m"], transmittance["complete"], signal["omegaA"], detector["qe_px"]) #electrons per second
    signal["cold_shield_es"] = F_thermal_background(cold_shield["b_m"], detector["qe_full_range"], cold_shield["omega"], \
                     detector["area_m2"], detector["nu_full_range_m"], cold_shield["emis"], cold_shield["trans"])
    signal["det_window_es"] = F_thermal_background(det_window["b_m"], detector["qe_full_range"], det_window["omega"], \
                    detector["area_m2"], detector["nu_full_range_m"], det_window["emis"], det_window["trans"])
    signal["spectrometer_es"] = F_thermal_background(spectrometer["b_m"], detector["qe_full_range"], spectrometer["omega"], \
                    detector["area_m2"], detector["nu_full_range_m"], spectrometer["emis"], spectrometer["trans"])
    



    signal["warmsection_es"] = 0.0
    
    # signal["tb_total_es"] = signal["cold_shield_es"] + signal["det_window_es"] + signal["spectrometer_es"] + signal["warmsection_es"]
    
    #TO DO: do this calculation properly rather than analytically
    signal["tb_total_es"] = signal["cold_shield_es"] + signal["det_window_es"] * det_window_tb_scalar +\
        signal["spectrometer_es"] * cold_section_tb_scalar + signal["warmsection_es"]
    # signal["tb_total_es"] = np.exp(0.0811523259 * t["cold_section"] - 11.9073906)
    
    signal["tb_total_e"] = signal["tb_total_es"] * detector["it"]
    
    # """adc noise"""
    signal["adc_venus"] = F_adc(signal["venus_es"], detector["dc_es"], signal["tb_total_es"], detector["it"], detector["n_bits"])
    signal["adc_dark"] = F_adc(0.0, detector["dc_es"], signal["tb_total_es"], detector["it"], detector["n_bits"]) #no signal on dark
    
    
    signal["venus_e"] = signal["venus_es"] * detector["it"]
    signal["noise_e"] = np.sqrt(signal["venus_e"] + signal["adc_venus"]**2. + signal["adc_dark"]**2. + 2. * \
              ((detector["dc_es"] + signal["tb_total_es"]) * detector["it"] + detector["readout_noise"]**2 ))
    
    signal["snr"] = signal["venus_e"] / signal["noise_e"]
    
    signal["snr_binned"] = signal["snr"] * np.sqrt(detector["n_px"])
    
    signal["S+TB+DC_e"] = signal["venus_e"] + signal["tb_total_e"] + detector["dc_e"]
    signal["TB+DC_e"] = signal["tb_total_e"] + detector["dc_e"]
    
    signal["S_LSB"] = signal["venus_e"] * 2**detector["n_bits"] / detector["full_well"]
    signal["TB+DC_LSB"] = signal["TB+DC_e"] * 2**detector["n_bits"] / detector["full_well"]
    
    return signal






"""find best scalars"""
# bands_dict = spectral_cal()
# print("Window, cold, % diff")
# for det_window_tb_scalar in np.arange(0.2, 0.4, 0.02):
#     for cold_section_tb_scalar in np.arange(0.7, 0.9, 0.02):
        
#         diff_percent = 0.0

#         for t_cold in [213., 228., 253., 273.]:
            
#             for t_window in [233., 243., 253., 263., 268., 273.]:
        
#                 band = "1"
            
#                 daynight = "n"
                
#                 band_dict = bands_dict[band]
        
        
#                 #temperatures in Kelvin
#                 t = {"cold_section":t_cold, "det_window":t_window, "detector":130.}
#                 t["cold_shield"] = t["detector"]
                
                
                
#                 transmittance = {"complete":0.37, "cold":0.46}
        
                
                
#                 detector = detector_dict(band_dict, daynight)
#                 venus = venus_dict(daynight, band)
#                 cold_shield = cold_shield_dict(detector, t)
#                 det_window = det_window_dict(detector, t)
#                 spectrometer = spectrometer_dict(detector, t)
#                 signal = signal_dict()
                
#                 d = {
#                     (213, 253):49.,
        
#                     (228, 233):47.,
#                     (228, 243):47.,
#                     (228, 253):46.,
#                     (228, 263):44.5,
#                     (228, 268):44.,
#                     (228, 273):42.5,
        
#                     (253, 253):28.,
        
#                     (273, 253):13.,
#                 }
#                 if (t_cold, t_window) in d.keys():
#                     snr_oip = d[(t_cold, t_window)]
#                     snr = signal["snr"][0]
#                     diff_percent += np.abs(100.*(snr_oip-snr)/snr_oip)
                    
#         if np.abs(diff_percent) < 12.0:
#             print(det_window_tb_scalar, cold_section_tb_scalar, diff_percent)
                    
               
"""output results like in OIP report tables"""
bands = [
    # ["2a", "d"],
    ["2b", "d"],
    # ["4", "d"],
    # ["1", "n"],
    # ["2a", "n"],
    # ["2b", "n"],
    # ["3", "n"],
]



bands_dict = spectral_cal()
# for t_cold in [213., 228., 253., 273.]:
for t_cold in [273.]:
    
    table = np.zeros((17, 7))

    for i, (band, daynight) in enumerate(bands):
            
            band_dict = bands_dict[band]
    

            #temperatures in Kelvin
            t = {"cold_section":t_cold, "det_window":253., "detector":130.}
            t["cold_shield"] = t["detector"]
            transmittance = {"complete":0.37, "cold":0.46}
    
            
            detector = detector_dict(band_dict, daynight, t)
            venus = venus_dict(daynight, band)
            cold_shield = cold_shield_dict(detector, t)
            det_window = det_window_dict(detector, t)
            spectrometer = spectrometer_dict(detector, t)
            signal = signal_dict()
            
            table[:, i] = [
                signal["snr"][0],
                signal["snr"][1],
                signal["snr"][2],

                detector["n_px"],

                signal["snr_binned"][0],
                signal["snr_binned"][1],
                signal["snr_binned"][2],

                detector["it"],
                detector["gain"],
                detector["full_well"],

                signal["S+TB+DC_e"][0],
                signal["S+TB+DC_e"][1],
                signal["S+TB+DC_e"][2],
                
                signal["S_LSB"][0],
                signal["S_LSB"][1],
                signal["S_LSB"][2],
                
                signal["TB+DC_LSB"],

            ]


"""plot QE"""
plt.plot(detector["nu_full_range_um"], detector["qe_full_range"])
plt.grid()
plt.xlabel("Wavelength um")
plt.ylabel("QE")
plt.title("Detector QE estimated from figure 3-8")
plt.savefig("qe.png")