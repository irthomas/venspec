# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 09:48:14 2023

@author: iant

MODEL IN PYTHON TO COMPARE TO OIP MODEL
"""


import numpy as np
import scipy.constants as spc
import os
import shutil
import json
import subprocess
import time
from datetime import datetime
import matplotlib.pyplot as plt

from instrument.snr_oip.config import BAND_DATA, MODEL_DATA, MODEL_RESULTS, MODEL_EXE, MODEL_ROOT, MODEL_OUTPUT, GAIN_DATA


bands = BAND_DATA.keys()

# temps_detector = np.arange(118., 119., 1.0) #K
# temps_warmsec = np.arange(273.15, 278.15, 5.0)
# temps_coldsec = np.arange(228.15, 233.15, 5.0) #K
# gains = ["low", "high"]


# for temp_detector in temps_detector:
#     for temp_warmsec in temps_warmsec:
#         for temp_coldsec in temps_coldsec:

temp_detector = 118.
temp_warmsec = 273.15
temp_coldsec = 228.15

#remove existing conf files
for band in BAND_DATA.keys():
    if os.path.exists(BAND_DATA[band]["conf"]):
        # print("Deleting %s config file" %band)
        os.remove(BAND_DATA[band]["conf"])

for band in bands:
    
    #read in config file
    conf_path = BAND_DATA[band]["conf_temp"]
    with open(conf_path, "r") as f:
        config = json.load(f)

    #make changes to conf file
    config["TempOffset"] = 0.0 #don't vary by temp
    config["TempDet"] = temp_detector #don't vary by temp
    config["TempIll"] = temp_warmsec #don't vary by temp
    config["TempSpec"] = temp_coldsec #don't vary by temp
    
    
        
    #write config file to dir to run snr model
    # print("Making band %s config file" %band)
    conf_run = BAND_DATA[band]["conf"]
    with open(conf_run, "w") as f:
        json_str = json.dumps(config, indent=4)
        f.write(json_str)
    


""""read in config data"""
beta_angles = np.loadtxt(os.path.join(MODEL_DATA, "beta.dat"), skiprows=1)
alpha_angle = config["AlphaAngle"]
grooves = config["Grooves"]
centre_order = config["Order"]
fno = config["FN"]
pixel_size = config["PixelPitch"]
full_well = config["Capacity"]
n_bits = config["ADC"]
readout_noise = config["RON"]
slf = config["SRF"]

trans_warm = config["TransmittanceIll"]



detector_area_m2 = pixel_size**2
order = centre_order


#spectral cal
#mλ= d (sinα + sinβ) 
# λ= d (sinα + sinβ)/m
px_nm = grooves * (np.sin(alpha_angle) + np.sin(beta_angles[:, 1])) / order * 1.0e9 #nm


"""get cold section / grating efficiency"""
trans_coldsec_path = BAND_DATA[band]["data"]
trans_coldsec_nm_all, trans_coldsec_all = np.loadtxt(trans_coldsec_path, unpack=True)
#interpolate to pixel wavelengths
trans_coldsec = np.interp(px_nm, trans_coldsec_nm_all, trans_coldsec_all)




"""get detector QE"""
qe_filename = os.path.split(config["QE"])[-1]
qe_nm_all, qe_all = np.loadtxt(os.path.join(MODEL_DATA, qe_filename), skiprows=1, unpack=True)
#interpolate to pixel wavelengths
qe = np.interp(px_nm, qe_nm_all, qe_all)
qe_full = np.interp(trans_coldsec_nm_all, qe_nm_all, qe_all)


"""get detector dc"""
dc_filename = os.path.split(config["DarkCurrent"])[-1]
dc_t_all, dc_all = np.loadtxt(os.path.join(MODEL_DATA, dc_filename), skiprows=3, unpack=True)
#interpolate to detector temperature
dc_nAcm2 = np.interp(temp_detector, dc_t_all, dc_all)

detector_dc_es = dc_nAcm2 * 36000.0




"""get venus spectrum"""
venus_filename = os.path.split(config["Radiance"])[-1]
venus_nm, venus_spectrum_Wm2sr1um1 = np.loadtxt(os.path.join(MODEL_DATA, venus_filename), skiprows=6, unpack=True)
#reverse order for searching
venus_nm = np.append(venus_nm[::-1], 2701.0)
venus_spectrum_Wm2sr1um1 = np.append(venus_spectrum_Wm2sr1um1[::-1], 0.0)
#simple interpolation to pixel wavelengths
venus_interp = np.interp(px_nm, venus_nm, venus_spectrum_Wm2sr1um1)


#TODO: check
nu_full_range_m = trans_coldsec_nm_all/1.0e9
integration_time = 14.126
n_rows = 200


transmittance_complete = trans_coldsec * trans_warm

px_sampling = np.diff(px_nm)
px_sampling = np.append(px_sampling, px_sampling[-1])



def oip_isrf_conv(px_nm, px_sampling, venus_nm, venus_spectrum):
    """convolve venus spectrum to local average i.e. pixel centre nm +- 1/3 delta_nm"""
    venus_conv = np.zeros_like(px_nm)
    
    #loop through pixels
    for i in np.arange(len(venus_conv)):
        #find nm bounds
        nm_start = px_nm[i] - px_sampling[i] / 3.0
        nm_end = px_nm[i] + px_sampling[i] / 3.0
    
        ixs = np.searchsorted(venus_nm, [nm_start, nm_end])
        venus_conv[i] = np.mean(venus_spectrum[ixs[0]:ixs[1]])

    return venus_conv


#get venus spectrum
venus = oip_isrf_conv(px_nm, px_sampling, venus_nm, venus_spectrum_Wm2sr1um1)

# plt.plot(venus_nm, venus_spectrum_Wm2sr1um1)
# plt.plot(px_nm, venus_conv)
# stop()

def F_solid_angle(fno):
    omega = 2. * np.pi * (1. - np.sqrt(1. - (1./(4. * fno**2.))))
    return omega

def F_omegaA(fno, px_m2):
    omegaA = F_solid_angle(fno) * px_m2
    return omegaA

def F_signal_det(px_Wm2m, px_m, px_fwhm_m, trans, omegaA, qe):
    signal = px_Wm2m * px_m * px_fwhm_m * trans * omegaA * qe / (spc.c * spc.h)
    return signal

def F_adc(full_well, n_bits):
    adc = full_well / (np.sqrt(12.) * 2.**n_bits)
    return adc



def F_signal_det2(px_Wm2m, px_m, px_delta_m, trans, qe, fno, px_m2):
    """𝑆(𝐿)=𝐿∙𝜆ℎ𝑐∙𝑇(𝜆)∙𝑄𝐸(𝜆)∙𝑝2∙2𝜋(1−√1−14∙𝐹𝑁2)∙𝜆𝑅"""
    signal = px_Wm2m * px_m * trans * qe * px_m2 * 2. * np.pi * (1. - np.sqrt(1. - (1./(4. * fno**2.)))) * px_delta_m  / (spc.c * spc.h)
    return signal

def F_planck(nu_m, T): #W/m2/sr/m
    a = 2. * spc.h * spc.c**2.
    b =  spc.h* spc.c /(nu_m * spc.k * T)
    planck = a / ( (nu_m**5.) * (np.exp(b) - 1.) )
    return planck

def F_thermal_background(planck, qe, omega, px_m2, nu_m, emis, trans):
    tb = emis * planck * qe * omega * px_m2 * nu_m * trans / (spc.c * spc.h)
    # plt.figure()
    # plt.plot(nu_m, tb)
    integrated_tb = np.trapz(tb, x=nu_m) #W/m2/sr
    return integrated_tb


def F_thermal_background2(qe, omega, px_m2, emis, trans, nu_full_range_m, temp):
    
    integral = (2 * spc.c / nu_full_range_m**4.0) * (1.0 / (np.exp((spc.h * spc.c)/(spc.k * temp * nu_full_range_m)) - 1.0)) * trans * qe
    
    integrated = np.trapz(integral, x=nu_full_range_m) #W/m2/sr
    
    tb = emis * px_m2 * omega * integrated
    # plt.figure()
    # plt.plot(nu_m, tb)
    return tb






def cold_shield_d(fno, t_detector, nu_full_range_m):
    cold_shield_omega = 2. * np.pi - F_solid_angle(fno)
    cold_shield = {"trans":1.0, "emis":1.0, "fno":fno, "omega":cold_shield_omega}
    cold_shield["b_m"] = F_planck(nu_full_range_m, t_detector)
    return cold_shield
    


        
def detector_window_d(fno, t_detector_window, nu_full_range_m):
    #TODO: replace with table values. At 228K thermal background changes by just 2 e-/s"""
    reflectance = 0.01
    thickness = 0.1 #cm
    abs_coeff = 0.002 #cm-1 #see page 182 for full table https://link.springer.com/content/pdf/10.1007/978-0-387-85695-7.pdf
    
    #note: typo in OIP formula
    #see equation 1 here: https://avs.scitation.org/doi/am-pdf/10.1116/1.4954211#:~:text=For%20light%20impinging%20normally%20on,1%20%E2%80%93%20R2exp(%2D2%CE%B1d)%5D.
    #for thermal purposes, transmission from window to detector = 1
    det_window_trans = ((1. - reflectance)**2 * np.exp(-abs_coeff * thickness))/(1. - np.exp(-abs_coeff * thickness) * reflectance**2)
    det_window_omega = F_solid_angle(fno)
    det_window = {"trans":1., "emis":1. - det_window_trans, "omega":det_window_omega}
    det_window["b_m"] = F_planck(nu_full_range_m, t_detector_window)
    return det_window


        
def cold_section_d(fno, t_cold_section, nu_full_range_m, transmittance_cold_section):
    cold_section_omega = F_solid_angle(fno)
    # spectrometer = {"trans":1., "emis":1., "omega":cold_section_omega}
    cold_section = {"trans":transmittance_cold_section, "emis":1., "omega":cold_section_omega}
    cold_section["b_m"] = F_planck(nu_full_range_m, t_cold_section)
    return cold_section        




"""make dicts for certain parts"""

cold_shield = cold_shield_d(fno, temp_detector, nu_full_range_m)
detector_window = detector_window_d(fno, temp_warmsec, nu_full_range_m)
# cold_section = cold_section_d(fno, temp_coldsec, nu_full_range_m, trans_coldsec_all)

#TODO: check transmittance
cold_section = cold_section_d(fno, temp_coldsec, nu_full_range_m, 0.4)



"""do signal calculations"""

signal = {}
omegaA = F_omegaA(fno, detector_area_m2)


#TODO: check this
px_delta_nm = np.diff(px_nm)
px_delta_nm = np.append(px_delta_nm, px_delta_nm[-1])


# these give the same values
# signal_venus_es = F_signal_det(venus*1.0e6, px_nm/1.0e9, px_delta_nm/1.0e9, transmittance_complete, omegaA, qe)
signal_venus_es = F_signal_det2(venus*1.0e6, px_nm/1.0e9, px_delta_nm/1.0e9, transmittance_complete, qe, fno, detector_area_m2)



# assuming cold shield is 2pi - solid angle, these give the same values
# signal_cold_shield_es = F_thermal_background(cold_shield["b_m"], qe_full, cold_shield["omega"], \
#                   detector_area_m2, nu_full_range_m, cold_shield["emis"], cold_shield["trans"])
cold_shield_omega = 2. * np.pi - F_solid_angle(fno)
signal_cold_shield_es = F_thermal_background2(qe_full, cold_shield_omega, detector_area_m2, \
                                              cold_shield["emis"], cold_shield["trans"], nu_full_range_m, temp_detector)

    
    
signal_det_window_es = F_thermal_background(detector_window["b_m"], qe_full, detector_window["omega"], \
                detector_area_m2, nu_full_range_m, detector_window["emis"], detector_window["trans"])

𝜀𝑤𝑖𝑛𝑑𝑜𝑤= 1−𝑇𝑤𝑖𝑛𝑑𝑜𝑤=1− (1−𝑅)2∙𝑒−𝜎∙𝑑1−𝑅2∙𝑒−2𝜎∙𝑑≈1−2𝑅
    
    
signal_spectrometer_es = F_thermal_background(cold_section["b_m"], qe_full, cold_section["omega"], \
            detector_area_m2, nu_full_range_m, cold_section["emis"], cold_section["trans"])
#𝐹𝑠𝑙𝑖𝑡(𝜆)= 𝑇𝑠𝑝𝑒𝑐𝑡𝑟𝑜𝑚𝑒𝑡𝑒𝑟∙Δ𝜆𝐵𝑊𝑜𝑟𝑑𝑒𝑟(𝜆)=𝑇𝑠𝑝𝑒𝑐𝑡𝑟𝑜𝑚𝑒𝑡𝑒𝑟∙𝜆𝑅⁄𝐹𝑆𝑅∙𝜆2=𝑇𝑠𝑝𝑒𝑐𝑡𝑟𝑜𝑚𝑒𝑡𝑒𝑟𝑅∙𝐹𝑆𝑅∙𝜆

#TODO: do warm section contribution
signal_warmsection_es = 0.0




signal_tb_total_es = signal_cold_shield_es + signal_det_window_es + signal_spectrometer_es + signal_warmsection_es

print("signal_cold_shield_es", signal_cold_shield_es)
print("signal_det_window_es", signal_det_window_es)
print("signal_spectrometer_es", signal_spectrometer_es)
print("signal_warmsection_es", signal_warmsection_es)
print("signal_tb_total_es", signal_tb_total_es)

print("signal_venus_es", signal_venus_es[192])
print("detector_dc_es", detector_dc_es)
    
signal_tb_total_e = signal_tb_total_es * integration_time

signal_adc = F_adc(full_well, n_bits)



# snr_1px = (signal_venus_es * integration_time) / \
#     np.sqrt((signal_venus_es + 2.0 * detector_dc_es + 2.0 * signal_tb_total_es) * integration_time + \
#             2.0 * readout_noise**2.0 + 2.0 * signal_adc**2.0)

snr_1px = (signal_venus_es * integration_time) / \
    np.sqrt(((slf + signal_venus_es) + 2.0 * detector_dc_es + 2.0 * signal_tb_total_es) * integration_time + \
            2.0 * readout_noise**2.0 + 2.0 * signal_adc**2.0)



snr_binned = snr_1px * np.sqrt(n_rows)


# plt.plot(cold_shield["b_m"])
# plt.plot(detector_window["b_m"])
# plt.plot(cold_section["b_m"])



#read in OIP equivalent
#get output directory to reflect temperatures
dir_name = "out_d%0.1f_w%0.1f_c%0.1f" %(temp_detector, temp_warmsec, temp_coldsec)
dir_path = os.path.join(MODEL_RESULTS, dir_name)

#read in output file
res_path = os.path.join(dir_path, BAND_DATA[band]["snr"])
with open(res_path, "r") as f:
    lines = f.readlines()
    
oip_snr = {"px":[], "nm":[], "d_nm":[], "bol_sat":[], "bol_snr":[], "eol_sat":[], "eol_snr":[]}
for line in lines:
    if line[0] == "#":
        continue
    else:
        split = line.split("\t")
        for i, key in enumerate(oip_snr.keys()):
            oip_snr[key].append(float(split[i]))
            
for key in oip_snr.keys():
    oip_snr[key] = np.asarray(oip_snr[key])

oip_snr["bol_snr_binned"] = oip_snr["bol_snr"] * np.sqrt(n_rows)

plt.plot(px_nm, snr_binned)
plt.plot(oip_snr["nm"], oip_snr["bol_snr_binned"])

# plt.plot(oip_snr["nm"], 2* oip_snr["bol_snr_binned"] / snr_binned)