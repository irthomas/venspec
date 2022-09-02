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

import streamlit as st




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



name_dict = {
"gain":"gain setting",
"readout_noise":"readout noise (e-)",
"full_well":"full well capacity (e-)",
"qe":"quantum efficiency (simplified)",
"n_px":"number of binned rows",
"it":"integration time (s)",
"fno":"f/#",
"n_bits":"ADC number of bits",

"px_um":"pixel wavelengths (um)",
"px_fwhm_m":"pixel FWHM (m)",
"px_fwhm_um":"pixel FWHM (um)",


"dc_nAcm2":"dark current (nAmps per cm2)",
"size_um":"pixel size (um)",
"size_m":"pixel size (m)",
"dc_nAum2":"dark current (nAmps per um2)",
"dc_Aum2":"dark current (A per um2)",
"area_um2":"pixel area (um2)",
"area_m2":"pixel area (m2)",
"dc_A":"dark current (A)",
"dc_es":"dark current (e- per s)",
"dc_e":"dark current (electrons)",
"t":"temperature (K)",
"trans":"transmission from surface",
"emis":"surface emissivity",
"omega":"solid angle (sr)",


"b_wm2um":"radiance (W/m2/sr/um)",

"cold_section":"spectrometer (cold) section (K)",
"warm_section":"warm baseplate section (K)",
"det_window":"detector window (K)",
"detector":"detector (K)",
"cold_shield":"cold shield around detector (K)",
# "b_m":"planck function",

"qe_px":"detector QE for each pixel",
"qe_full_range":"detector QE full wavelength range",

"omegaA":"etendue (sr m2)",
"venus_es":"Venus (e- per s)",
"cold_shield_es":"cold shield thermal background (e- per s)",
"det_window_es":"detector window thermal background (e- per s)",
"spectrometer_es":"spectrometer (cold) section thermal background (e- per s)",
"warmsection_es":"warm baseplate thermal background (not yet implemented) (e- per s)",
"tb_total_es":"total thermal background (e- per s)",
"tb_total_e":"total thermal background (e-)",
"adc_venus":"ADC noise Venus + background (e-)",
"adc_dark":"ADC noise background (e-)",
"venus_e":"Venus (e-)",
"noise_e":"Noise (e-)",
"snr":"SNR",
"TB+DC_e":"TB+DC (e-)",
}


plot_dict = {
"b_wm2um":{"title":"radiance", "ylabel":"Radiance W/m2/sr/um", "x":"px_um"},
"b_m":{"title":"Planck function", "ylabel":"Radiance W/m2/sr/m", "ylog":True, "x":"nu_full_range_um"},
"snr":{"title":"SNR", "ylabel":"SNR", "x":"px_um"},
"snr_binned":{"title":"SNR detector rows binned", "ylabel":"SNR binned", "x":"px_um"},
"venus_es":{"title":"Venus signal electrons per second", "ylabel":"electrons per s", "x":"px_um"},
"adc_venus":{"title":"ADC noise Venus + background (e-)", "ylabel":"electrons", "x":"px_um"},
"venus_e":{"title":"Venus signal on detector (e-)", "ylabel":"electrons", "x":"px_um"},
"noise_e":{"title":"Noise signal on detector (e-)", "ylabel":"electrons", "x":"px_um"},
"qe_px":{"title":"QE on each detector pixel", "ylabel":"QE", "x":"px_um"},
"qe_full_range":{"title":"QE across full detector range", "ylabel":"QE", "x":"nu_full_range_um"},
"S+TB+DC_e":{"title":"S+TB+DC (e-)", "ylabel":"electrons", "x":"px_um"},
}


hide_dict = {
# "b_m":[],
"nu_full_range_m":[],
"nu_full_range_um":[],
"px_m":[],
}


def write_variables(d, d_name, detector):
    
    st.subheader(f"{d_name} variables")
    for key in d.keys():
        
        if key in hide_dict:
            continue
        
        if key in plot_dict:

            fig1, ax1 = plt.subplots()
            
            if "x" in plot_dict[key]:
                if plot_dict[key]["x"] == "px_um":
                    x = detector["px_um"]
                elif plot_dict[key]["x"] == "nu_full_range_um":
                    x = detector["nu_full_range_um"]
                ax1.plot(x, d[key])
            else:
                ax1.plot(d[key])

            
            if "ylog" in plot_dict[key]:
                if plot_dict[key]["ylog"]:
                    ax1.set_yscale("log")
            title = plot_dict[key]["title"]
            fig1.suptitle(f"{d_name} {title}")
            # ax1.set_xlabel("Wavelength (um)")
            ax1.set_ylabel(plot_dict[key]["ylabel"])
            
            st.pyplot(fig=fig1)
            
        else:
        
            if key in name_dict:
                text = name_dict[key]
                st.write(f"{d_name} {text}: {d[key]}")
            else:
                st.write(f"{d_name} {key}: {d[key]}")
         
            
         
def t_dict():
    #temperatures in Kelvin
    t = {"cold_section":228., "warm_section":273., "det_window":253., "detector":130.}
    t["cold_shield"] = t["detector"]
    return t


def transmittance_dict():
    

    transmittance = {"complete":0.37, "warm":0.82, "cold":0.46, "detector":0.98}
    return transmittance


            
         
def detector_dict(band_dict, daynight):
    
    px_um = band_dict["detector_ums"]
    px_m = px_um / 1.0e6
    px_fwhm_um = band_dict["delta_lambda_um"]
    px_fwhm_m = px_fwhm_um / 1.0e6
    
    fno = 3.
    n_bits = 14.


    if daynight == "d": #gain 2
        integration_time = 0.982
        
        if band == "4":
            integration_time = 0.37
        
        detector = {"gain":2, "readout_noise":235., "full_well":1.0e6, "it":integration_time, "fno":fno, "n_bits":n_bits} # electons

    if daynight == "n": #gain 1
        integration_time = 14.4

        detector = {"gain":1, "readout_noise":125., "full_well":340000., "it":integration_time, "fno":fno, "n_bits":n_bits} # electons

    

    
    
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



def venus_dict(daynight):
    
    #venus signal
    nadir = hr_nadir_spectra(1)
    if daynight == "d": #gain 2
        px_Wm2um = np.interp(detector["px_m"], nadir["day_um"]/1.0e6, nadir["day_Wm2um"]) #interpolate venus radiance to pixel wavelengths
    if daynight == "n": #gain 1
        px_Wm2um = np.interp(detector["px_m"], nadir["night_um"]/1.0e6, nadir["night_Wm2um"])
    venus = {"b_wm2um":px_Wm2um}
    #manually select OIP radiances
    if daynight == "d":
        venus["b_wm2um"][0:50] = {"2a":23, "2b":18, "4":176}[band]
        venus["b_wm2um"][50:150] = {"2a":21, "2b":17, "4":165}[band]
        venus["b_wm2um"][150:] = {"2a":11, "2b":13, "4":126}[band]
    if daynight == "n":
        venus["b_wm2um"][0:50] = {"1":0.1051, "2a":0.1507, "2b":0.0451, "3":0.2635}[band]
        venus["b_wm2um"][50:150] = {"1":0.0314, "2a":0.0702, "2b":0.0252, "3":0.0667}[band]
        venus["b_wm2um"][150:] = {"1":0.0022, "2a":0.0298, "2b":0.0099, "3":0.0133}[band]

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
        

def warmsection_dict():
    #TODO: add contribution through slit (small)
    warmsection = {}

    return warmsection



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
    
    return signal










st.title("VenSpec-H SNR model")
st.write("Calculate the signal and noise contributions and plot the SNR for each band")



bands_dict = spectral_cal()

for band in bands_dict.keys():
    
    for daynight in bands_dict[band]["daynight"]:
        
        band_dict = bands_dict[band]
        
        header = "SNR model for band %s %s" %(band, {"d":"dayside nadir", "n":"nightside nadir"}[daynight])
        header2 = header.lower().replace(" ","-")
        st.header(header)
        st.sidebar.markdown(f"[{header}](#{header2})", unsafe_allow_html=True)
        
        
        t = t_dict()
        transmittance = transmittance_dict()
        
        
        detector = detector_dict(band_dict, daynight)
        venus = venus_dict(daynight)
        cold_shield = cold_shield_dict(detector, t)
        det_window = det_window_dict(detector, t)
        spectrometer = spectrometer_dict(detector, t)
        warmsection = warmsection_dict()
        signal = signal_dict()
        
        


        
        
        write_variables(t, "Temperature", detector)
        write_variables(detector, "Detector", detector)
        write_variables(venus, "Venus", detector)
        write_variables(cold_shield, "Cold shield", detector)
        write_variables(det_window, "Detector window", detector)
        write_variables(spectrometer, "Spectrometer (cold) section", detector)
        write_variables(warmsection, "Warm section", detector)
        write_variables(signal, "Signal", detector)
        
        indices = np.array([0, 50, 150])
        print(band, daynight, signal["snr"][indices])
        
        
        # fig1, ax1 = plt.subplots()
        
        # ax1.plot(px_um, signal["snr"])
        # # ax1.set_yscale("log")
        # fig1.suptitle("SNR band %s (%s)" %(band, daynight))
        # ax1.set_xlabel("Wavelength (um)")
        # ax1.set_ylabel("SNR")
        
        # st.pyplot(fig=fig1)
        
        
        
    
