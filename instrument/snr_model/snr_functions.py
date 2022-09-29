# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 10:44:58 2022

@author: iant
"""


import numpy as np
import scipy.constants as spc




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

            
