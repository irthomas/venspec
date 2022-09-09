# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 11:29:14 2022

@author: iant
"""

import numpy as np

from snr_functions import F_solid_angle, F_planck



def cold_shield(self):
    cold_shield_omega = 2. * np.pi - F_solid_angle(self.fno)
    self.cold_shield = {"trans":1.0, "emis":1.0, "fno":self.fno, "omega":cold_shield_omega}
    
    self.cold_shield["b_m"] = F_planck(self.nu_full_range_m, self.t_detector)
    
    


        
def detector_window(self):
    #TODO: replace with table values. At 228K thermal background changes by just 2 e-/s"""
    reflectance = 0.01
    thickness = 0.1 #cm
    abs_coeff = 0.002 #cm-1 #see page 182 for full table https://link.springer.com/content/pdf/10.1007/978-0-387-85695-7.pdf
    
    #note: typo in OIP formula
    #see equation 1 here: https://avs.scitation.org/doi/am-pdf/10.1116/1.4954211#:~:text=For%20light%20impinging%20normally%20on,1%20%E2%80%93%20R2exp(%2D2%CE%B1d)%5D.
    #for thermal purposes, transmission from window to detector = 1
    det_window_trans = ((1. - reflectance)**2 * np.exp(-abs_coeff * thickness))/(1. - np.exp(-abs_coeff * thickness) * reflectance**2)
    det_window_omega = F_solid_angle(self.fno)
    self.det_window = {"trans":1., "emis":1. - det_window_trans, "omega":det_window_omega}
    self.det_window["b_m"] = F_planck(self.nu_full_range_m, self.t_detector_window)
    


        
def cold_section(self):
    cold_section_omega = F_solid_angle(self.fno)
    # spectrometer = {"trans":1., "emis":1., "omega":cold_section_omega}
    self.cold_section = {"trans":self.transmittance_cold_section, "emis":1., "omega":cold_section_omega}
    
    self.cold_section["b_m"] = F_planck(self.nu_full_range_m, self.t_cold_section)
        

