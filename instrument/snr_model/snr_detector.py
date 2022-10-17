# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 10:56:16 2022

@author: iant
"""


import numpy as np



def detector_settings(self, band, daynight):
    
    max_it = self.max_integration_time
    
    if self.t_cold_section == 273:
        self.gain = 2
        
        if daynight == "d":
            self.integration_time = {"2a":0.948, "2b":0.948, "4":0.365}[band]
        if daynight == "n":
            self.integration_time = {"1":max_it, "2a":max_it, "2b":max_it, "3":max_it}[band]
    
    elif self.t_cold_section == 253:
        
        if daynight == "d":
            self.gain = 2
            self.integration_time = {"2a":0.976, "2b":0.976, "4":0.369}[band]
        if daynight == "n":
            self.gain = 1
            self.integration_time = {"1":max_it, "2a":max_it, "2b":max_it, "3":max_it}[band]
    
    elif self.t_cold_section == 228:
        
        if daynight == "d":
            self.gain = 2
            self.integration_time = {"2a":0.982, "2b":0.982, "4":0.370}[band]
        if daynight == "n":
            self.gain = 1
            self.integration_time = {"1":max_it, "2a":max_it, "2b":max_it, "3":max_it}[band]
    
    elif self.t_cold_section == 213:
        
        if daynight == "d":
            self.gain = 2
            self.integration_time = {"2a":0.982, "2b":0.982, "4":0.370}[band]
        if daynight == "n":
            self.gain = 1
            self.integration_time = {"1":max_it, "2a":max_it, "2b":max_it, "3":max_it}[band]

    else:
        #if cold section temperature is not a given OIP setpoint, just use the cold settings
        
        if daynight == "d":
            self.gain = 2
            self.integration_time = {"2a":0.982, "2b":0.982, "4":0.370}[band]
        if daynight == "n":
            self.gain = 1
            self.integration_time = {"1":max_it, "2a":max_it, "2b":max_it, "3":max_it}[band]
    

    self.readout_noise = {1:125., 2:235.}[self.gain]
    self.full_well = {1:340.0e3, 2:1.0e6}[self.gain]
    
    
    self.n_rows = {
        "1":222,
    	"2a":36,
        "2b":186,
        "3":222,
        "4":222,
        }[band]




def detector_qe():
    """get QE versus wavelength (microns)"""
    qe_plot = np.array([[0.5, 0.], [0.75, 0.], [0.78, 0.1], [0.87, 0.75], [2.27, 0.75], [2.5, 0.69], [2.64, 0.02], [2.66, 0.], [3., 0.]])
    qe_um = qe_plot[:, 0]
    qe = qe_plot[:, 1]
    # plt.plot(qe_plot[:, 0], qe_plot[:, 1])
    return qe_um, qe




def detector_qe_interp(self):
    
    qe_um, qe = detector_qe()
    
    self.qe_px = np.interp(self.px_um, qe_um, qe)
    self.qe_full_range = np.interp(self.nu_full_range_um, qe_um, qe)

    
            
    
    
         
def detector(self):


    self.px_m = self.px_um / 1.0e6
    self.px_delta_lambda_m = self.px_delta_lambda_um / 1.0e6
    

    self.detector_size_um = 24 #um
    self.detector_size_m = self.detector_size_um / 1.0e6 #um

    self.detector_area_um2= self.detector_size_um ** 2.
    self.detector_area_m2 = self.detector_size_m ** 2.


    self.detector_dc_nAcm2 = 0.023 # nA/cm2
    self.detector_dc_nAum2 = self.detector_dc_nAcm2 * (100.)**2. / (1.0e6)**2
    self.detector_dc_Aum2 = self.detector_dc_nAum2 * 1.0e-9
    self.detector_dc_A = self.detector_dc_Aum2 * self.detector_area_um2
    self.detector_dc_es = 6.242e18 * self.detector_dc_A
    self.detector_dc_e = self.detector_dc_es * self.integration_time


    #define spectral range of detector for thermal background
    self.nu_full_range_um = np.arange(0.8, 2.6, 0.0001)
    self.nu_full_range_m = self.nu_full_range_um * 1.0e-6

    

