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
        
        # if daynight == "d":
        #     self.gain = 2
        #     self.integration_time = {"2a":0.982, "2b":0.982, "4":0.370}[band]
        # if daynight == "n":
        #     self.gain = 1
        #     self.integration_time = {"1":max_it, "2a":max_it, "2b":max_it, "3":max_it}[band]

        #June 2023 new values
        if daynight == "d":
            self.gain = 2
            self.integration_time = {"2a":1.17, "2b":1.17, "4":0.44}[band]
        if daynight == "n":
            self.gain = 1
            self.integration_time = {"1":11.74, "2a":9.42, "2b":9.42, "3":9.69}[band]


    
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
    

    # self.readout_noise = {1:125., 2:235.}[self.gain]
    # self.full_well = {1:340.0e3, 2:1.0e6}[self.gain]

    self.readout_noise = {1:150., 2:500.}[self.gain]
    self.full_well = {1:0.296e6, 2:0.893e6}[self.gain]
    
    
    # self.n_rows = {
    #     "1":222,
    # 	"2a":36,
    #     "2b":186,
    #     "3":222,
    #     "4":222,
    #     }[band]

    self.n_rows = {
        "1":200,
    	"2a":19,
        "2b":181,
        "3":200,
        "4":200,
        }[band]







def detector_qe_interp(self):
    
    qe_um, qe = detector_qe()
    
    self.qe_px = np.interp(self.px_um, qe_um, qe)
    self.qe_full_range = np.interp(self.nu_full_range_um, qe_um, qe)

    
            
    
def dc(t):
    return (4.29e8) * np.exp(-2788.0 / t) + (1.35e14) * np.exp(-5247.0 / t)

    
    
         
def detector(self):


    self.px_m = self.px_um / 1.0e6
    self.px_delta_lambda_m = self.px_delta_lambda_um / 1.0e6
    

    self.detector_size_um = 24 #um
    self.detector_size_m = self.detector_size_um / 1.0e6 #um

    self.detector_area_um2= self.detector_size_um ** 2.
    self.detector_area_m2 = self.detector_size_m ** 2.


    # self.detector_dc_nAcm2 = 0.023 # nA/cm2
    # self.detector_dc_nAum2 = self.detector_dc_nAcm2 * (100.)**2. / (1.0e6)**2
    # self.detector_dc_Aum2 = self.detector_dc_nAum2 * 1.0e-9
    # self.detector_dc_A = self.detector_dc_Aum2 * self.detector_area_um2
    # self.detector_dc_es = 6.242e18 * self.detector_dc_A
    
    self.detector_dc_nAcm2 = dc(self.t_detector) # nA/cm2
    self.detector_dc_es = self.detector_dc_nAcm2 * 36000.0
    self.detector_dc_e = self.detector_dc_es * self.integration_time


    #define spectral range of detector for thermal background
    self.nu_full_range_um = np.arange(0.8, 2.6, 0.0001)
    self.nu_full_range_m = self.nu_full_range_um * 1.0e-6


# def detector_qe():
#     """get QE versus wavelength (microns)"""
#     qe_plot = np.array([[0.5, 0.], [0.75, 0.], [0.78, 0.1], [0.87, 0.75], [2.27, 0.75], [2.5, 0.69], [2.64, 0.02], [2.66, 0.], [3., 0.]])
#     qe_um = qe_plot[:, 0]
#     qe = qe_plot[:, 1]
#     # plt.plot(qe_plot[:, 0], qe_plot[:, 1])
#     return qe_um, qe

def detector_qe():
    """get QE versus wavelength (microns) inc detector window trans"""
    qe_plot = np.array([
        [0.40, 0.0000000],
        [0.45, 0.0000000],
        [0.50, 0.0000000],
        [0.55, 0.0000000],
        [0.60, 0.0000000],
        [0.65, 0.0000000],
        [0.70, 0.0000000],
        [0.75, 0.0000000],
        [0.80, 0.1803750],
        [0.85, 0.4870725],
        [0.90, 0.7056000],
        [0.95, 0.7276800],
        [1.00, 0.7517500],
        [1.05, 0.7449600],
        [1.10, 0.7391400],
        [1.15, 0.7129500],
        [1.20, 0.7197400],
        [1.25, 0.7003400],
        [1.30, 0.7090700],
        [1.35, 0.7042200],
        [1.40, 0.7061600],
        [1.45, 0.7110100],
        [1.50, 0.7129500],
        [1.55, 0.7129500],
        [1.60, 0.7061600],
        [1.65, 0.7110100],
        [1.70, 0.7071300],
        [1.75, 0.7119800],
        [1.80, 0.7090700],
        [1.85, 0.7187700],
        [1.90, 0.7158600],
        [1.95, 0.7226500],
        [2.00, 0.7226500],
        [2.05, 0.7226500],
        [2.10, 0.7197400],
        [2.15, 0.7236200],
        [2.20, 0.7207100],
        [2.25, 0.7168300],
        [2.30, 0.7090700],
        [2.35, 0.6940800],
        [2.40, 0.6835200],
        [2.45, 0.6700800],
        [2.50, 0.6381189],
        [2.55, 0.4947627],
        [2.60, 0.1567500],
        [2.65, 0.0390000],
        [2.70, 0.0000000],
    ])
    
    qe_um = qe_plot[:, 0]
    qe = qe_plot[:, 1]
    # plt.plot(qe_plot[:, 0], qe_plot[:, 1])
    return qe_um, qe


