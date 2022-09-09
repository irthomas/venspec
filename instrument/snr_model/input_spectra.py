# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 16:28:02 2022

@author: iant
"""

import numpy as np

def load_asimut_spectra(self):
    
    
    if self.day_input_filepath != "":
    
        day_nm, day_Wm2um = np.loadtxt(self.day_input_filepath, skiprows=4, unpack=True)
        self.day_um = np.flip(day_nm) / 1000.
        self.day_Wm2um = np.flip(day_Wm2um)

    if self.night_input_filepath != "":
        night_nm, night_Wm2um = np.loadtxt(self.night_input_filepath, skiprows=4, unpack=True)
        self.night_um = np.flip(night_nm) / 1000.
        self.night_Wm2um = np.flip(night_Wm2um)

