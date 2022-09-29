# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 16:28:02 2022

@author: iant
"""

import numpy as np
import matplotlib.pyplot as plt
import os



def load_asimut_spectra(self):
    
    
    if self.day_input_filepath != "":
        
        self.day_um, self.day_Wm2um = load_asimut_file(self.day_input_filepath)
        
        #scale input radiance
        self.day_Wm2um *= self.input_radiance_scalar
    
        # day_nm, day_Wm2um = np.loadtxt(self.day_input_filepath, skiprows=4, unpack=True)
        # self.day_um = np.flip(day_nm) / 1000.
        # self.day_Wm2um = np.flip(day_Wm2um)

    if self.night_input_filepath != "":

        self.night_um, self.night_Wm2um = load_asimut_file(self.night_input_filepath)

        #scale input radiance
        self.night_Wm2um *= self.input_radiance_scalar

        # night_nm, night_Wm2um = np.loadtxt(self.night_input_filepath, skiprows=4, unpack=True)
        # self.night_um = np.flip(night_nm) / 1000.
        # self.night_Wm2um = np.flip(night_Wm2um)




def load_asimut_file(filepath, plot=False):
    

    #load numerical data
    data = np.loadtxt(filepath, comments = "%", unpack=True)
    
    #get header units 
    x_unit, y_unit = np.fromregex(filepath, r"%+\s*@\s*(\S+)\s+\((.*?)\)", dtype={'names': ('x', 'y'), 'formats': ((np.str_,50), (np.str_,50))})[0]
    
    # [A-Za-zα-ωΑ-Ω]
    #convert to correct units
    if x_unit == "cm-1":
        um = 10000. / data[0, :]
    elif x_unit == "nm":
        um = data[0, :] / 1000.
    else:
        print("Error: units of input file %s unknown" %x_unit)
        
    #[W/cm^2/sr/um] = [W/cm^2/sr/cm^-1] * 10^4 / λ^2 where λ in um
    #[W/m^2/sr/um] = [W/cm^2/sr/um] * 100^2 = [W/cm^2/sr/cm^-1] *10^8 / λ^2 where λ in um
    if y_unit == "W.cm-2.cm.sr-1":
    
        #convert to W/m2/sr/um
        Wm2um = data[1, :] * 1.0e8 / um**2
        
    elif y_unit in ["W m-2 sr-1 µm-1", "W m-2 sr-1 Âµm-1"]:
        
        Wm2um = data[1, :]
    else:
        print("Error: units of input file %s unknown" %y_unit)
        
    # plot
    if plot:
        # plt.figure()
        plt.title(os.path.basename(filepath))
        plt.plot(um, Wm2um)
        plt.xlabel("um")
        plt.ylabel("W/m2/sr/um")
        
    return np.flip(um), np.flip(Wm2um)



def load_asimut_spectra_cm1(self):
    
    
    if self.day_input_filepath != "":
        
        self.day_cm1, self.day_Wcm2cm1 = load_asimut_file_cm1(self.day_input_filepath)
    
        #scale input radiance
        self.day_Wcm2cm1 *= self.input_radiance_scalar

        # day_nm, day_Wm2um = np.loadtxt(self.day_input_filepath, skiprows=4, unpack=True)
        # self.day_um = np.flip(day_nm) / 1000.
        # self.day_Wm2um = np.flip(day_Wm2um)

    if self.night_input_filepath != "":

        self.night_cm1, self.night_Wcm2cm1 = load_asimut_file_cm1(self.night_input_filepath)

        #scale input radiance
        self.night_Wcm2cm1 *= self.input_radiance_scalar

        # night_nm, night_Wm2um = np.loadtxt(self.night_input_filepath, skiprows=4, unpack=True)
        # self.night_um = np.flip(night_nm) / 1000.
        # self.night_Wm2um = np.flip(night_Wm2um)



def load_asimut_file_cm1(filepath, plot=False):
    

    #load numerical data
    data = np.loadtxt(filepath, comments = "%", unpack=True)
    
    return data[0, :], data[1, :]
    
    

# for testing
# plt.figure()
# filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\nominal_test2_LidortG_AER_Haus_fact050_contXS_3E8_fact2_HITRAN2020_FBord_125_final_specR_030_nightside_SP1_FEN1_rad_forw.dat"
# um, Wm2um = load_asimut_file(filepath, plot=True)
# filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\Simulations_LIDORT_nightside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"
# um, Wm2um = load_asimut_file(filepath, plot=True)

# filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\nominal_test2_LidortG_AER_Haus_fact050_contXS_3E8_fact2_HITRAN2020_FBord_125_final_specR_040_dayside_SP1_FEN1_rad_forw.dat"
# um, Wm2um = load_asimut_file(filepath, plot=True)
# filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\Simulations_LIDORT_dayside_scattering_molecule_aer_Rayleigh_entirerange_AERNIgn_alt540_Nstreams24_sp_nonoise.dat"
# um, Wm2um = load_asimut_file(filepath, plot=True)

# plt.figure()
# day_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\nominal_test2_LidortG_AER_Haus_fact050_contXS_3E8_fact2_EMlinelist_specR_040_dayside_SP1_FEN1_rad_forw.dat"
# um, Wm2um = load_asimut_file(day_input_filepath, plot=True)
# night_input_filepath = r"C:\Users\iant\Dropbox\VenSpec\Python\reference_files\nominal_test2_LidortG_AER_Haus_fact050_contXS_3E8_fact2_EMlinelist_specR_030_nightside_SP1_FEN1_rad_forw.dat"
# um, Wm2um = load_asimut_file(night_input_filepath, plot=True)
