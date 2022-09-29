# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 16:08:12 2022

@author: iant
"""

import os
import matplotlib.pyplot as plt


# def save_generic_output(self, band, daynight):

#     #append band and daynight
#     output_filename = os.path.basename(self.output_filepath)
#     output_split = os.path.splitext(output_filename)
#     output_dir = os.path.dirname(self.output_filepath)

#     title = output_split[0] + "_band_%s_%s_all_columns_t_cold_section=%0.0fK_gain=%i_it=%0.3fs_nbits=%i_nrows=%i_real_RP=%i_ils_convolution=%s" \
#             %(band, {"d":"day", "n":"night"}[daynight], self.t_cold_section, self.gain, self.integration_time, \
#               self.n_bits, self.n_rows, self.real_resolving_power, self.ils_convolution)

#     output_filepath = os.path.join(output_dir, title + ".txt")

#     with open(output_filepath, "w") as f:
#         f.write("# VenSpec-H SNR Model\n")
#         f.write("# Dayside input:%s\n" %os.path.basename(self.day_input_filepath))
#         f.write("# Nightside input:%s\n" %os.path.basename(self.night_input_filepath))
        
#         f.write("# gain=%i, integration time=%0.3f, n_bits=%i, n_rows=%i, t_cold_section=%0.0f, real_RP=%i, ils_convolution=%s\n" \
#                 %(self.gain, self.integration_time, self.n_bits, self.n_rows, self.t_cold_section, self.real_resolving_power, self.ils_convolution))
#         f.write("# Wavelength (um), Radiance (W.m-2.sr-1.um-2), SNR (1 pixel), SNR (binned), Radiance error 1 pixel (W.m-2.sr-1.um-2), Radiance error binned (W.m-2.sr-1.um-2)\n")
        
#         for i in range(len(self.px_um)):
#             venus_rad = self.venus_wm2um[i]
#             venus_rad_error_px = venus_rad / self.signal["snr"][i]
#             venus_rad_error_binned = venus_rad / self.signal["snr_binned"][i]
            
            
#             f.write("%0.5f, %0.3e, %0.3f, %0.5f, %0.5e, %0.5e\n" \
#                     %(self.px_um[i], venus_rad, self.signal["snr"][i], self.signal["snr_binned"][i], venus_rad_error_px, venus_rad_error_binned))








def prepare_plot(self):

    if len(self.bands) == 7:
        self.fig1, self.axes1 = plt.subplots(nrows=4, ncols=2, figsize=(15, 10), tight_layout=True)
    elif len(self.bands) > 4:
        self.fig1, self.axes1 = plt.subplots(nrows=3, ncols=2, figsize=(15, 10), tight_layout=True)
    else:
        self.fig1, self.axes1 = plt.subplots(nrows=len(self.bands), figsize=(15, 10), tight_layout=True)

    if len(self.bands) == 1:
        self.axes1 = [self.axes1]
    else:
        self.axes1 = self.axes1.flatten()







def save_plot(self, band, daynight, band_ix):

    output_filename = os.path.basename(self.output_filepath)
    output_split = os.path.splitext(output_filename)
    output_dir = os.path.dirname(self.output_filepath)

    self.fig1.suptitle(output_split[0])

            
    self.axes1[band_ix].plot(self.px_um, self.signal["snr"], label="SNR 1 pixel")
    self.axes1[band_ix].plot(self.px_um, self.signal["snr_binned"], label="SNR binned")
    self.axes1[band_ix].legend(loc="lower right")
    self.axes1[band_ix].set_xlabel("Wavelength (um)")
    self.axes1[band_ix].set_ylabel("SNR")
    self.axes1[band_ix].grid()
    self.axes1[band_ix].set_title("Band %s %s" %(band, {"d":"day", "n":"night"}[daynight]))
    
    self.axes1[band_ix].set_yscale("log")
    
    title = output_split[0] + "_radiance_scalar=%0.2f_t_cold_section=%0.0fK_real_RP=%i_ils_convolution=%s" \
            %(self.input_radiance_scalar, self.t_cold_section, self.real_resolving_power, self.ils_convolution)
    output_filepath = os.path.join(output_dir, title + ".png")
    self.fig1.savefig(output_filepath)



    plt.subplots(figsize=(13, 5), tight_layout=True)
    
    title = output_split[0] + "_band_%s_%s_radiance_scalar=%0.2f_t_cold_section=%0.0fK_gain=%i_it=%0.3fs_nbits=%i_nrows=%i_real_RP=%i_ils_convolution=%s" \
            %(band, {"d":"day", "n":"night"}[daynight], self.input_radiance_scalar, self.t_cold_section, self.gain, self.integration_time, \
              self.n_bits, self.n_rows, self.real_resolving_power, self.ils_convolution)
    output_filepath = os.path.join(output_dir, title + ".png")

    
    # plt.plot(self.px_um, self.venus_wm2um)
    plt.errorbar(self.px_um, self.venus_wm2um, ecolor="k", capsize=2, yerr=(self.venus_wm2um / self.signal["snr_binned"]))
    plt.xlabel("Wavelength (um)")
    plt.ylabel("Radiance (W.m-2.sr-1.um-2)")
    plt.grid()

    plt.title("Band %s %s pixel radiance and error (detector rows binned)\ngain=%i, integration time=%0.3fs, n_bits=%i, n_rows=%i, scalar=%0.2f, t_cold_section=%0.0fK, real_RP=%i, ils_convolution=%s" \
            %(band, {"d":"day", "n":"night"}[daynight], self.gain, self.integration_time, self.n_bits, self.n_rows, self.input_radiance_scalar, self.t_cold_section, self.real_resolving_power, self.ils_convolution))
    plt.savefig(output_filepath)
    plt.close()