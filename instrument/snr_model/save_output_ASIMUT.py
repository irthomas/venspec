# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 16:08:12 2022

@author: iant
Adapted to be used directly with ASIMUT by srobert
Change of plans, we save clean radiances on one file and noise in another file in ASIMUT units
"""

import os
import numpy as np
# import matplotlib.pyplot as plt





def convert_units_add_noise(venus_um, venus_wm2um, venus_snr):
    #convert um to cm-1
    venus_cm1 = 10000.0 / venus_um
    
    #convert W/m2/um to W/cm2/cm-1
    #[W/cm^2/sr/um] = [W/cm^2/sr/cm^-1] * 10^4 / λ^2 where λ in um
    #[W/m^2/sr/um] = [W/cm^2/sr/um] * 100^2 = [W/cm^2/sr/cm^-1] *10^8 / λ^2 where λ in um
    venus_wcm2um = venus_wm2um / 100.0**2.0
    venus_wcm2cm1 = venus_wcm2um * 1.0e4 / venus_cm1**2.0
    
    #error (snr is same in cm-1 or microns)
    venus_wcm2cm1_error = venus_wcm2cm1 / venus_snr
    #generate random error
    random_noise = np.random.normal(0.0, venus_wcm2cm1_error, len(venus_wcm2cm1_error))
    
    #add random error to radiance
    venus_wcm2cm1_with_noise = venus_wcm2cm1 + random_noise

    #test plotting
    # plt.figure()
    # plt.plot(venus_cm1, venus_wcm2cm1)
    # plt.plot(venus_cm1, venus_wcm2cm1 - venus_wcm2cm1_error)
    # plt.plot(venus_cm1, venus_wcm2cm1 + venus_wcm2cm1_error)
    # plt.plot(venus_cm1, venus_wcm2cm1_with_noise)
    # plt.yscale("log")
    # plt.figure()
    # plt.plot(venus_cm1, venus_wcm2cm1_error)
    # plt.plot(venus_cm1, -venus_wcm2cm1_error)
    # plt.plot(venus_cm1, random_noise)
    
    #output: wavenumbers, venus radiance, radiance error, venus radiance with noise added
    return np.flip(venus_cm1), np.flip(venus_wcm2cm1), np.flip(venus_wcm2cm1_error), np.flip(venus_wcm2cm1_with_noise)






# def save_cm1_radiance_with_noise(self, band, daynight):
#     """Convert to wavenumbers and save file with 2 columns: wavenumbers and radiances with added random noise"""

#     #append band and daynight and model input parameters to filename
#     output_filename = os.path.basename(self.output_filepath)
#     output_split = os.path.splitext(output_filename)
#     output_dir = os.path.dirname(self.output_filepath)

#     #find input spectrum filename    
#     if daynight == "d":
#         input_filename = os.path.basename(self.day_input_filepath)
#     elif daynight == "n":
#         input_filename = os.path.basename(self.night_input_filepath)

#     """Radiance error and noise for 1 pixel"""
#     n_rows = 1
#     snr = self.signal["snr"]

#     #convert output to correct units and add noise
#     venus_cm1, venus_wcm2cm1, venus_wcm2cm1_error, venus_wcm2cm1_with_noise = convert_units_add_noise(self.px_um, self.venus_wm2um, snr)

#     title = output_split[0] + "_band_%s_%s_NoisyRadiance1px_t_cold_section=%0.0fK_gain=%i_it=%0.3fs_nbits=%i_nrows=%i_real_RP=%i_ils_convolution=%s" \
#             %(band, {"d":"day", "n":"night"}[daynight], self.t_cold_section, self.gain, self.integration_time, \
#               self.n_bits, n_rows, self.real_resolving_power, self.ils_convolution)

#     output_filepath = os.path.join(output_dir, title + ".txt")

#     #write file
#     with open(output_filepath, "w") as f:
        
#         #write file header
#         f.write("%% Simulated radiances including noise based on the %s and using the SNR model for T= %0.1fK\n" %(input_filename, self.t_cold_section))
#         f.write("%% gain=%i, integration time=%0.3f, n_bits=%i, n_rows=%i, t_cold_section=%0.0f, real_RP=%i, ils_convolution=%s\n" \
#                 %(self.gain, self.integration_time, self.n_bits, n_rows, self.t_cold_section, self.real_resolving_power, self.ils_convolution))
#         f.write("%% Wavenumber (cm-1)\tNoisy radiance (W.cm-2.cm.sr-1)\n")
        
#         #write wavenumbers and radiance spectrum with added noise, line by line
#         for i in range(len(venus_cm1)):
#             f.write("%0.5f\t%0.5e\n" %(venus_cm1[i], venus_wcm2cm1_with_noise[i]))

#     """Radiance error and noise for binned pixels"""
#     n_rows = self.n_rows
#     snr = self.signal["snr_binned"]

#     #convert output to correct units and add noise
#     venus_cm1, venus_wcm2cm1, venus_wcm2cm1_error, venus_wcm2cm1_with_noise = convert_units_add_noise(self.px_um, self.venus_wm2um, snr)

#     title = output_split[0] + "_band_%s_%s_NoisyRadianceBinned_t_cold_section=%0.0fK_gain=%i_it=%0.3fs_nbits=%i_nrows=%i_real_RP=%i_ils_convolution=%s" \
#             %(band, {"d":"day", "n":"night"}[daynight], self.t_cold_section, self.gain, self.integration_time, \
#               self.n_bits, n_rows, self.real_resolving_power, self.ils_convolution)

#     output_filepath = os.path.join(output_dir, title + ".txt")

#     #write file
#     with open(output_filepath, "w") as f:
        
#         #write file header
#         f.write("%% Simulated radiances including noise based on the %s and using the SNR model for T= %0.1fK\n" %(input_filename, self.t_cold_section))
#         f.write("%% gain=%i, integration time=%0.3f, n_bits=%i, n_rows=%i, t_cold_section=%0.0f, real_RP=%i, ils_convolution=%s\n" \
#                 %(self.gain, self.integration_time, self.n_bits, n_rows, self.t_cold_section, self.real_resolving_power, self.ils_convolution))
#         f.write("%% Wavenumber (cm-1)\tNoisy radiance (W.cm-2.cm.sr-1)\n")
        
#         #write wavenumbers and radiance spectrum with added noise, line by line
#         for i in range(len(venus_cm1)):
#             f.write("%0.5f\t%0.5e\n" %(venus_cm1[i], venus_wcm2cm1_with_noise[i]))





def save_cm1_clean_radiance(self, band, daynight):
    """Convert to wavenumbers and save file with 2 columns: wavenumbers and clean radiances
    The radiances are the same for 1 detector row or when binned, so only 1 file is written"""

    #append band and daynight and model input parameters to filename
    output_filename = os.path.basename(self.output_filepath)
    output_split = os.path.splitext(output_filename)
    output_dir = os.path.dirname(self.output_filepath)

    #find input spectrum filename    
    if daynight == "d":
        input_filename = os.path.basename(self.day_input_filepath)
    elif daynight == "n":
        input_filename = os.path.basename(self.night_input_filepath)

    """Radiance error and noise for 1 pixel"""
    snr = self.signal["snr"]


    #convert output to correct units and add noise
    venus_cm1, venus_wcm2cm1, venus_wcm2cm1_error, venus_wcm2cm1_with_noise = convert_units_add_noise(self.px_um, self.venus_wm2um, snr)

    title = output_split[0] + "_band_%s_%s_Radiance_scalar=%0.2f_t_cold_section=%0.0fK_gain=%i_it=%0.3fs_nbits=%i_real_RP=%i_ils_convolution=%s" \
            %(band, {"d":"day", "n":"night"}[daynight], self.input_radiance_scalar, self.t_cold_section, self.gain, self.integration_time, \
              self.n_bits, self.real_resolving_power, self.ils_convolution)

    output_filepath = os.path.join(output_dir, title + ".txt")

    #write file
    with open(output_filepath, "w") as f:
        
        #write file header
        f.write("%% Simulated clean radiances based on the input file %s and using the SNR model for T= %0.1fK\n" %(input_filename, self.t_cold_section))
        f.write("%% gain=%i, integration time=%0.3f, n_bits=%i, radiance_scalar=%0.2f, t_cold_section=%0.0f, real_RP=%i, ils_convolution=%s\n" \
                %(self.gain, self.integration_time, self.n_bits, self.input_radiance_scalar, self.t_cold_section, self.real_resolving_power, self.ils_convolution))
        f.write("%% Wavenumber (cm-1)\tRadiance (W.cm-2.cm.sr-1)\n")
        
        #write wavenumbers and radiance spectrum with added noise, line by line
        for i in range(len(venus_cm1)):
            f.write("%0.5f\t%0.5e\n" %(venus_cm1[i], venus_wcm2cm1[i]))







def save_cm1_radiance_error(self, band, daynight):
    """Convert to wavenumbers and save file with 2 columns: wavenumbers and random noise on the radiances"""

    #append band and daynight and model input parameters to filename
    output_filename = os.path.basename(self.output_filepath)
    output_split = os.path.splitext(output_filename)
    output_dir = os.path.dirname(self.output_filepath)

    if daynight == "d":
        input_filename = os.path.basename(self.day_input_filepath)
    elif daynight == "n":
        input_filename = os.path.basename(self.night_input_filepath)

    """Radiance error for 1 pixel"""
    n_rows = 1
    snr = self.signal["snr"]

    #convert output to correct units and add noise
    venus_cm1, venus_wcm2cm1, venus_wcm2cm1_error, venus_wcm2cm1_with_noise = convert_units_add_noise(self.px_um, self.venus_wm2um, snr)

    title = output_split[0] + "_band_%s_%s_RadianceError1px_scalar=%0.2f_t_cold_section=%0.0fK_gain=%i_it=%0.3fs_nbits=%i_nrows=%i_real_RP=%i_ils_convolution=%s" \
            %(band, {"d":"day", "n":"night"}[daynight], self.input_radiance_scalar, self.t_cold_section, self.gain, self.integration_time, \
              self.n_bits, n_rows, self.real_resolving_power, self.ils_convolution)

    output_filepath = os.path.join(output_dir, title + ".txt")

    #write file
    with open(output_filepath, "w") as f:
        
        #write file header
        f.write("%% Noise values based on the %s and using the SNR model for T= %0.1fK\n" %(input_filename, self.t_cold_section))
        f.write("%% gain=%i, integration time=%0.3f, n_bits=%i, n_rows=%i, scalar=%0.2f, t_cold_section=%0.0f, real_RP=%i, ils_convolution=%s\n" \
                %(self.gain, self.integration_time, self.n_bits, n_rows, self.input_radiance_scalar, self.t_cold_section, self.real_resolving_power, self.ils_convolution))
        f.write("%% Wavenumber (cm-1)\tRadiance Error (W.cm-2.cm.sr-1)\n")
        
        
        #write wavenumbers and radiance spectrum with added noise, line by line
        for i in range(len(venus_cm1)):
            f.write("%0.5f\t%0.5e\n" %(venus_cm1[i], venus_wcm2cm1_error[i]))

    """Radiance error for binned pixels"""
    n_rows = self.n_rows
    snr = self.signal["snr_binned"]

    #convert output to correct units and add noise
    venus_cm1, venus_wcm2cm1, venus_wcm2cm1_error, venus_wcm2cm1_with_noise = convert_units_add_noise(self.px_um, self.venus_wm2um, snr)

    title = output_split[0] + "_band_%s_%s_RadianceErrorBinned_scalar=%0.2f_t_cold_section=%0.0fK_gain=%i_it=%0.3fs_nbits=%i_nrows=%i_real_RP=%i_ils_convolution=%s" \
            %(band, {"d":"day", "n":"night"}[daynight], self.input_radiance_scalar, self.t_cold_section, self.gain, self.integration_time, \
              self.n_bits, n_rows, self.real_resolving_power, self.ils_convolution)

    output_filepath = os.path.join(output_dir, title + ".txt")

    #write file
    with open(output_filepath, "w") as f:
        
        #write file header
        f.write("%% Noise values based on the %s and using the SNR model for T= %0.1fK\n" %(input_filename, self.t_cold_section))
        f.write("%% gain=%i, integration time=%0.3f, n_bits=%i, n_rows=%i, scalar=%0.2f, t_cold_section=%0.0f, real_RP=%i, ils_convolution=%s\n" \
                %(self.gain, self.integration_time, self.n_bits, n_rows, self.input_radiance_scalar, self.t_cold_section, self.real_resolving_power, self.ils_convolution))
        f.write("%% Wavenumber (cm-1)\tRadiance Error (W.cm-2.cm.sr-1)\n")
        
        #write wavenumbers and radiance spectrum with added noise, line by line
        for i in range(len(venus_cm1)):
            f.write("%0.5f\t%0.5e\n" %(venus_cm1[i], venus_wcm2cm1_error[i]))



