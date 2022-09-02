# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 17:10:45 2022

@author: iant

VENSPEC SNR CALCULATIONS
"""

import numpy as np
import scipy.constants as spc
import matplotlib.pyplot as plt



def planck(wav, T): #W/m2/sr/m

    a = 2. * spc.h * spc.c**2.
    b =  spc.h* spc.c /(wav * spc.k * T)
    intensity = a / ( (wav**5.) * (np.exp(b) - 1.) )
    return intensity


def trapezium(x, y):
    return np.sum((x[1:] - x[:-1]) * (y[1:] + y[:-1]) / 2)

n_detector_rows = 288
n_detector_pixels = 384

t_cold = 273. - 43. #kelvin
t_warm = 273.
t_detector = 135

um, Wm2um_earth = np.loadtxt("reference_files/astm_solar_spectrum.tsv", skiprows=1, unpack=True)


day_nm, day_Wm2um = np.loadtxt("reference_files/Simulations_Haus_LIDORT_dayside_scattering_molecule_aer_Rayleigh_sp_nonoise.dat", skiprows=4, unpack=True)
night_nm, night_Wm2um = np.loadtxt("reference_files/Simulations_Haus_LIDORT_nightside_scattering_molecule_aer_Rayleigh_sp_nonoise.dat", skiprows=4, unpack=True)

day_um = np.flip(day_nm) / 1000.
day_Wm2um = np.flip(day_Wm2um)

night_um = np.flip(night_nm) / 1000.
night_Wm2um = np.flip(night_Wm2um)




d_sun_venus = 108.53e9 #m
d_sun_earth = 1.5e11 #m


Wm2_earth = trapezium(um, Wm2um_earth)
Wm2_venus = Wm2_earth * (d_sun_earth / d_sun_venus)**2.

Wm2um_venus = Wm2um_earth * (d_sun_earth / d_sun_venus)**2.



um_range = [1.165, 2.480]
um_range_ix = np.searchsorted(um, um_range)
um_range_ixs = np.arange(um_range_ix[0], um_range_ix[1])


# p_sun = 3.828e26 #W total flux
# Wm2_venus = p_sun / (4. * np.pi * d_sun_venus**2.)

d_venus_vsh = 400.0e3 #m

angle_fov = np.deg2rad([6.7/n_detector_rows, 0.13]) #per detector row

a_ifootprint = np.prod(2.0 * (np.tan(angle_fov/2.0) * d_venus_vsh)) #m2  #per detector row

p_ifootprint = Wm2_venus * a_ifootprint #per detector row
Wum_ifootprint = Wm2um_venus * a_ifootprint #per detector row

albedo_cloud = 0.3

p_vsh = p_ifootprint / (4. * np.pi * d_venus_vsh**2.) * albedo_cloud #per detector row
Wum_vsh = Wum_ifootprint / (4. * np.pi * d_venus_vsh**2.) * albedo_cloud #per detector row

e_photon = spc.h * spc.c / (um * 1.0e-6)

n_photons = Wum_vsh / e_photon #per detector row

# plt.plot(um[um_range_ixs], n_photons[um_range_ixs])


band_dict = {
    "1":{"daynight":["n"], "um_range_aim":[1.16, 1.18], "n_px":194+34},
    "2a":{"daynight":["d", "n"], "um_range_aim":[2.34, 2.422], "n_px":34},
    "2b":{"daynight":["d", "n"], "um_range_aim":[2.422, 2.51], "n_px":194},
    "3":{"daynight":["n"], "um_range_aim":[1.704, 1.747], "n_px":194+34},
    "4":{"daynight":["d"], "um_range_aim":[1.367, 1.394], "n_px":194+34},
         }


# May not be correct in new design
gr_alpha = np.deg2rad(63.43)
gr_beta_centre = np.deg2rad(-1.12)
gr_spacing_um = 79. #um
gr_blaze = np.deg2rad(31.12)

gr_alpha = np.deg2rad(68.57)
gr_beta_centre = np.deg2rad(2.75)
gr_spacing_um = 59.24 #um
gr_blaze = np.deg2rad(35.66)


gr_sin_alpha = np.sin(gr_alpha)
gr_sin_beta = np.sin(gr_beta_centre)

gr_d_sin_alpha_sin_beta = gr_spacing_um * (gr_sin_alpha + gr_sin_beta) #um


fsr = 144.88 #cm-1
resolving_power = 11000.
slit_width = 1. #pixels




plt.figure()
for band in band_dict.keys():
    
    #desired spectral range
    um_range_aim = band_dict[band]["um_range_aim"]
    
    #find corresponding diffraction order
    #order = d(sin alpha + sin beta)/lambda
    order_calc = (10000. / np.mean(um_range_aim)) / fsr
    # order_calc = gr_d_sin_alpha_sin_beta/np.mean(um_range_aim)
    order = int(np.round(order_calc))
    print(np.mean(um_range_aim), "um", order_calc, "->", order)
    
    #wavelength of centre of diffraction order
    # band_centre_calc = gr_d_sin_alpha_sin_beta / order_calc
    band_centre = (10000. / order) / fsr
    print("centre um=", band_centre)
    
    # band_centre = gr_d_sin_alpha_sin_beta / order
    band_dict[band]["order"] = order
    band_dict[band]["um_centre"] = band_centre
    band_dict[band]["cm-1_centre"] = 10000. / band_centre
    
    #diffration order wavenumber range: centre +- half the FSR
    order_range_cm = np.array([band_dict[band]["cm-1_centre"] + fsr/2, band_dict[band]["cm-1_centre"] - fsr/2])
    band_dict[band]["order_range_um"] = 10000. / order_range_cm
    
    
    # # print FSRs
    # if band != "1":
    #     fsr = (10000.0/band_centre_calc - 10000.0/1.17) / (order_calc - 59.07053614734678)
    #     print(fsr)
    
    
    # um_ix = np.searchsorted(um, um_range_aim)
    # um_ixs = np.arange(um_ix[0], um_ix[1]+1)
    # band_dict[band]["um_ix"] = um_ix
    # band_dict[band]["um_ixs"] = um_ixs
    
    #spectral resolution: width of gaussian encompassing wavelengths received by a pixel
    #resolving power = lambda/delta_lambda
    band_dict[band]["delta_lambda_um"] = band_centre / resolving_power
    
    #spectral sampling: wavelength between centre of adjacent pixels
    band_dict[band]["px_sampling_um"] = band_dict[band]["delta_lambda_um"] / slit_width
    # print("sampling=", band_dict[band]["px_sampling_um"]*1000.)
    
    #total wavelength range allowed on detector
    band_dict[band]["detector_um_width"] = band_dict[band]["px_sampling_um"] * n_detector_pixels
    band_dict[band]["detector_um_range"] = np.array([
        band_centre - band_dict[band]["detector_um_width"]/2., 
        band_centre + band_dict[band]["detector_um_width"]/2.
        ])
    
    #if FSR is bigger than detector, light cannot be detected => reduce band spectral range
    #if detector width is bigger than FSR, extra pixels will not be illuminated -> reduce band spectral range
    band_um_range = np.array([
        np.max([band_dict[band]["detector_um_range"][0], band_dict[band]["order_range_um"][0]]),
        np.min([band_dict[band]["detector_um_range"][1], band_dict[band]["order_range_um"][1]]),
        ])
    band_dict[band]["band_um_range"] = band_um_range
    
    print(band_dict[band]["detector_um_range"][0], band_dict[band]["order_range_um"][0], "->", band_dict[band]["band_um_range"][0])
    print(band_dict[band]["detector_um_range"][1], band_dict[band]["order_range_um"][1], "->", band_dict[band]["band_um_range"][1])

    #centre wavelengths on each pixel
    band_dict[band]["detector_ums"] = np.arange(band_um_range[0], band_um_range[1], band_dict[band]["px_sampling_um"])
    
    print(len(band_dict[band]["detector_ums"]))
    

    n_px_band = band_dict[band]["n_px"]
    
    
    if "d" in band_dict[band]["daynight"]:
        day_ixs = np.searchsorted(day_um, band_dict[band]["detector_ums"])
        plt.plot(day_um[day_ixs], day_Wm2um[day_ixs], label="%s day" %band)
    if "n" in band_dict[band]["daynight"]:
        night_ixs = np.searchsorted(night_um, band_dict[band]["detector_ums"])
        plt.plot(night_um[night_ixs], night_Wm2um[night_ixs], label="%s night" %band)
    

plt.plot(day_um, day_Wm2um, alpha=0.3, label="Dayside spectrum")
plt.plot(night_um, night_Wm2um, alpha=0.3, label="Nightside spectrum")

plt.xlabel("Wavelength um")
plt.ylabel("Radiance W/m2/sr/um")
    
    # plt.plot(um[um_ixs], n_photons[um_ixs] * n_px_band)
plt.grid()
plt.legend()
plt.yscale("log")
    
tb_um = um[um_range_ixs]
# tb_um = np.arange(0.1, 3.0, 0.01)
tb_m = tb_um * 1.0e-6

b_warm = planck(tb_m, t_warm)
b_cold = planck(tb_m, t_cold)
b_detector = planck(tb_m, t_detector)

fno = 3.
solid_angle_cold_shield = 2. * np.pi * np.sqrt(1. - 1./(4. * fno**2.))
solid_angle_ = np.pi / (4 * fno**2.)
solid_angle = np.pi / (4 * fno**2.)


# tb = 

e_photon = spc.h * spc.c / tb_m #energy per photon
r_photons_warm = b_warm / e_photon
r_photons_cold = b_cold / e_photon
r_photons_detector = b_detector / e_photon

plt.figure()
plt.plot(tb_um, b_warm) #W/m2/sr/m-1
plt.plot(tb_um, b_cold) #W/m2/sr/m-1
plt.plot(tb_um, b_detector) #W/m2/sr/m-1
plt.yscale("log")

plt.figure()
plt.plot(tb_um, r_photons_warm) #photons/s/m2/sr/m-1
plt.plot(tb_um, r_photons_cold) #photons/s/m2/sr/m-1
plt.plot(tb_um, r_photons_detector) #photons/s/m2/sr/m-1
plt.yscale("log")
