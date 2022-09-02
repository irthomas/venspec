# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 15:42:15 2022

@author: iant
"""

import matplotlib.pyplot as plt

import numpy as np
import streamlit as st
import mpld3
import streamlit.components.v1 as components

from instrument.snr.hr_nadir_spectra import hr_nadir_spectra, extend_to_2p6um
from instrument.snr.spectrometer_calculations import spectral_cal, fsr


st.title("VenSpec-H Signal on each Pixel")

st.write("Increase the spectral resolution of the simulated day and night radiances, then convolve to the pixel's ILS functions and a sinc2 grating blaze function")

def gaussian(x, a, b, c, d):
    return a * np.exp(-((x - b)/c)**2.0) + d

def F_blaze(x, blazef, blazew):
    """wavenumber grid, wavenumber peak, FSR wavenumber"""
    
    dx = x - blazef
    F = np.sinc((dx) / blazew)**2
    return F



def blaze_um(pixel_ums, blaze_centre, fsr):
    
    cm_grid = 10000. / pixel_ums
    
    F = F_blaze(cm_grid, blaze_centre, fsr)
    
    return F



nadir = hr_nadir_spectra(1)
nadir = extend_to_2p6um(nadir)

band_dict = spectral_cal()

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

for band in band_dict.keys():
    
    detector_ums = band_dict[band]["detector_ums"]
    blaze_centre = band_dict[band]["cm-1_centre"]
    
    blaze = blaze_um(detector_ums, blaze_centre, fsr)
    

    for daynight in band_dict[band]["daynight"]:
        px_rad = np.zeros_like(detector_ums)
        
        
        for px_ix, px_um in enumerate(detector_ums):
            
            gauss_width = band_dict[band]["delta_lambda_um"] / 2.355
    
            #calculate relative wavelength range and grid of gaussian on first loop
            if px_ix == 0:
                um_grid_start = -gauss_width * 3
                um_grid_end = gauss_width * 3
                # um_grid = np.arange(um_grid_start, um_grid_end, 0.000001)
    
            
            #find indices in hr grid covering +- 0.5 cm-1 of pixel wavenumber
            um_grid = nadir["day_um"]
                
            
            nu_hr_start_ix = np.searchsorted(um_grid, px_um + um_grid_start) #start index
            nu_hr_end_ix = np.searchsorted(um_grid, px_um + um_grid_end) #end index
    
            px_grid = um_grid[nu_hr_start_ix:nu_hr_end_ix] - px_um
            
            if daynight == "d":
                nadir_rad = nadir["day_Wm2um"][nu_hr_start_ix:nu_hr_end_ix]
                ax = ax1
            if daynight == "n":
                nadir_rad = nadir["night_Wm2um"][nu_hr_start_ix:nu_hr_end_ix]
                ax = ax2

            ils = gaussian(px_grid, 1.0, 0.0, gauss_width, 0.0)
        
            px_rad[px_ix] = np.sum(ils * nadir_rad) * blaze[px_ix] / np.sum(ils)

        ax.plot(detector_ums, px_rad, label="Band %s (%s)" %(band, daynight))

        # plt.yscale("log")
        
ax1.plot(nadir["day_um"], nadir["day_Wm2um"], "k", alpha=0.3)
ax2.plot(nadir["night_um"], nadir["night_Wm2um"], "k", alpha=0.3)

fig1.suptitle("Dayside Nadir")
fig2.suptitle("Nightside Nadir")

ax1.grid()
ax2.grid()
ax1.legend()
ax2.legend()

ax1.set_xlabel("Wavelength um")
ax2.set_xlabel("Wavelength um")

ax1.set_ylabel("Radiance W/m2/sr/um")
ax2.set_ylabel("Radiance W/m2/sr/um")

ax1.set_yscale("log")
ax2.set_yscale("log")

st.subheader("Dayside nadir")
fig1_html = mpld3.fig_to_html(fig1)
components.html(fig1_html, height=600)

st.subheader("Nightside nadir")
fig2_html = mpld3.fig_to_html(fig2)
components.html(fig2_html, height=600)
