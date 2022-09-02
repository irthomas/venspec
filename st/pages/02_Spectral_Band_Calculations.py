# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 15:35:45 2022

@author: iant
"""

import matplotlib.pyplot as plt

import numpy as np
import streamlit as st
import mpld3
import streamlit.components.v1 as components



from instrument.snr.hr_nadir_spectra import hr_nadir_spectra
from instrument.snr.spectrometer_calculations import spectral_cal

st.title("VenSpec-H Spectral Band Calculations: ILS Figures")

st.subheader("Pixel Wavelengths")
st.write("For each spectral band, calculate the diffraction order, min and max wavelengths within the order and hitting the detector")
st.write("Zoom in to see the Gaussian ILS of each pixel")

def gaussian(x, a, b, c, d):
    return a * np.exp(-((x - b)/c)**2.0) + d


nadir = hr_nadir_spectra(1)
band_dict = spectral_cal()

for band in band_dict.keys():

    fig2, ax2 = plt.subplots(constrained_layout=True)
    
    #for plotting
    label_1 = False
    label_2 = False
    for i, um in enumerate(band_dict[band]["detector_px_centres"]):
        
        gauss_width = band_dict[band]["delta_lambda_um"] / 2.355

        #calculate wavelength range of gaussian on first loop
        if i == 0:
            um_grid_start = -gauss_width * 3
            um_grid_end = gauss_width * 3
            
        #if within real detector
        if um in band_dict[band]["detector_ums"]:
            alpha = 0.3
            if not label_1:
                label = "Pixel within FSR"
                label_1 = True
            else:
                label = ""
        else:
            alpha = 0.1
            if not label_2:
                label = "Pixel outside FSR"
                label_2 = True
            else:
                label = ""
        
        um_grid = np.arange(um_grid_start, um_grid_end, 0.000001)
        gauss = gaussian(um_grid, 1.0, 0.0, gauss_width, 0.0)
        
        
        ax2.plot(um_grid + um, gauss, c="k", alpha=alpha, label=label)
        # ax2.plot([um, um], [0, 1], c="k", alpha=alpha)
        

    if "d" in band_dict[band]["daynight"]:
        day_ixs = np.searchsorted(nadir["day_um"], band_dict[band]["detector_ums"])
        ax2.plot(nadir["day_um"][day_ixs], nadir["day_Wm2um"][day_ixs]/np.max(nadir["day_Wm2um"][day_ixs]), label="Band %s day normalised" %band)
    if "n" in band_dict[band]["daynight"]:
        night_ixs = np.searchsorted(nadir["night_um"], band_dict[band]["detector_ums"])
        ax2.plot(nadir["night_um"][night_ixs], nadir["night_Wm2um"][night_ixs]/np.max(nadir["night_Wm2um"][night_ixs]), label="Band %s night normalised" %band)
    ax2.legend(loc="upper right")
    
    ax2.set_xlabel("Wavelength um")
    ax2.set_ylabel("Normalised radiance / Pixel responses")

    ax2.set_title("Band %s, diffraction order %i: %0.3f-%0.3fum (%i pixels used)" %(band, band_dict[band]["order"], *band_dict[band]["band_um_range"], len(band_dict[band]["detector_ums"])))
    
    fig2_html = mpld3.fig_to_html(fig2)
    components.html(fig2_html, height=600)
    # st.pyplot(fig=fig2)

