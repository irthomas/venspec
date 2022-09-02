# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 16:22:26 2022

@author: iant

VENSPEC INSTRUMENT MODEL STREAMLIT
"""

import sys
sys.path.append(r"C:\Users\iant\Dropbox\VenSpec\Python")


import streamlit as st

import matplotlib.pyplot as plt


import numpy as np

from instrument.snr.band_spectra import filter_dict, slit_dict, combined_dict


st.title("VenSpec-H Instrument Model")

st.subheader("Filter passbands")


fig1, (ax1a, ax1b, ax1c) = plt.subplots(figsize=(12, 7), nrows=3, constrained_layout=True)

for filter_name in filter_dict.keys():
    um = np.array(filter_dict[filter_name]["um"])
    um2 = np.zeros_like(um)
    ax1a.plot(um[:, 0], um[:, 1], label="Filter wheel band %s" %filter_name, alpha=1.0)
    ax1a.fill_between(um[:, 0], um[:, 1], um2[:,1], label="Filter wheel band %s" %filter_name, alpha=0.5)
    
for filter_name in slit_dict.keys():
    um = np.array(slit_dict[filter_name]["um"])
    um2 = np.zeros_like(um)
    ax1b.plot(um[:, 0], um[:, 1], linestyle="--", label="Slit band %s" %filter_name, alpha=1.0)
    ax1b.fill_between(um[:, 0], um[:, 1], um2[:,1], label="Slit band %s" %filter_name, alpha=0.2)

for filter_name in combined_dict.keys():
    um = np.array(combined_dict[filter_name]["um"])
    um2 = np.zeros_like(um)
    ax1c.plot(um[:, 0], um[:, 1], linestyle="--", label="Filter and slit band %s" %filter_name, alpha=1.0)
    ax1c.fill_between(um[:, 0], um[:, 1], um2[:,1], label="Filter and slit band %s" %filter_name, alpha=0.5)

    
ax1a.legend(loc="upper right")
ax1b.legend(loc="upper right")
ax1c.legend(loc="upper right")
ax1c.set_xlabel("Wavelength um")
ax1a.set_ylabel("Filter transmission on/off")
ax1b.set_ylabel("Filter transmission on/off")
ax1c.set_ylabel("Filter transmission on/off")

fig1.suptitle("VenSpec-H filter wheel and slit spectra")
st.pyplot(fig=fig1)

    
    
    
