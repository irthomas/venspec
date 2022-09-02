# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 15:35:45 2022

@author: iant
"""


import numpy as np
import streamlit as st
import pandas as pd



from instrument.snr.hr_nadir_spectra import hr_nadir_spectra
from instrument.snr.spectrometer_calculations import spectral_cal, pixel_illumination_table

st.title("VenSpec-H Spectral Band Calculations: Table")

st.subheader("Pixel Wavelengths (nm) of each band")
st.write("Zero = order not illuminated on this pixel")


nadir = hr_nadir_spectra(1)
band_dict = spectral_cal()


table = pixel_illumination_table()
df = pd.DataFrame(
    table, columns=["Band %s (order %i)" %(band, band_dict[band]["order"]) for band in band_dict.keys()])
st.table(df)