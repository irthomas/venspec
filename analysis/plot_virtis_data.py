# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 17:23:35 2023

@author: iant

PLOT VIRTIS-M-IR DATA

"""
import re
import os
import numpy as np
import matplotlib.pyplot as plt
import struct
import glob
import platform
import h5py
import scipy.stats

if platform.system() == "Windows":
    psa_data_root_path = r"C:\Users\iant\Documents\DATA\venus\virtis\VEX-V-VIRTIS-2-3-V3.0\DATA"
    h5_data_root_path = r"E:\DATA\venus\virtis\hdf5\hdf5_level_1p0a"
elif platform.system() == "Linux":
    psa_data_root_path = r"/bira-iasb/projects/work/EnVision/Data/VIRTIS/psa/VEX-V-VIRTIS-2-3-V3.0/DATA"
    h5_data_root_path = r"/bira-iasb/projects/work/EnVision/Data/VIRTIS/hdf5/hdf5_level_1p0a/"



data_filepaths = glob.glob(f'{h5_data_root_path}{os.sep}**{os.sep}*.h5', recursive=True)



degree_binning = 1


def read_h5_data(data_filepaths):
    band1_d = {"lons":[], "lats":[], "rads":[]}
    for data_filepath in data_filepaths:
        
        with h5py.File(data_filepath) as h5_f:
            
            h5_basename = os.path.basename(data_filepath)
            print(h5_basename)
            
            # print(h5_f["Geometry/Point0"].keys())
        
            sza = h5_f["Geometry/Point0/SZA"][...]
            
            
            night_ixs = np.where(sza > 100.0)[0]
            
            x = h5_f["Science/X"][...]
            y = h5_f["Science/Y"][night_ixs, :]
    
            lats = h5_f["Geometry/Point0/Lat"][night_ixs]
            lons = h5_f["Geometry/Point0/Lon"][night_ixs]
            
            band1_ixs = np.searchsorted(x, 1.17)
            
            # print(min(sza), max(sza))
            
            y_band1 = y[:, band1_ixs]
            
            # plt.figure()
            # plt.plot(sza[night_ixs], y_band1)
            
        band1_d["lons"].extend(lons.tolist())
        band1_d["lats"].extend(lats.tolist())
        band1_d["rads"].extend(y_band1.tolist())
        
        # plt.plot(y_band1)

    for key in band1_d.keys():
        band1_d[key] = np.asarray(band1_d[key])

    return band1_d


band1_d = read_h5_data(data_filepaths)
        
print("Binning data")
bins = [int(360/degree_binning), int(180/degree_binning)]


H, xedges, yedges, binnumber = scipy.stats.binned_statistic_2d(band1_d["lons"], band1_d["lats"], band1_d["rads"], statistic='mean', bins=bins)



vmin = np.min(H[np.isfinite(H)])
vmax = np.max(H[np.isfinite(H)])

fig1, ax1 = plt.subplots(figsize=(14, 7))
plot = ax1.imshow(np.flipud(H.T), alpha=1.0, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()])#, vmin=vmin, vmax=vmax)

ax1.set_xlabel("Longitude")
ax1.set_ylabel("Latitude")
ax1.set_xlim((0, 360))
ax1.set_ylim((-90, 90))
cb1 = fig1.colorbar(plot)
cb1.set_label("Radiance W/m2/sr/um", rotation=270, labelpad=10)
ax1.grid()
fig1.tight_layout()
ax1.set_title("VIRTIS-M-IR Radiance in VenSpec-H band 1")
