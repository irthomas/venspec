# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 13:53:44 2023

@author: iant

VENSPEC-H BINNING AND DETECTOR PLOTS
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import numpy as np


n_cols = 384
n_rows = 288

data_rate_with_margin = 181747.2 #bits per second
data_rate_science = (data_rate_with_margin / 1.2) - 500.*8. #bits per second


sections = {
    "a":{"rows":np.arange(0, 21), "rgb":(50/256, 168/256, 82/256), "colour":"green"},
    "d":{"rows":np.arange(21, 21+88), "rgb":(31/256, 168/256, 209/256), "colour":"blue"},
    "b":{"rows":np.arange(21+88, n_rows), "rgb":(222/256, 129/256, 16/256), "colour":"orange"},
}


it = {"2d":2.54, "4d":0.91, "1n":13.35, "2n":13.35, "3n":13.35}
it_rhythm = {"2d":3., "4d":1., "1n":15., "2n":15., "3n":15.}


bins_d = {}


n_bits_all = [16]

for n_bits in n_bits_all:

    data_bins = int(np.floor(data_rate_science / n_cols / n_bits))
    print(n_bits, data_bins)
    
    
    if n_bits == 16:
        bins = [np.arange(0, 5),np.arange(5, 10),np.arange(10, 15),np.arange(15, 20),np.arange(60, 65)]
    elif n_bits == 20:
        bins = [np.arange(0, 5),np.arange(5, 10),np.arange(10, 15),np.arange(15, 20),np.arange(60, 65)]
    elif n_bits == 24:
        bins = [np.arange(0, 7),np.arange(7, 14),np.arange(14, 21),np.arange(60, 67)]

    n_bins_section2 = data_bins - len(bins)
    bin_rows_section2 = int(np.floor((n_rows - 21 - 88 - 1) / n_bins_section2))
    n_bins_unused = (n_rows - 21 - 88 - 1) - (n_bins_section2*bin_rows_section2)
    
    for i in range(21+88+n_bins_unused, n_rows-bin_rows_section2+n_bins_unused, bin_rows_section2):
        bins.append(np.arange(i, i+bin_rows_section2))
    bins_d[n_bits] = bins





array = np.zeros((n_rows, n_cols, 3))
for section, section_d in sections.items():
    rows = section_d["rows"]
    array[rows, :, :] = section_d["rgb"]


for bits, bins in bins_d.items():
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    im = ax.imshow(array, aspect=1, alpha=1.0)
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    ax.set_xlabel("Detector column (spectral direction)")
    ax.set_ylabel("Detector row (spatial direction)")
    ax.set_title("VenSpec-H detector binning with dark row (typical data rate; %i bit readout per pixel)" %bits)
    
    for i, bin_ in enumerate(bins):
        # print(bin_[0])
        rect = patches.Rectangle((0, bin_[0]), 384, bin_[-1]-bin_[0]+1, linewidth=1, edgecolor='r', facecolor='none')
    
        ax.add_patch(rect)
        if len(bin_) > 5:
            ax.text(192, bin_[-1], "Detector bin %i" %(i+1))
        else:
            ax.text(192, bin_[-1]+1, "Detector bin %i" %(i+1), fontsize=8)

    plt.savefig("venspech_detector_binning_%ibits.png" %bits, dpi=300)

