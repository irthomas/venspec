# -*- coding: utf-8 -*-
"""
Created on Thu May 25 13:41:44 2023

@author: iant

GENERATE EXAMPLE VENSPEC-H DETECTOR FRAMES FROM NOMAD DATA

1) GET RANDOM NOMAD H5
2) EXTRACT RELEVANT SPECTRA, ONE PER VSH BIN
3) INTERPOLATE TO NEW FRAME SIZE
4) MAKE NON-ILLUMINATED REGION AND CUTOFF
5) OUTPUT FOR DIFFERENT NBITS AND BIN SIZES

"""

import h5py
import numpy as np
import struct


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FixedLocator

from tools.file.hdf5_functions import get_filepath, make_filelist


#code to check SO linescans
# import re
# h5fs, h5s, _ = make_filelist(re.compile(".*SO_._CL.*"), "hdf5_level_1p0a")

# for h5f, h5 in zip(h5fs, h5s):
#     print(h5)
#     orders = h5f["Channel/DiffractionOrder"][...]
#     print(min(orders), max(orders))

# stop()



#detector specs
# fovA = 34
# fovB = 194
# fovDark = 256 - fovA - fovB

n_cols = 384

saturation_level = 0.8

SAVE_FRAMES = True
# SAVE_FRAMES = False

for name in ["solar_cal_type1", "solar_cal_type2", "solar_cal_type3"]:



    #solar calibration
    if name == "solar_cal_type1":
        #type 1
        n_bits = 16
        n_bits_signed = n_bits - 1
        n_bins = [256]*100 #all pixels
        # n_frames_per_second = 1
    
    
    if name == "solar_cal_type2":
        #type 2
        n_bits = 16
        n_bits_signed = n_bits - 1
        n_bins = [40]*100 #all illuminated
        # n_frames_per_second = 10
    
    
    if name == "solar_cal_type3":
        #type 3
        n_bits = 16
        n_bits_signed = n_bits - 1
        n_bins = [34, 194]*50 #alternate top and bottom
        # n_frames_per_second = 10
    
    
    
    # h5 = "20221009_141204_1p0a_SO_1_CF" #fullscan slow
    # h5 = "20230119_115108_1p0a_SO_1_CL" #linescan fast order 164
    # h5 = "20230119_214056_1p0a_SO_1_CL"
    h5 = "20180828_223824_1p0a_SO_1_CL" #linescan slow order 167 - decent solar lines
    
    def get_nomad_spectra(h5):
        if h5 == "20221009_141204_1p0a_SO_1_CF":
            order = 169
            bins = [124, 125, 126, 127, 130, 131]
            
        elif h5 in ["20180828_223824_1p0a_SO_1_CL"]:
            order = 167
            bins = [127]
            min_sig = 180000.
        
        h5_filepath = get_filepath(h5)
        
        
        with h5py.File(h5_filepath, "r") as h5_f:
        
        
            orders_in = h5_f["Channel/DiffractionOrder"][...]
            bins_in = h5_f["Science/Bins"][:, 0]
            
            bin_ixs = [i for i,v in enumerate(bins_in) if v in bins]
            orders_ixs = [i for i,v in enumerate(orders_in) if v == order]
            
            frame_ixs = [x for x in bin_ixs if x in orders_ixs]
            
            y = h5_f["Science/Y"][frame_ixs, 50:] #ignore low signal on first pixels
            
            # plt.figure(); plt.plot(np.mean(y, axis=1))
            sig_ixs = np.where(np.mean(y, axis=1) > min_sig)[0]
            
            y = y[sig_ixs, :]
        
        
        
        #normalise to the max
        y_max = np.max(y, axis=1)
        
        #add some random variability
        y_rand = np.random.rand(len(y_max)) / 200.0
        
        y = ((y.T / y_max) + y_rand).T
        
        return y
    
    
    y = get_nomad_spectra(h5)
    
    if np.mean(n_bins) > 100:
        y = np.vstack((y, np.flipud(y), y, np.flipud(y))) #if fullscan, copy a few times
    # stop()
    #expand to fill 384 pixels
    xold = np.arange(y.shape[1])
    xnew = np.arange(n_cols)
    
    
    ynew = np.zeros((y.shape[0], n_cols))
    for i in range(y.shape[0]):
        ynew[i, :] = np.interp(xnew, xold, y[i, :])
    
    #add noise
    rand_array = (np.random.rand(ynew.shape[0], ynew.shape[1]) - 0.5) / 1000.0
    ynew += rand_array
    
    #normalise to saturation level
    ynew *= 2**n_bits_signed * saturation_level
    
    ynew = ynew.astype(int)
    
    
    n_frames = int(np.floor(ynew.shape[0]/np.mean(n_bins)))
    
    
    #gradually increase and decrease frame levels like a solar calibration
    frame_scalar = np.polyval(np.polyfit([0, n_frames/2, n_frames], [0.8, 1.0, 0.8], 2), np.arange(n_frames))
    
    # plt.figure()
    # plt.plot(frame_scalar)
    
    print("name")
    
    start_row = 0
    with PdfPages("venspech_%s_%icols_%ibits_frames.pdf" %(name, n_cols, n_bits)) as pdf: #open pdf
        for frame_ix in np.arange(n_frames):
            end_row = start_row + n_bins[frame_ix] + 1
            if end_row > ynew.shape[0]:
                continue
            print(frame_ix, start_row, end_row)
            
            y_frame = (ynew[start_row:end_row, :] * frame_scalar[frame_ix]).astype(int)
            
            fig, ax = plt.subplots(figsize=(18, 10))
            
            
            im = ax.imshow(y_frame, aspect="auto")
            ax.xaxis.set_minor_locator(FixedLocator(np.arange(0, y_frame.shape[1], 1.0)))
            ax.yaxis.set_minor_locator(FixedLocator(np.arange(0, y_frame.shape[0], 1.0)))
            ax.grid(which="minor", color='k', linewidth=0.1)
            ax.grid(which="major", color='k', linewidth=0.1)
            # plt.minorticks_on()
            plt.colorbar(im)
            ax.set_xlabel("Spectral pixel number")
            ax.set_ylabel("Bin number")
            
            #convert to bytes and write to file
            
            name_out = "venspech_%s_%irows_%icols_%ibits_frame%03i.dat" %(name, n_bins[frame_ix], n_cols, n_bits, frame_ix)
            ax.set_title(name_out)
            
            if SAVE_FRAMES:
                with open(name_out, "wb") as f:
                    for i in range(y_frame.shape[0]):
                        out = struct.pack('>'+"h"*n_cols, *y_frame[i, :])
                        f.write(out)
                    
                #read in file to check all is good
                with open(name_out, "rb") as f:
                    bytes_ = f.read()
                    
                
                y_in = np.zeros_like(y_frame)
                for i in range(y_frame.shape[0]):
                    
                    start_byte = int(i*(n_bits/8)*n_cols)
                    end_byte = int((i+1)*(n_bits/8)*n_cols)
                    
                    out = struct.unpack('>'+"h"*n_cols, bytes_[start_byte:end_byte])
                    y_in[i, :] = out
                    
                if not np.array_equal(y_frame, y_in):
                    print("Error: frame %i failed" %frame_ix)
            
            pdf.savefig()
            plt.close()
            
            start_row = end_row
    
