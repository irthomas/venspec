# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 15:01:59 2023

@author: iant

ANALYSE OUTPUT OF OIP SNR MODEL
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from tools.general.get_nearest_index import get_nearest_index

from instrument.snr_oip.config import BAND_DATA, MODEL_RESULTS, MODEL_EXE, MODEL_ROOT, MODEL_OUTPUT

bands = BAND_DATA.keys()

temps_detector = [118, 120, 135] #K
temps_warmsec = [273.15] #K
temps_coldsec = [228.15] #K


res_d = {}
for temp_detector in temps_detector:
    for temp_warmsec in temps_warmsec:
        for temp_coldsec in temps_coldsec:

            #get output directory to reflect temperatures
            dir_name = "out_d%0.1f_w%0.1f_c%0.1f" %(temp_detector, temp_warmsec, temp_coldsec)
            dir_path = os.path.join(MODEL_RESULTS, dir_name)

            res_d[dir_name] = {}
            for band in bands:
                #read in output file if measured
                if BAND_DATA[band]["measured"]:
                    res_path = os.path.join(dir_path, BAND_DATA[band]["snr"])
                    with open(res_path, "r") as f:
                        lines = f.readlines()
                        
                    res_d[dir_name][band] = {"px":[], "nm":[], "d_nm":[], "bol_sat":[], "bol_snr":[], "eol_sat":[], "eol_snr":[]}
                    for line in lines:
                        if line[0] == "#":
                            continue
                        else:
                            split = line.split("\t")
                            for i, key in enumerate(res_d[dir_name][band].keys()):
                                res_d[dir_name][band][key].append(float(split[i]))

band = "1_NIGHT"
# band = "2A_NIGHT"
# band = "2B_NIGHT"
# band = "3_NIGHT"
ref_nm = BAND_DATA[band]["ref_wavel"]*1000.0
n_detector_rows = BAND_DATA[band]["rows"]
                            
plt.figure(figsize=(10, 6))           
for key in res_d.keys():
    
    td = float(key.split("_")[1][1:])
    tw = float(key.split("_")[2][1:])
    tc = float(key.split("_")[3][1:])
    
    snr_1px = np.asfarray(res_d[key][band]["eol_snr"])
    snr_allpx = snr_1px * np.sqrt(n_detector_rows)
    
    plt.plot(res_d[key][band]["nm"], snr_allpx, label="DetT=%0.1fK, WarmT=%0.1fK, ColdT=%0.1fK" %(td, tw, tc), alpha=0.7)

    ix = get_nearest_index(ref_nm, res_d[key][band]["nm"])
    plt.scatter(res_d[key][band]["nm"][ix], snr_allpx[ix])
    print(key, "SNR=", snr_allpx[ix])
    
    plt.text(res_d[key][band]["nm"][ix]+1, snr_allpx[ix], "SNR=%0.1f" %snr_allpx[ix])

plt.legend()
plt.ylabel("SNR Binned Detector Rows")
plt.xlabel("Pixel wavelength (nm)")
plt.title("VenSpec-H SNRs Band %s" %band)
plt.grid()
plt.axvline(ref_nm, color="k", linestyle="--")


# res_d[key][band]["nm"]
    
            
                
            
            
