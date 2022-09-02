# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 12:17:58 2022

@author: iant

SPECTROMETER CALCULATIONS
"""

import numpy as np
# import scipy.constants as spc
import matplotlib.pyplot as plt

from instrument.snr.hr_nadir_spectra import hr_nadir_spectra

def gaussian(x, a, b, c, d):
    return a * np.exp(-((x - b)/c)**2.0) + d


band_dict = {
    "1":{"daynight":["n"], "um_range_aim":[1.16, 1.18], "n_px":222},
    "2a":{"daynight":["d", "n"], "um_range_aim":[2.34, 2.422], "n_px":36},
    "2b":{"daynight":["d", "n"], "um_range_aim":[2.422, 2.51], "n_px":186},
    "3":{"daynight":["n"], "um_range_aim":[1.704, 1.747], "n_px":222},
    "4":{"daynight":["d"], "um_range_aim":[1.367, 1.394], "n_px":222},
         }


# n_detector_rows = 288
n_detector_pixels = 384


fsr = 144.88 #cm-1
resolving_power = 11000.
slit_width = 1. #pixels

# resolving_power = 5500.
# slit_width = 2. #pixels


nadir = hr_nadir_spectra(1)


def spectral_cal():

    for band in band_dict.keys():
        
        #desired spectral range
        um_range_aim = band_dict[band]["um_range_aim"]
        
        #find corresponding diffraction order
        #order = d(sin alpha + sin beta)/lambda
        order_calc = (10000. / np.mean(um_range_aim)) / fsr
        # order_calc = gr_d_sin_alpha_sin_beta/np.mean(um_range_aim)
        order = int(np.round(order_calc))
        # print(np.mean(um_range_aim), "um", order_calc, "->", order)
        
        #wavelength of centre of diffraction order
        # band_centre_calc = gr_d_sin_alpha_sin_beta / order_calc
        band_centre = (10000. / order) / fsr
        # print("centre um=", band_centre)
        
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
        
        #spectral resolution: width of gaussian encompassing wavelengths received by a pixel * 2.355
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
        
        # print(band_dict[band]["detector_um_range"][0], band_dict[band]["order_range_um"][0], "->", band_dict[band]["band_um_range"][0])
        # print(band_dict[band]["detector_um_range"][1], band_dict[band]["order_range_um"][1], "->", band_dict[band]["band_um_range"][1])


    
        #centre wavelengths on each pixel even if outside FSR
        detector_px_centres = np.arange(band_dict[band]["detector_um_range"][0], band_dict[band]["detector_um_range"][-1]-0.000001, band_dict[band]["px_sampling_um"])
    
        #centre wavelengths on each pixel (real detector)
        detector_ixs = np.where((detector_px_centres >= band_um_range[0]) & (detector_px_centres <= band_um_range[1]))
        band_dict[band]["detector_ums"] = detector_px_centres[detector_ixs]
        band_dict[band]["detector_px_centres"] = detector_px_centres
        # band_dict[band]["detector_ums"] = np.arange(band_um_range[0], band_um_range[1], band_dict[band]["px_sampling_um"])
        print("Band", band, "pixels used:", len(band_dict[band]["detector_ums"]))
        
    return band_dict


def pixel_illumination_table():
    
    band_dict = spectral_cal()

    table = np.zeros((n_detector_pixels, len(band_dict.keys())))
    for band_ix, band in enumerate(band_dict.keys()):
        px_real = band_dict[band]["detector_ums"]
        for px_ix, px in enumerate(band_dict[band]["detector_px_centres"]):
            if px in px_real:
                table[px_ix, band_ix] = px
                
    return table
    


def plot_spectral_cal():

    band_dict = spectral_cal()
    

    for band in band_dict.keys():

        fig2, ax2 = plt.subplots(figsize=(14, 6))
        
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



if __name__ == "__main__":
    plot_spectral_cal()