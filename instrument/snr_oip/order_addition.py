# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 12:59:07 2023

@author: iant


READ IN OIP INPUTS:
    DO SPECTRAL CALIBRATION
    CHECK ORDER ADDITION
    IN-BAND STRAYLIGHT

"""


import os
import json
import numpy as np
import matplotlib.pyplot as plt



from instrument.snr_oip.config import MODEL_DATA, BAND_DATA, MODEL_RESULTS

band = "1_NIGHT"
# band = "2A_NIGHT"
# band = "2B_NIGHT"


RELOAD = True


def F_sinc2(x, a, b): #inverted sinc2
    a = a + 0.0001 #fudge to stop infinity at peak
    return ((np.sin((x - a) / b * 2.783)**2.0) / (((x - a) / b * 2.783)**2.0))

def gaussian(x, a, b, c, d):
    #amplitude is a variable (a)
    return a * np.exp(-(x - b)**2.0/(2.0 * c**2.0)) + d




def convolve_isrf(px_nm, dpx_nm, isrf_gaussian_width, venus_nm, venus_spectrum, ax=None):

    #convert spectral resolution (FWHM) to gaussian width (sigma)
    gauss_width = dpx_nm / 2.355
    
    #calculate wavelength range of gaussian and make grid
    nm_grid_start = -1.0 * gauss_width * isrf_gaussian_width
    nm_grid_end = gauss_width * isrf_gaussian_width
    
    #find indices in hr grid covering +- wavelength range of gaussian
    nu_hr_start_ix = np.searchsorted(venus_nm, px_nm + nm_grid_start) #start index
    nu_hr_end_ix = np.searchsorted(venus_nm, px_nm + nm_grid_end) #end index
    
    if nu_hr_start_ix == nu_hr_end_ix:
        #pixel ILS is beyond the range of the input spectrum -> using nearest value instead
        px_venus_convolved = venus_spectrum[-1]

        warning = True
        
    else:
        
        px_grid = venus_nm[nu_hr_start_ix:nu_hr_end_ix] - px_nm
        px_venus_signal = venus_spectrum[nu_hr_start_ix:nu_hr_end_ix]
        
        isrf_gauss = gaussian(px_grid, 1.0, 0.0, gauss_width, 0.0)
        isrf_venus_signal = px_venus_signal * isrf_gauss
        
        px_venus_convolved = np.sum(isrf_venus_signal) / np.sum(isrf_gauss)
        
        warning = False

        # test ILS
        if ax:
            ax[0].plot(px_grid + px_nm, isrf_gauss, "k")
            ax[1].plot(px_grid + px_nm, px_venus_signal, "orange")
            # ax[1].plot(px_grid + lambda_um, ils_venus_signal)
            # ax[1].plot([px_grid[0] + lambda_um, px_grid[-1] + lambda_um], [px_venus_convolved, px_venus_convolved])  
        # stop()
    
    return px_venus_convolved, warning





def detector_isrf(pxs_nm, dpxs_nm, isrf_gaussian_width, venus_nm, venus_spectrum, plot_isrf=False):
    
    venus_px = np.zeros_like(pxs_nm)
    
    warning_shown = False
    if plot_isrf:
        fig_ils, ax_ils = plt.subplots(figsize=(17, 5), constrained_layout=True)
        ax_ils2 = ax_ils.twinx()
        ax_ils.set_title("ISRF vs radiance")
        ax_ils.set_xlabel("Wavelength (nm)")
        ax_ils.set_ylabel("Normalised ISRF")
        ax_ils2.set_yscale("log")
        ax_ils2.set_ylabel("Venus radiance (log scale)")
        ax = [ax_ils, ax_ils2]
    for px_ix, (px_nm, dpx_nm) in enumerate(zip(pxs_nm, dpxs_nm)):
        if plot_isrf:
            venus_px[px_ix], warning = convolve_isrf(px_nm, dpx_nm, isrf_gaussian_width, venus_nm, venus_spectrum, ax=ax)
        else:
            venus_px[px_ix], warning = convolve_isrf(px_nm, dpx_nm, isrf_gaussian_width, venus_nm, venus_spectrum)

        #display error just once
        if warning and not warning_shown:
            print("Warning: pixel ISRF is beyond the range of the input spectrum -> using nearest value instead")
            warning_shown = True

    
    return venus_px








#read in night spectrum
if "night_spectrum" not in globals() or RELOAD:
    night_nm, night_spectrum = np.loadtxt(os.path.join(MODEL_DATA, "ENVIS-VS-VEH-SN0016-iss0rev0+AER-NightSpectrum_20221116.dat"), skiprows=6, unpack=True)
    #reverse order for searching
    night_nm = np.append(night_nm[::-1], 2701.0)
    night_spectrum = np.append(night_spectrum[::-1], 0.0)
    


#read in beta angles
beta_angles = np.loadtxt(os.path.join(MODEL_DATA, "beta.dat"), skiprows=1)

#read in band cutoffs
filter_path = BAND_DATA[band]["data"]
band_nm, band_trans = np.loadtxt(filter_path, unpack=True)



# # plot all band grating transmittances
# plt.figure()
# for band in BAND_DATA.keys():
#     filter_path = BAND_DATA[band]["data"]
#     band_nm, band_trans = np.loadtxt(filter_path, unpack=True)
#     plt.plot(band_nm, band_trans, label="Band %s" %band)   
#     plt.title("Spectrometer section transmittances")
#     plt.xlabel("Wavelength (nm)")
#     plt.ylabel("Transmittance")
# plt.grid()
# plt.legend()



# convert to cm-1 for equal widths
# plt.figure()
# for band in BAND_DATA.keys():
#     filter_path = BAND_DATA[band]["data"]
#     band_nm, band_trans = np.loadtxt(filter_path, unpack=True)
#     x = 10000000.0/band_nm
#     x_centre_ix = np.argmax(band_trans)
#     plt.plot(x - x[x_centre_ix], band_trans, label="Band %s" %band, alpha=0.7)      
#     plt.title("Spectrometer section transmittances")
#     plt.xlabel("Wavenumber normalised to peak (cm-1)")
#     plt.ylabel("Transmittance")
# plt.grid()

# #plot sinc in wavenumber space
# # x = 10000000.0 / np.arange(-250, 250, 0.1)
# a = 0.0
# b = 200
# sinc2 = F_sinc2(x - x[x_centre_ix], a, b)
# plt.plot(x - x[x_centre_ix], sinc2 * 0.4, "k--", label="Sinc squared function")
# plt.legend()
# plt.xlim((-750, 750))

# stop()




#read in config file
conf_path = BAND_DATA[band]["conf"]
with open(conf_path, "r") as f:
    config = json.load(f)

min_wavel = config["MinWave"]
max_wavel = config["MaxWave"]

alpha_angle = config["AlphaAngle"]
grooves = config["Grooves"]
centre_order = config["Order"]

orders_d = {}
orders = np.arange(centre_order-1, centre_order+2)
# orders = np.arange(centre_order, centre_order+1)

filteroff_ixs = np.searchsorted(night_nm, [min_wavel, max_wavel])
spectrum_filteroff = np.copy(night_spectrum)
spectrum_filteroff[:filteroff_ixs[0]] = 0.0
spectrum_filteroff[filteroff_ixs[1]:] = 0.0

for order in orders:
    #mλ= d (sinα + sinβ) 
    # λ= d (sinα + sinβ)/m
    
    wavels = grooves * (np.sin(alpha_angle) + np.sin(beta_angles[:, 1])) / order * 1.0e9
    
    #ISRF convolution
    ixs = np.searchsorted(night_nm, wavels)
    isrf_conv = detector_isrf(wavels, wavels/11000.0, 3.0, night_nm, night_spectrum, plot_isrf=False)
    #ISRF convolution with filter
    isrf_filter_conv = detector_isrf(wavels, wavels/11000.0, 3.0, night_nm, spectrum_filteroff, plot_isrf=False)


    orders_d[order] = {"wavels":wavels, "isrf_conv":isrf_conv, "isrf_filter_conv":isrf_filter_conv}

#get the band transmittance and copy to each order
for order in orders:
    
    d_wavel = np.mean(orders_d[order]["wavels"] - orders_d[centre_order]["wavels"])
    orders_d[order]["band_nm"] = band_nm + d_wavel
    orders_d[order]["band_trans"] = band_trans
    
    #interpolate band onto detector wavelengths and multiply together
    band_trans_det = np.interp(orders_d[order]["wavels"], orders_d[order]["band_nm"], orders_d[order]["band_trans"])
    band_trans_isrf = orders_d[order]["isrf_filter_conv"] * band_trans_det
    orders_d[order]["band_trans_det"] = band_trans_det
    orders_d[order]["band_trans_isrf"] = band_trans_isrf
    
    # plt.figure()
    # plt.plot(orders_d[order]["band_nm"], orders_d[order]["band_trans"])
    # plt.plot(orders_d[order]["wavels"], band_trans_det)


#plot raw radiances
fig1, ax1a = plt.subplots()
ax1b = ax1a.twinx()
ax1a.plot(band_nm, band_trans, "k--")
# ax1a.set_yscale("log")
fig1.suptitle("Spectrometer section transmittance band %s" %band)
ax1a.set_xlabel("Wavelength (nm)")
ax1a.set_ylabel("Transmittance")
ax1b.set_ylabel("Radiance (grey=unconvolved,\ncoloured=unconvolved to detector ISRF and filter)")
ax1a.grid()

conv_sum = np.zeros_like(beta_angles[:, 1]) #sum of order addition

for i, order in enumerate(orders_d.keys()):
    
    colour = "C%i" %i
    
    wavels = orders_d[order]["wavels"]
    ax1a.fill_between(wavels, y1=np.zeros_like(wavels), y2=np.zeros_like(wavels)+0.4, alpha=0.3, label="Detector wavelength range for order %i" %order, color=colour)
    ixs = np.searchsorted(night_nm, wavels)
    # ax1b.plot(night_nm[ixs], night_spectrum[ixs], color=colour)
    
    # ax1b.plot(wavels, orders_d[order]["isrf_conv"])
    ax1b.plot(wavels, orders_d[order]["isrf_filter_conv"], color=colour)
    # ax1b.plot(wavels, orders_d[order]["band_trans_isrf"], color=colour)
    ax1a.plot(orders_d[order]["band_nm"], orders_d[order]["band_trans"], color=colour, label="Spectrometer section transmittance order %i" %order)
    
    # plt.figure()
    # plt.plot(wavels, orders_d[order]["isrf_filter_conv"], label="isrf_filter_conv")
    # plt.plot(wavels, orders_d[order]["band_trans_isrf"], label="band_trans_isrf")
    # plt.plot(wavels, orders_d[order]["band_trans_det"], label="band_trans_det")
    # plt.plot(wavels, orders_d[order]["band_trans_isrf"]/orders_d[order]["isrf_filter_conv"], label="band_trans_isrf/isrf_filter_conv")
    # plt.legend()
    
    
    conv_sum += orders_d[order]["band_trans_isrf"]
    
ax1a.legend()

ax1a.axvline(min_wavel, color="k", linestyle=":")
ax1a.axvline(max_wavel, color="k", linestyle=":")

ax1b.plot(night_nm, night_spectrum, "k", alpha=0.3)
ax1a.set_xlim((1100, 1240))

plt.figure()
plt.title("Radiance as measured by detector, with order addtion, ISRF and filter convolutions")
plt.xlabel("Pixel number")
plt.ylabel("Radiance after order addtion")
plt.plot(conv_sum)
plt.grid()





#band1 order addition: plot each order separately


#make band1 order addition plot






#plot sinc2 reduced radiances
# fig1, ax1a = plt.subplots()
# ax1b = ax1a.twinx()
# ax1a.plot(band_nm, band_trans, "k--")
# # ax1a.set_yscale("log")
# fig1.suptitle("TransmittanceSpec band %s" %band)
# ax1a.set_xlabel("Wavelength (nm)")
# ax1a.set_ylabel("Transmittance")
# ax1a.grid()


# #interpolate band onto input radiance and multiply together
# band_trans_hr = np.interp(night_nm, band_nm, band_trans)
# night_spec_band = night_spectrum * band_trans_hr

# for order in orders_d.keys():
#     wavels = orders_d[order]["wavels"]
#     ax1a.fill_between(wavels, y1=np.zeros_like(wavels), y2=np.zeros_like(wavels)+0.4, alpha=0.3, label="Detector wavelength range for order %i" %order)
#     ixs = np.searchsorted(night_nm, wavels)
#     ax1b.plot(night_nm[ixs], night_spec_band[ixs])

#     conv = detector_isrf(wavels, wavels/11000.0, 3.0, night_nm, night_spectrum, plot_isrf=False)
#     ax1b.plot(wavels, conv)

# ax1a.legend()

# ax1a.axvline(min_wavel, color="k", linestyle=":")
# ax1a.axvline(max_wavel, color="k", linestyle=":")

# ax1b.plot(night_nm, night_spec_band, "k", alpha=0.3)



res_dir = "out_d120.0_w268.1_c228.2"
res_out = np.loadtxt(os.path.join(MODEL_RESULTS, res_dir, "VEH%s.txt" %(band.replace("_", "-"))), skiprows=20)



# plot OIP SNR standalone
# snr = res_out[:, 6]
# snr_x = res_out[:, 1]
# plt.figure()
# plt.plot(snr_x, snr)
# plt.axvline(min_wavel, color="k", linestyle=":")
# plt.axvline(max_wavel, color="k", linestyle=":")
# plt.grid()
# plt.title("Band %s SNR (%s)" %(band, res_dir))
# plt.xlabel("Wavelength (nm)")
# plt.ylabel("SNR for one pixel")



# # plot OIP SNR scaled to radiances
# snr = res_out[:, 6]
# snr_x = res_out[:, 1]
# ax1b.plot(snr_x, snr/700.0)





#plot resolving power
# plt.figure()
# plt.plot(res_out[:, 1], res_out[:, 1]/res_out[:,2])
# plt.xlabel("Wavelength nm")
# plt.ylabel("Resolving power")
# plt.title("Band %s" %band)

# plot all sinc2 to see if shape is the same





