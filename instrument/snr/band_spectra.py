# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 12:41:23 2022

@author: iant
"""

import numpy as np
# import scipy.constants as spc
import matplotlib.pyplot as plt



filter_dict = {
    "1":{"daynight":["n"], "um":[[1.1, 0.0], [1.16, 0.0], [1.16, 1.0], [1.18, 1.0], [1.18, 0.0], [2.9, 0.0]], "n_px":194+34},
    "2":{"daynight":["d", "n"], "um":[[1.1, 0.0], [2.34, 0.0], [2.34, 1.0], [2.51, 1.0], [2.51, 0.0], [2.9, 0.0]], "n_px":34},
    "3":{"daynight":["n"], "um":[[1.1, 0.0], [1.704, 0.0], [1.704, 1.0], [1.747, 1.0], [1.747, 0.0], [2.9, 0.0]], "n_px":194+34},
    "4":{"daynight":["d"], "um":[[1.1, 0.0], [1.367, 0.0], [1.367, 1.0], [1.394, 1.0], [1.394, 0.0], [2.9, 0.0]], "n_px":194+34},
         }

slit_dict = {
    "2a":{"um":[[1.1, 1.0], [2.422, 1.0], [2.422, 0.0], [2.9, 0.0]]},
    "2b":{"um":[[1.1, 1.0], [1.747, 1.0], [1.747, 0.0], [2.422, 0.0], [2.422, 1.0], [2.9, 1.0]]},
    }

combined_dict = {
    "1":{"daynight":["n"], "um":[[1.1, 0.0], [1.16, 0.0], [1.16, 1.0], [1.18, 1.0], [1.18, 0.0], [2.9, 0.0]], "n_px":194+34},
    "2a":{"daynight":["d", "n"], "um":[[1.1, 0.0], [2.34, 0.0], [2.34, 1.0], [2.422, 1.0], [2.422, 0.0], [2.9, 0.0]], "n_px":34},
    "2b":{"daynight":["d", "n"], "um":[[1.1, 0.0], [2.422, 0.0], [2.422, 1.0], [2.51, 1.0], [2.51, 0.0], [2.9, 0.0]], "n_px":34},
    "3":{"daynight":["n"], "um":[[1.1, 0.0], [1.704, 0.0], [1.704, 1.0], [1.747, 1.0], [1.747, 0.0], [2.9, 0.0]], "n_px":194+34},
    "4":{"daynight":["d"], "um":[[1.1, 0.0], [1.367, 0.0], [1.367, 1.0], [1.394, 1.0], [1.394, 0.0], [2.9, 0.0]], "n_px":194+34},
         }


def plot_filters():

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
    
        
    # plt.grid()
    ax1a.legend(loc="upper right")
    ax1b.legend(loc="upper right")
    ax1c.legend(loc="upper right")
    ax1c.set_xlabel("Wavelength um")
    ax1a.set_ylabel("Filter transmission on/off")
    ax1b.set_ylabel("Filter transmission on/off")
    ax1c.set_ylabel("Filter transmission on/off")
    
    fig1.suptitle("VenSpec-H filter wheel and slit spectra")
    fig1.savefig("bands.png")
    
