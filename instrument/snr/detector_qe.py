# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 13:43:40 2022

@author: iant

DETETCTOR QE VS WAVELENGTH
"""


import numpy as np
# import matplotlib.pyplot as plt


def detector_qe():
    """get QE versus wavelength (microns)"""
    qe_plot = np.array([[0.5, 0.], [0.75, 0.], [0.78, 0.1], [0.87, 0.75], [2.27, 0.75], [2.5, 0.69], [2.64, 0.02], [2.66, 0.], [3., 0.]])
    qe_um = qe_plot[:, 0]
    qe = qe_plot[:, 1]
    # plt.plot(qe_plot[:, 0], qe_plot[:, 1])
    return qe_um, qe



