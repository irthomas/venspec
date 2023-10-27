# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 13:44:13 2022

@author: iant

Blaze function
"""
import numpy as np

def F_blaze(x, blazef, blazew):
    """wavenumber grid, wavenumber peak, FSR wavenumber"""
    
    dx = x - blazef
    F = np.sinc((dx) / blazew)**2
    return F
