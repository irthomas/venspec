# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 16:06:05 2022

@author: iant


https://pds-geosciences.wustl.edu/mgn/mgn-v-rss-5-gravity-l2-v1/mg_5201/images/topogrd.lbl

  LINES                        = 180
  LINE_SAMPLES                 = 360
  SAMPLE_TYPE                  = UNSIGNED_INTEGER 
  SAMPLE_BITS                  = 8
  UNIT                         = "KILOMETER"
  OFFSET                       = 2.492
  SCALING_FACTOR               = 0.0478
  DESCRIPTION                  = "A byte-scaled image showing the
  topography of Venus.  The following algorithm translates topography
  to an image DN:
         DN = (TOPOGRAPHY <KILOMETERS> + 2.492) / 0.0478"

DN * 0.0478 - 2.492

"""

import os
import numpy as np

from tools.file.paths import paths


def make_venus_magellan_map():

    import matplotlib.pyplot as plt

    venus_map = np.zeros((180, 360))
    with open(os.path.join(paths["REFERENCE_DIRECTORY"], "venus_topogrd.img"), "rb") as f:
        bytes_ = f.read()
        
    for i, byte in enumerate(bytes_):
        venus_map[int(np.floor(i/360.)), np.mod(i, 360)] = np.float32(byte) * 0.0478 - 2.492
    
    np.savetxt(os.path.join(paths["REFERENCE_DIRECTORY"], "venus_magellan_map.txt"), venus_map, fmt="%0.3f")    

    im = plt.imshow(venus_map, extent=(-180, 180, -90, 90))
    plt.colorbar(im)  #in kilometers



def get_venus_magellan_data():
    
    path = os.path.join(paths["REFERENCE_DIRECTORY"], "venus_magellan_map.txt")
    
    return np.loadtxt(path)
