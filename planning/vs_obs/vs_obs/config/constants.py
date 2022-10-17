# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:27:02 2022

@author: iant

"""

import spiceypy as sp


#SPICE constants
SPICE_FORMATSTR = "C"
SPICE_PRECISION = 3
SPICE_DATETIME_FMT = "%Y %b %d %H:%M:%S.%f"
# SHORT_DATETIME_FMT = "%d %b %Y"

DATETIME_FMT_SECONDS = "%Y %b %d %H:%M:%S"

SPICE_METHOD = "INTERCEPT/ELLIPSOID"
SPICE_SHAPE_MODEL_METHOD = "Ellipsoid"
SPICE_ABCORR = "NONE"
# SPICE_OBSRVR = "-668"

SPICE_LONGITUDE_FORM = "PLANETOCENTRIC"
SPICE_PLANET_ID = 299 # venus

SPICE_TARGET = "VENUS"
SPICE_PLANET_REFERENCE_FRAME = "IAU_VENUS"

SPICE_DREF = "ENVISION_VENSPEC_H"

SP_DPR = sp.dpr()

CHANNEL_ID = sp.bods2c(SPICE_DREF) #find channel id number
SPICE_VENUS_AXES = sp.bodvrd("VENUS", "RADII", 3)[1] #no ellipsoid for Venus -> circular so only one radius value
SPICE_VENUS_RADIUS = SPICE_VENUS_AXES[0]



#instrument constants
DETECTOR_READOUT_TIME = 0.5
STANDBY_DURATION = 10.
PRECOOLING_DURATION = 600.


WHEEL_ROTATION_SPEED = 60.0 #degrees per second
FILTER_WHEEL_LAYOUT = ["2v", "2h", "2", "dk", "4", "4h", "1", "3"]
ANGLE_BETWEEN_FILTERS = 360. / len(FILTER_WHEEL_LAYOUT)



SPICE_OBSERVER = "-668"