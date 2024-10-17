# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 10:17:26 2024

@author: iant

ALTITUDE ABOVE SURFACE
"""

import spiceypy as sp

et = 1117896296.9732976
et = 1117899165.535671

subpnt = sp.subpnt('NEAR POINT/ELLIPSOID', 'VENUS', et, 'IAU_VENUS', 'NONE', '-668')

dist = sp.vnorm(subpnt[2])

print(dist)
