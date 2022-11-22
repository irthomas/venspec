# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 13:50:46 2022

@author: iant

DIFFRACTION GRATING
"""


import numpy as np

import shapely
from shapely.geometry import LineString, Point, Polygon
from shapely import affinity

import matplotlib.pyplot as plt


def azimuth(point1, point2):
    '''azimuth between 2 shapely points (interval 0 - 360)'''
    angle = np.arctan2(point2.x - point1.x, point2.y - point1.y)
    return np.degrees(angle) if angle >= 0 else np.degrees(angle) + 360

def reflection(startP, length, angle):
    end = Point(startP.x - length, startP.y)
    line = LineString([startP, end])
    line = affinity.rotate(line, angle, origin=startP)   

    return line        

def dist(P1, P2):
    
    if isinstance(P1, Point):
        x1 = P1.x
        y1 = P1.y
    else:
        x1 = P1[0]
        y1 = P1[1]
         
    if isinstance(P2, Point):
        x2 = P2.x
        y2 = P2.y
    else:
        x2 = P2[0]
        y2 = P2[1]
    
    dist = np.sqrt( (x2 - x1)**2 + (y2 - y1)**2 )
    return dist

def closest_point(point, points):
    
    dists = [dist(i, point) for i in points]
    idx = np.array(dists).argmin()
    print(idx)
    return Point(points[int(idx)])


#slit origin 0,0
scx = 0.
scy = 0.

#grating centre A,0
gcx = 200. #mm
gcy = 0. #mm
gl = 200. ##mm
ga = 30. #degrees

#line start
lp0 = Point(scx, scy)
lp1 = Point(gcx+50, gcy)





fig, ax = plt.subplots(figsize=(14, 8))
ax.set_aspect("equal")


#grating is a mirror
gp0 = Point(gcx, -gl/2)
gp1 = Point(gcx, gl/2)

grating = LineString([gp0, gp1])
grating = affinity.rotate(grating, ga, 'center')
# plt.plot(*np.array(grating.coords).T)

gpoints = [[gcx + 5, -gl/2]]
for i, y in enumerate(np.arange(-gl/2, gl/2+3, 3)):
    if i % 2 == 0:
        gpoints.append([gcx - 0.5, y])
    else:
        gpoints.append([gcx + 0.5, y])
gpoints.append([gcx + 5, gl/2])

grating_poly = Polygon(gpoints)
grating_poly = affinity.rotate(grating_poly, ga, origin="center")

plt.plot(*grating_poly.exterior.xy)




for sca in np.arange(-5, 6)[:]:
    print(sca)
    line1 = LineString([lp0, lp1])
    line1 = affinity.rotate(line1, sca, origin=lp0)
    plt.plot(*np.array(line1.coords).T)
    


    
    line1grating = line1.intersection(grating_poly)
    plt.plot(*np.array(line1grating.coords).T)
    line1gratingP = Point(line1grating.coords)
    plt.scatter(*np.array(line1gratingP))
    

    if sca == 0.0:
        centre_line = line1
        centre_lgP = line1gratingP

    centre_rangle = ga*2.0


    
    P1 = Point(lp0)
    # P2 = Point(np.array(np.array(grating.coords))[0, :])
    P2 = closest_point(line1gratingP, grating_poly.exterior.coords)

    if P2.y > line1gratingP.y:
        angle = 2.* (90.0 -(azimuth(line1gratingP, P2) - azimuth(line1gratingP, P1)))
    # else:
    #     angle = azimuth(P2, line1gratingP) - azimuth(P1, line1gratingP)
        
    
    
    # angle = 20
    line2 = reflection(line1gratingP, 100, angle)
    
    plt.scatter(*np.array(P1.coords)[0])
    plt.scatter(*np.array(P2.coords)[0])

    # ax.set_xlim(200, 220)
    # ax.set_ylim(-30, -10)
    
    
    for gra in np.arange(-5, 5):
        line3 = reflection(line1gratingP, 100, centre_rangle + gra)
        plt.plot(*np.array(line3.coords).T, "--")
        
centre_rline = reflection(centre_lgP, 100, centre_rangle)
plt.plot(*np.array(centre_rline.coords).T, "k--")

