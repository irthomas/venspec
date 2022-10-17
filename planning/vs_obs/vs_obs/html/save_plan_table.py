# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 11:07:32 2022

@author: iant
"""


def save_plan_table(h_rows, filename):
    
    h = "<!DOCTYPE html><html><body>\n"
    h += "<h1>Operations plan</h1>\n"
    h += "<table border=1><tr><th>Time</th><th>Sequence</th></tr>\n"
    h += h_rows
    h += "</table>\n"
    
    h += "</body></html>\n"
    
    with open("%s.html" %filename, "w") as f:
        f.write(h)
    
