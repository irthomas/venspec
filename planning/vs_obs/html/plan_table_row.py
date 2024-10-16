# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 10:50:54 2022

@author: iant
"""

from datetime import datetime


from planning.vs_obs.config.constants import SPICE_DATETIME_FMT


def plan_table_row(dt, string):

    dt_str = datetime.strftime(dt, SPICE_DATETIME_FMT)
    h = "<tr><td>%s</td><td>%s</td></tr>\n" % (dt_str, string)

    return h
