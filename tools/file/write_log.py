# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 12:07:41 2018

@author: iant


WRITE_LOG

"""
import os

from tools.file.paths import paths




def write_log(line_to_write, filename, path=""):
    """function to append log file"""
    if path == "":
        log_file = open(os.path.join(paths["BASE_DIRECTORY"], filename), 'a')
    else:
        log_file = open(os.path.join(path, filename), 'a')

    log_file.write(line_to_write+'\n')
    log_file.close()
