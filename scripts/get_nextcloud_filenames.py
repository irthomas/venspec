# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 08:59:14 2024

@author: iant

RUN THROUGH NEXTCLOUD COLLECTING A LIST OF ALL FILENAMES
"""


import glob
import os


ROOT_DIR = r"C:\Users\iant\Nextcloud"

dirs = [
    "EnVision.General",
    "VenSpecH.CoworkingSpace",
    "VenSpecH.Internal",
    "VenSpecH.Released",
    "VenSpecH.Staging"]

exts = [".doc", ".docx", ".xls", ".xlsx", ".ppt", ".pptx", ".pdf"]

lines = []
for dir_ in dirs:
    path = os.path.join(ROOT_DIR, dir_)

    filepaths = sorted(glob.glob(path + os.sep + "**" + os.sep + "*.*", recursive=True))
    filedirs = [os.path.abspath(os.path.join(filepath, os.pardir)) for filepath in filepaths]

    filenames = [os.path.basename(filepath) for filepath in filepaths]

    extensions = [os.path.splitext(filename)[1] for filename in filenames]

    for i in range(len(filepaths)):
        if extensions[i] in exts:
            lines.append("%s\t%s\n" % (filenames[i], filedirs[i]))


lines = sorted(lines)

with open("Nextcloud_list.txt", "wb") as f:
    for line in lines:
        f.write(line.encode())
