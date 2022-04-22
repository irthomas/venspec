# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 21:12:50 2022

@author: iant
"""

import os
import paramiko
import glob
import numpy as np
from datetime import datetime

from tools.file.passwords import passwords

user = "iant"
password = passwords["hera"]
remote_root_path = "/ae/projects4/EnVision/DOCS_TO_CLOUD/VS/VEH/"
local_root_path = r"C:\Users\iant\Nextcloud\EnVision\DOCUMENTS\VS\VEH\*"


SYNC_LIMIT = 60.0 #seconds limit for checking last modified times
DT_FORMAT = "%Y-%m-%d %H:%M:%S"

def get_linux_tree_filenames(user, password, remote_path, host="hera.oma.be"):
    """Get a list of all h5 filenames present in HDF5 level subdirectories on a remote linux server. Search month by month""" 
    
    print("Connecting to %s" %host)
    p = paramiko.SSHClient()
    p.set_missing_host_key_policy(paramiko.AutoAddPolicy())   # This script doesn't work for me unless this line is added!
    p.connect(host, port=22, username=user, password=password, banner_timeout=30)
    
    filename_dict = {}
    
    # find_fmt = 'find %s/*.* -exec md5sum {} +'
    find_cmd = f'find "{remote_path}" -type f -printf "%p %T@ %s\\n"'
    
        
    
    stdin, stdout, stderr = p.exec_command(find_cmd)
    
    output = stdout.read().decode().split("\n")
    for s in output:
        if s != "":
            rsplit = s.rsplit(" ", 2)
            filename = os.path.basename(rsplit[0])
            dt_seconds = int(np.round(np.float32(rsplit[1]))) #seconds since Jan. 1, 1970, 00:00 GMT
            dt = datetime.utcfromtimestamp(dt_seconds)
            file_size = int(rsplit[2]) #bytes
            
            if filename in filename_dict.keys():
                print("Error: %s already found" %filename)
            else:
                filename_dict[filename] = {"dt_seconds":dt_seconds, "dt":dt, "size":file_size, "path":rsplit[0]}
                


    
    p.close()

    return filename_dict



def get_windows_tree_filenames(local_path):

    filename_dict = {}
    file_paths = sorted(glob.glob(local_path+"/**/*.*", recursive=True))
    for file_path in file_paths:
        filename = os.path.basename(file_path)
        dt_seconds = os.path.getmtime(file_path)
        dt = datetime.utcfromtimestamp(dt_seconds)
        file_size = os.path.getsize(file_path)
    
        if filename in filename_dict.keys():
            print("Error: %s already found" %filename)
        else:
            filename_dict[filename] = {"dt_seconds":dt_seconds, "dt":dt, "size":file_size, "path":file_path}
            
        #checking Md5 is not recommended as it needs to download all files first
    
    return filename_dict






local_paths = glob.glob(local_root_path, recursive = True)

remote_paths = [remote_root_path + os.path.basename(s) for s in local_paths]

h = "<!DOCTYPE html><html><head><title>VenSpec-H Document Tree Comparison</title></head><body>\n"
h += "<h1>VenSpec-H Document Management</h1>"

for local_path, remote_path in zip(local_paths, remote_paths):

    h += "<h2>VS/VEH/%s</h2>\n" %os.path.basename(local_path)
    h += "Comparing path <strong>%s</strong> on the BIRA internal server to path <strong>%s</strong> on NextCloud<br>\n" %(remote_path, local_path)

    d2c_dict = get_linux_tree_filenames(user, password, remote_path)
    nc_dict = get_windows_tree_filenames(local_path)
    
    
    def returnNotMatches(a, b):
        return [[x for x in a if x not in b], [x for x in b if x not in a]]
    
    
    not_in_nc, not_in_d2c = returnNotMatches(d2c_dict.keys(), nc_dict.keys())
    
    
    if len(not_in_nc) > 0 or len(not_in_d2c) > 0:
        h += "<h3>Missing files</h3>\n"
        
        if len(not_in_nc) > 0:
            h += "<h4>Missing from NextCloud</h4>\n"
            for filename in not_in_nc:
                h += filename + "<br>\n"
            h += "<br>\n"

        if len(not_in_d2c) > 0:
            h += "<h4>Missing from DOCS_TO_CLOUD</h4>\n"
            for filename in not_in_d2c:
                h += filename + "<br>\n"
            h += "<br>\n"

    else:
        h += "<h3>No missing files</h3>"
    


    
    table = ""
    for filename, info_d2c in d2c_dict.items():
        if filename in nc_dict.keys():
            info_nc = nc_dict[filename]
            dt_d2c = info_d2c["dt"]
            dt_nc = info_nc["dt"]
            
            dt_delta = (dt_d2c - dt_nc).total_seconds()
            
            if np.abs(dt_delta) > SYNC_LIMIT:
                
                size_d2c = info_d2c["size"]
                size_nc = info_nc["size"]
                
                table += "<tr><td>%s</td>" %filename
                if dt_d2c > dt_nc:
                    table += "<td><strong>%s</strong></td><td>%s</td>" %(datetime.strftime(dt_d2c, DT_FORMAT), datetime.strftime(dt_nc, DT_FORMAT))
                else:
                    table += "<td>%s</td><td><strong>%s</strong></td>" %(datetime.strftime(dt_d2c, DT_FORMAT), datetime.strftime(dt_nc, DT_FORMAT))
                
                if size_d2c > size_nc:
                    table += "<td><strong>%0.1fkb</strong></td><td>%0.1fkb</td>" %(size_d2c/1024.0, size_nc/1024.0)
                elif size_d2c < size_nc:
                    table += "<td>%0.1fkb</td><td><strong>%0.1fkb</strong></td>" %(size_d2c/1024.0, size_nc/1024.0)
                else:
                    table += "<td>%0.1fkb</td><td>%0.1fkb</td>" %(size_d2c/1024.0, size_nc/1024.0)
                    
                table += "</tr>"
                
    
    if table != "":
        h += "<h3>Out of sync files</h3>\n"
        h += "<table border=1><tr><th>Filename</th><th>DOCS_TO_CLOUD last modified (UTC)</th><th>NextCloud last modified (UTC)</th><th>DOCS_TO_CLOUD size</th><th>DOCS_TO_CLOUD size</th></tr>\n"
        h += table
        h += "</table><br><br>\n"
    else:
        h += "<h3>All files in sync</h3><br><br><br>"

h += "</body></html>\n"
with open("doc_comparison.html", "w") as f:
    f.write(h)