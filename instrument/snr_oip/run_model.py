# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 10:42:47 2023

@author: iant

RUN OIP SNR MODEL EXE AND SAVE OUTPUT TO DIFFERENT FOLDERS

ALL CONF INPUT FILES NEED TO BE PRESENT TO RUN SUCCESSFULLY
"""
import numpy as np
import os
import shutil
import json
import subprocess
import time
from datetime import datetime
from instrument.snr_oip.config import BAND_DATA, MODEL_RESULTS, MODEL_EXE, MODEL_ROOT, MODEL_OUTPUT, GAIN_DATA

bands = BAND_DATA.keys()

temps_detector = np.arange(111., 140., 1.) #K
# temps_detector = []
# temps_detector = []

temps_warmsec = np.arange(263.15, 268.16, 5.0)
temps_coldsec = np.arange(213.15, 243.16, 5.0) #K
gains = ["low", "high"]

#make results dir if doesn't yet exist
if not os.path.exists(MODEL_RESULTS):
    print("Making directory %s" %MODEL_RESULTS)
    os.mkdir(MODEL_RESULTS)


for temp_detector in temps_detector:
    for temp_warmsec in temps_warmsec:
        for temp_coldsec in temps_coldsec:

    
            #remove existing conf files
            for band in BAND_DATA.keys():
                if os.path.exists(BAND_DATA[band]["conf"]):
                    # print("Deleting %s config file" %band)
                    os.remove(BAND_DATA[band]["conf"])
            
            for band in bands:
                
                #read in config file
                conf_path = BAND_DATA[band]["conf_temp"]
                with open(conf_path, "r") as f:
                    config = json.load(f)
            
                #make changes to conf file
                config["TempOffset"] = 0.0 #don't vary by temp
                config["TempDet"] = temp_detector #don't vary by temp
                config["TempIll"] = temp_warmsec #don't vary by temp
                config["TempSpec"] = temp_coldsec #don't vary by temp
                
                
                    
                #write config file to dir to run snr model
                # print("Making band %s config file" %band)
                conf_run = BAND_DATA[band]["conf"]
                with open(conf_run, "w") as f:
                    json_str = json.dumps(config, indent=4)
                    f.write(json_str)
                
            
            
            #delete output directory and files
            if os.path.exists(MODEL_OUTPUT):
                print("Deleting %s" %MODEL_OUTPUT)
                shutil.rmtree(MODEL_OUTPUT)
            
            
            #run model, remaking output directory
            print("Running model Td=%0.1f Tw=%0.1f_Tc=%0.1f" %(temp_detector, temp_warmsec, temp_coldsec))
            subprocess.run([MODEL_EXE])
            
            output = subprocess.Popen(MODEL_EXE, cwd=MODEL_ROOT, stdout=subprocess.PIPE).communicate()[0].decode()
                
            time.sleep(1)
            #rename output directory to reflect temperatures
            new_name = "out_d%0.1f_w%0.1f_c%0.1f" %(temp_detector, temp_warmsec, temp_coldsec)
            new_dir = os.path.join(MODEL_RESULTS, new_name)
            
            print("Renaming output directory to %s (%s)" %(new_name, str(datetime.now())[:-7]))
            os.rename(MODEL_OUTPUT, new_dir)