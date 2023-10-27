# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 10:27:34 2023

@author: iant

RUN OIP SNR MODEL
"""

import os

MODEL_ROOT = r"C:\Users\iant\Documents\PROGRAMS\venspech_snr_model"

MODEL_EXE = os.path.join(MODEL_ROOT, "VENSPEC-H SNR v1.00.exe")
MODEL_OUTPUT = os.path.join(MODEL_ROOT, "output")
MODEL_DATA = os.path.join(MODEL_ROOT, "data")
MODEL_TEMP = os.path.join(MODEL_ROOT, "temp")

MODEL_RESULTS = os.path.join(MODEL_ROOT, "results")

GAIN_DATA = {
    "low":{
        "Capacity": 875000.0, "RON": 500,
        },
    "high":{
        "Capacity": 296000.0, "RON": 150,
        }}

BAND_DATA = {
    "1_NIGHT":{
        "measured":True,
        "rows":200,
        "ref_wavel":1.170,
        "conf":os.path.join(MODEL_ROOT, "VEH1_NIGHT.conf"),
        "conf_temp":os.path.join(MODEL_TEMP, "VEH1_NIGHT.conf"),
        "data":os.path.join(MODEL_DATA, "TSPEC_VEH1.dat"),
        "spectrum":os.path.join(MODEL_DATA, "ENVIS-VS-VEH-SN0016-iss0rev0+AER-NightSpectrum_20221116.dat"),
        "snr":"VEH1-NIGHT.txt",
        "snr_vs_t":os.path.join(MODEL_OUTPUT, "VEH1-NIGHT_VarT.txt"),
    },

    # "1_DAY":{
    #     "measured":False,
    #     "conf":os.path.join(MODEL_ROOT, "VEH1_DAY.conf"),
    #     "conf_temp":os.path.join(MODEL_TEMP, "VEH1_DAY.conf"),
    #     "data":os.path.join(MODEL_DATA, "TSPEC_VEH1.dat"),
    #     "spectrum":os.path.join(MODEL_DATA, "ENVIS-VS-VEH-SN-0017-iss0rev0+AER-DaySpectrum_20221116.dat"),
    #     "snr":"VEH1-DAY.txt",
    #     "snr_vs_t":os.path.join(MODEL_OUTPUT, "VEH1-DAY_VarT.txt"),
    # },



    # "2A_DAY":{
    #     "measured":True,
    #     "rows":18,
    #     "ref_wavel":2.380,
    #     "conf":os.path.join(MODEL_ROOT, "VEH2A_DAY.conf"),
    #     "conf_temp":os.path.join(MODEL_TEMP, "VEH2A_DAY.conf"),
    #     "data":os.path.join(MODEL_DATA, "TSPEC_VEH2A.dat"),
    #     "spectrum":os.path.join(MODEL_DATA, "ENVIS-VS-VEH-SN-0017-iss0rev0+AER-DaySpectrum_20221116.dat"),
    #     "snr":"VEH2A-DAY.txt",
    #     "snr_vs_t":os.path.join(MODEL_OUTPUT, "VEH2A-DAY_VarT.txt"),
    # },
    # "2A_NIGHT":{
    #     "measured":True,
    #     "rows":18,
    #     "ref_wavel":2.380,
    #     "conf":os.path.join(MODEL_ROOT, "VEH2A_NIGHT.conf"),
    #     "conf_temp":os.path.join(MODEL_TEMP, "VEH2A_NIGHT.conf"),
    #     "data":os.path.join(MODEL_DATA, "TSPEC_VEH2A.dat"),
    #     "spectrum":os.path.join(MODEL_DATA, "ENVIS-VS-VEH-SN0016-iss0rev0+AER-NightSpectrum_20221116.dat"),
    #     "snr":"VEH2A-NIGHT.txt",
    #     "snr_vs_t":os.path.join(MODEL_OUTPUT, "VEH2A-NIGHT_VarT.txt"),
    # },



    # "2B_DAY":{
    #     "measured":True,
    #     "rows":182,
    #     "ref_wavel":2.460,
    #     "conf":os.path.join(MODEL_ROOT, "VEH2B_DAY.conf"),
    #     "conf_temp":os.path.join(MODEL_TEMP, "VEH2B_DAY.conf"),
    #     "data":os.path.join(MODEL_DATA, "TSPEC_VEH2B.dat"),
    #     "spectrum":os.path.join(MODEL_DATA, "ENVIS-VS-VEH-SN-0017-iss0rev0+AER-DaySpectrum_20221116.dat"),
    #     "snr":"VEH2B-DAY.txt",
    #     "snr_vs_t":os.path.join(MODEL_OUTPUT, "VEH2B-DAY_VarT.txt"),
    # },
    # "2B_NIGHT":{
    #     "measured":True,
    #     "rows":182,
    #     "ref_wavel":2.460,
    #     "conf":os.path.join(MODEL_ROOT, "VEH2B_NIGHT.conf"),
    #     "conf_temp":os.path.join(MODEL_TEMP, "VEH2B_NIGHT.conf"),
    #     "data":os.path.join(MODEL_DATA, "TSPEC_VEH2B.dat"),
    #     "spectrum":os.path.join(MODEL_DATA, "ENVIS-VS-VEH-SN0016-iss0rev0+AER-NightSpectrum_20221116.dat"),
    #     "snr":"VEH2B-NIGHT.txt",
    #     "snr_vs_t":os.path.join(MODEL_OUTPUT, "VEH2B-NIGHT_VarT.txt"),
    # },



    # "3_DAY":{
    #     "measured":False,
    #     "conf":os.path.join(MODEL_ROOT, "VEH3_DAY.conf"),
    #     "conf_temp":os.path.join(MODEL_TEMP, "VEH3_DAY.conf"),
    #     "data":os.path.join(MODEL_DATA, "TSPEC_VEH3.dat"),
    #     "spectrum":os.path.join(MODEL_DATA, "ENVIS-VS-VEH-SN-0017-iss0rev0+AER-DaySpectrum_20221116.dat"),
    #     "snr":"VEH3-DAY.txt",
    #     "snr_vs_t":os.path.join(MODEL_OUTPUT, "VEH3-DAY_VarT.txt"),
    # },
    # "3_NIGHT":{
    #     "measured":True,
    #     "rows":200,
    #     "ref_wavel":1.730,
    #     "conf":os.path.join(MODEL_ROOT, "VEH3_NIGHT.conf"),
    #     "conf_temp":os.path.join(MODEL_TEMP, "VEH3_NIGHT.conf"),
    #     "data":os.path.join(MODEL_DATA, "TSPEC_VEH3.dat"),
    #     "spectrum":os.path.join(MODEL_DATA, "ENVIS-VS-VEH-SN0016-iss0rev0+AER-NightSpectrum_20221116.dat"),
    #     "snr":"VEH3-NIGHT.txt",
    #     "snr_vs_t":os.path.join(MODEL_OUTPUT, "VEH3-NIGHT_VarT.txt"),
    # },



    # "4_DAY":{
    #     "measured":True,
    #     "rows":200,
    #     "ref_wavel":1.380,
    #     "conf":os.path.join(MODEL_ROOT, "VEH4_DAY.conf"),
    #     "conf_temp":os.path.join(MODEL_TEMP, "VEH4_DAY.conf"),
    #     "data":os.path.join(MODEL_DATA, "TSPEC_VEH4.dat"),
    #     "spectrum":os.path.join(MODEL_DATA, "ENVIS-VS-VEH-SN-0017-iss0rev0+AER-DaySpectrum_20221116.dat"),
    #     "snr":"VEH4-DAY.txt",
    #     "snr_vs_t":os.path.join(MODEL_OUTPUT, "VEH4-DAY_VarT.txt"),
    # },
    # "4_NIGHT":{
    #     "measured":False,
    #     "conf":os.path.join(MODEL_ROOT, "VEH4_NIGHT.conf"),
    #     "conf_temp":os.path.join(MODEL_TEMP, "VEH4_NIGHT.conf"),
    #     "data":os.path.join(MODEL_DATA, "TSPEC_VEH4.dat"),
    #     "spectrum":os.path.join(MODEL_DATA, "ENVIS-VS-VEH-SN0016-iss0rev0+AER-NightSpectrum_20221116.dat"),
    #     "snr":"VEH4-NIGHT.txt",
    #     "snr_vs_t":os.path.join(MODEL_OUTPUT, "VEH4-NIGHT_VarT.txt"),
    # },
    
}