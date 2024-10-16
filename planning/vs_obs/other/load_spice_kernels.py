# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:30:37 2022

@author: iant
"""
import os
# from datetime import datetime
import spiceypy as sp


from planning.vs_obs.config.paths import paths


def load_spice_kernels():
    kernel_root_dir = paths["KERNEL_ROOT_DIRECTORY"]

    # load each kernel manually

    kernels = [
        'fk/envision_v02.tf',
        'fk/estrack_v04.tf',
        'fk/earth_topo_050714.tf',
        'fk/earthfixeditrf93.tf',

        'ik/envision_venspec_v02.ti',

        'lsk/naif0012.tls',

        'sclk/envision_200715_fict.tsc',

        'pck/pck00010.tpc',
        'pck/de-403-masses.tpc',
        'pck/earth_200101_990628_predict.bpc',

        'spk/de432s.bsp',
        'spk/estrack_v04.bsp',
        'spk/earthstns_itrf93_050714.bsp',
        'spk/EnVision_ET1_2031_NorthVOI_ML008_v01.bsp',

        'ck/EnVision_ET1_2031_NorthVOI_ML008_v01.bc',
    ]

    kernel_name = "envision_study_et1_north_voi_ml008_v051.tm"
    print("KERNEL_DIRECTORY=%s, METAKERNEL_NAME=%s" % (kernel_root_dir, kernel_name))

    for kernel in kernels:
        sp.furnsh(os.path.join(kernel_root_dir, os.path.normcase(kernel)))

    # sp.furnsh(os.path.join(kernel_root_dir, "mk", kernel_name))
