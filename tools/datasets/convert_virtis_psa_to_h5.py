# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 08:58:35 2023

@author: iant

READ VIRTIS CALIBRATED PSA QUBE AND CONVERT TO NOMAD-LIKE H5 FILE
CURRENTLY WORKS FOR VIRTIS M IR ONLY I.E. FILES STARTING WITH VI


The RECORD_BYTES data element identifies the number of bytes in each physical
record of the data product file. Records length is always equal to 512 bytes.= 512
The FILE_RECORDS data element identifies the number of physical records in the file. = 222574
The LABEL_RECORDS data element identifies the number of physical records that make up the PDS
product label. = 12




A QUBE core is a 3-dimension structure containing science measurements
The axes of the QUBE are called
BAND (spectral band defined by the wavelength) 
SAMPLE (spatial direction along the slit)
LINE (acquisitions in successive time steps)

wavelength x row x frame





single precision floating point, IEEE encoding, BIP interleave mode
BIP interleave is best for spectral processing, as the values for all the bands of a given pixel are adjacent

For calibrated data QUBES, the backplane is used to store the original acquisition time
(SCET). Individual parameters are stored as 2-bytes integers, but always appear in even number

The M calibrated files contain one qube with a backplane and a bottomplane.
All files first contain a HISTORY object, 512 bytes long, filled with ASCII 0. This is a relic from
previous requirements that are no longer relevant.

The QUBE contains the data calibrated in radiance. Dark frames are removed from the qube
The x dimension of the qube provides spectral measurements in increasing wavelength order. Data are scaled in W/m2/sr/ìm and are
stored as binary floating points on 4 bytes, MSB encoding. A single row of the backplane is used to store the
SCET associated to the whole frame, encoded on the first 3 4-bytes integers of this row (with only two bytes
used out of 4) (see section 3.2.2). This SCET is reconstructed at mid-exposure.
The bottomplane has 3 frames containing reference data for each spectel: wavelength and FWHM
(in microns), then absolute uncertainty on the signal (1-sigma deviation). There is no sideplane associated to
this QUBE.

 AXIS_NAME                   = (BAND,SAMPLE,LINE)
 CORE_ITEMS                  = (432,64,1025)

 SUFFIX_BYTES                = 4
 SUFFIX_ITEMS                = (1,0,3)

 BAND_SUFFIX_NAME            = "SCET"
 BAND_SUFFIX_UNIT            = DIMENSIONLESS
 BAND_SUFFIX_ITEM_BYTES      = 2

 LINE_SUFFIX_NAME            = ("WAVELENGTH","FWHM","UNCERTAINTY")
 LINE_SUFFIX_UNIT            = ("MICRON","MICRON","W/m**2/sr/micron")
 LINE_SUFFIX_ITEM_BYTES      = 4

"""
import re
import os
import numpy as np
# import matplotlib.pyplot as plt
import struct
import glob
import platform


from tools.file.read_write_hdf5 import write_hdf5_from_dict

if platform.system() == "Windows":
    psa_data_root_path = r"C:\Users\iant\Documents\DATA\venus\virtis\VEX-V-VIRTIS-2-3-V3.0\DATA"
    h5_data_root_path = r"C:\Users\iant\Documents\DATA\venus\virtis"
elif platform.system() == "Linux":
    psa_data_root_path = r"/bira-iasb/projects/work/EnVision/Data/VIRTIS/psa/VEX-V-VIRTIS-2-3-V3.0/DATA"
    h5_data_root_path = r"/bira-iasb/projects/work/EnVision/Data/VIRTIS/hdf5/hdf5_level_1p0a/"


# H: H image transfer mode (backup observation mode)
# S: H single spectrum transfer mode (including dark current files in nominal mode)
# T: H "64-spectra frame" transfer mode (nominal mode)
# I: M-IR data
# V: M-Vis data

data_filepaths = glob.glob(f'{psa_data_root_path}{os.sep}**{os.sep}VI*.CAL', recursive=True)
# data_filepaths = glob.glob(f'{psa_data_root_path}{os.sep}**{os.sep}VT*.CAL', recursive=True)





CORE_NULL = -1004
CORE_HIGH_REPR_SATURATION = -1001


def get_label(regex, header, dtype="int"):
    """run regex on the header text to find a value or multiple values on one line. Returns a list"""
    
    result = [re.findall(regex, line) for line in header if re.findall(regex, line)]
    if len(result) == 1:
        if isinstance(result[0][0], tuple):
            #if tuple returned
            if dtype == "int":
                return [int(i) for i in list(result[0][0])]
        
        else:
            #if integer return
            if dtype == "int":
                return [int(result[0][0])]
    else:
        print("Regex error: %i results found" %len(result))
        return []
    



def unpack(bytes_, start, n_vals, c_type="f4"):
    """unpack bytes to 4 byte floats or 4 byte integers"""
    
    if c_type == "f4":
        n_bytes = 4
        letter = "f"
    if c_type == "i4":
        n_bytes = 4
        letter = "i"
    
    if n_vals == -1:
        n_vals = int((len(bytes_) - start) / n_bytes)
    
    end = start + n_vals * n_bytes
    
    out = struct.unpack('>'+letter*n_vals, bytes_[start:end])
    if len(out) == 1:
        out = out[0]
    return out, end



#map M-IR data to NOMAD-like geometry fields
#second value is the scaler to convert stored int to real value
mappings = [
    ["Geometry/Point1/Lon", 10000.],
    ["Geometry/Point2/Lon", 10000.],
    ["Geometry/Point3/Lon", 10000.],
    ["Geometry/Point4/Lon", 10000.],

    ["Geometry/Point1/Lat", 10000.],
    ["Geometry/Point2/Lat", 10000.],
    ["Geometry/Point3/Lat", 10000.],
    ["Geometry/Point4/Lat", 10000.],

    ["Geometry/Point0/Lon", 10000.],
    ["Geometry/Point0/Lat", 10000.],

    ["Geometry/Point0/SZA", 10000.],
    ["Geometry/Point0/EmissionAngle", 10000.],
    ["Geometry/Point0/PhaseAngle", 10000.],

    ["Geometry/Point0/SurfaceAltEllipsoid", 1000.],
    ["Geometry/Point0/SlantDistance", 1000.],
    ["Geometry/Point0/LST", 100000./24.],

    ["Geometry/Point1/Lon60km", 10000.],
    ["Geometry/Point2/Lon60km", 10000.],
    ["Geometry/Point3/Lon60km", 10000.],
    ["Geometry/Point4/Lon60km", 10000.],

    ["Geometry/Point1/Lat60km", 10000.],
    ["Geometry/Point2/Lat60km", 10000.],
    ["Geometry/Point3/Lat60km", 10000.],
    ["Geometry/Point4/Lat60km", 10000.],

    ["Geometry/Point0/Lon60km", 10000.],
    ["Geometry/Point0/Lat60km", 10000.],

    ["Geometry/Point0/SurfaceAltEllipsoid60km", 1000.],
    ]




#loop through files
for data_filepath in data_filepaths:
    
    print(data_filepath)
    
    geom_filepath = data_filepath.replace("CALIBRATED", "GEOMETRY").replace(".CAL", ".GEO")


    data_in = open(data_filepath, "rb").read()
    
    
    #read in first part of file, decode to strings
    header = data_in[:1000].decode("utf-8").splitlines()
        
    
    
    
    #find file size data
    record_bytes = get_label("RECORD_BYTES\s+=\s+(\d+)", header)[0]
    label_records = get_label("LABEL_RECORDS\s+=\s+(\d+)", header)[0]
    qube_record = get_label("\^QUBE\s+=\s+(\d+)", header)[0]
    
    #get full header now we know the size
    header = data_in[:record_bytes*label_records].decode("utf-8") .splitlines()
    
    
    #now search header for start time, qube size and aux data size
    dt_list = get_label("START_TIME\s+=\s+(\d+)-(\d+)-(\d+)T(\d+):(\d+):(\d+).\d+", header)
    
    qube_data_size = get_label("\sCORE_ITEMS\s+=\s+\((\d+),(\d+),(\d+)\)", header)
    aux_data_size = get_label("\sSUFFIX_ITEMS\s+=\s+\((\d+),(\d+),(\d+)\)", header)
    
    #real dimension size is data qube + suffix
    dim_size = [qube_data_size[0]+aux_data_size[0], qube_data_size[1], qube_data_size[2]]
    
    #get bytes from file from start byte number to start byte + size, then unpack
    start_byte = record_bytes*(qube_record-1)
    qube_data = data_in[start_byte:(start_byte + np.prod(dim_size)*4)]
    data_qube, b = unpack(qube_data, 0, np.prod(dim_size))
    
    #convert to 3d for checking/plotting and replace invalid values with nans
    data_qube_arr = np.asfarray(data_qube).reshape((-1, qube_data_size[1], qube_data_size[0]+aux_data_size[0]))
    data_qube_arr[data_qube_arr < -999.0] = np.nan
    
    #convert to 2d for h5 output
    data_qube_arr_2d = np.reshape(data_qube_arr, (-1, dim_size[0]))
    
    #remove suffix column
    data_qube_arr_2d = np.delete(data_qube_arr_2d, -1, axis=1)
    
    
    # qube_binned = np.nanmean(qube_arr, axis=1)
    # for i in range(qube_binned.shape[0]):
    #     plt.plot(qube_binned[i, :])

    #repeat for the 3 suffixes: wavelength, fwhm, and yerror    
    start_byte = start_byte + np.prod(dim_size)*4
    um_data = data_in[start_byte:(start_byte + qube_data_size[0]*qube_data_size[1]*4)]
    um_qube, b = unpack(um_data, 0, -1)
    um_qube_arr = np.asfarray(um_qube).reshape((qube_data_size[1], qube_data_size[0]))
    
    # plt.imshow(um_qube_arr)
    
    #get bytes from file from start byte number to start byte + size, then unpack
    start_byte = start_byte + qube_data_size[0]*qube_data_size[1]*4
    fwhm_data = data_in[start_byte:(start_byte + qube_data_size[0]*qube_data_size[1]*4)]
    fwhm_qube, b = unpack(fwhm_data, 0, -1)
    fwhm_qube_arr = np.asfarray(fwhm_qube).reshape((qube_data_size[1], qube_data_size[0]))
    
    # plt.figure()
    # plt.imshow(fwhm_qube_arr)
    
    
    #get bytes from file from start byte number to start byte + size, then unpack
    start_byte = start_byte + qube_data_size[0]*qube_data_size[1]*4
    yerror_data = data_in[start_byte:(start_byte + qube_data_size[0]*qube_data_size[1]*4)]
    yerror_qube, b = unpack(yerror_data, 0, -1)
    yerror_qube_arr = np.asfarray(yerror_qube).reshape((qube_data_size[1], qube_data_size[0]))
    
    #ignore padding zeros at end of file
    
    
    
    #now repeat for the geometry
    
    geom_in = open(geom_filepath, "rb").read()
    
    #read in first part of file, decode to strings
    header = geom_in[:1000].decode("utf-8").splitlines()
    
    #find file size data
    record_bytes = get_label("RECORD_BYTES\s+=\s+(\d+)", header)[0]
    label_records = get_label("LABEL_RECORDS\s+=\s+(\d+)", header)[0]
    qube_record = get_label("\^QUBE\s+=\s+(\d+)", header)[0]
    
    #get full header now we know the size
    header = geom_in[:record_bytes*label_records].decode("utf-8") .splitlines()
    
    
    #now search header for geom qube size and aux data size
    qube_data_size = get_label("\sCORE_ITEMS\s+=\s+\((\d+),(\d+),(\d+)\)", header)
    aux_data_size = get_label("\sSUFFIX_ITEMS\s+=\s+\((\d+),(\d+),(\d+)\)", header)
    
    dim_size = [qube_data_size[0]+aux_data_size[0], qube_data_size[1], qube_data_size[2]]
    
    #get bytes from file from start byte number to start byte + size, then unpack
    start_byte = record_bytes*(qube_record-1)
    qube_data = geom_in[start_byte:(start_byte + np.prod(dim_size)*4)]
    qube, b = unpack(qube_data, 0, np.prod(dim_size), c_type="i4")
    
    #convert to 3d for checking/plotting
    qube_arr = np.asfarray(qube).reshape((-1, qube_data_size[1], qube_data_size[0]+aux_data_size[0]))
    
    #convert to 2d for h5 output
    qube_arr_2d = np.reshape(qube_arr, (-1, qube_data_size[0]))
    
    
    
    
    #make dictionary for writing to h5 output file
    h5_d = {}
    for i, (h5_field, scaling_factor) in enumerate(mappings):
        h5_d[h5_field] = qube_arr_2d[:, i] / scaling_factor
    
    h5_d["Science/X"] = um_qube_arr[0, :]
    h5_d["Science/Y"] = data_qube_arr_2d
    h5_d["Channel/FWHM"] = fwhm_qube_arr[0, 0]
    
    year = "%04i" %dt_list[0]
    month = "%02i" %dt_list[1]
    day = "%02i" %dt_list[2]
    
    h5_filename = "%04i%02i%02i_%02i%02i%02i_1p0a_VIRTISMIR_D" %(*dt_list,)
    
    #make year/month/day directory structure for output file
    year_path = os.path.join(h5_data_root_path, year)
    os.makedirs(year_path, exist_ok=True)
    month_path = os.path.join(h5_data_root_path, year, month)
    os.makedirs(month_path, exist_ok=True)
    day_path = os.path.join(h5_data_root_path, year, month, day)
    os.makedirs(day_path, exist_ok=True)
    
    h5_filepath = os.path.join(day_path, h5_filename)
    
    #write to h5
    write_hdf5_from_dict(h5_filepath, h5_d, {}, {}, {})


# code for dumping hex data
# hex_out = geom_in[:10000].hex()
# with open("text.txt", "w") as f:
#     # f.write(hex_out)
#     i = 0
    
#     for j in range(100):
#         step = 33
#         print(i*8, (i+step)*8)
#         f.write(hex_out[i*8:(i+step)*8]+"\n")
#         i += step
        
        # step = 1
        # print(i*8, (i+step)*8)
        # f.write(hex_out[i*8:(i+step)*8]+"\n")
        # i += step


    


    
