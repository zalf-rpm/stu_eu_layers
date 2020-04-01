#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

import time
import os
import math
import json
import csv
import copy
#from io import StringIO
from datetime import date, timedelta
from collections import defaultdict
import sys
import h5py
import numpy as np
import sqlite3
from pyproj import Proj, Transformer


def create_hdf5():

    top_files = [
        ["STU_EU_T_TEXT_CLS", "texture_class", np.uint8, np.uint8, None],
        ["STU_EU_T_TAWC", "tawc", float, np.uint16, lambda v: int(round(v, 2)*100)],
        ["STU_EU_T_SILT", "silt", np.uint8, np.uint8, None],
        ["STU_EU_T_SAND", "sand", np.uint8, np.uint8, None],
        ["STU_EU_T_CLAY", "clay", np.uint8, np.uint8, None],
        ["STU_EU_T_GRAVEL", "gravel", np.uint8, np.uint8, None],
        ["STU_EU_T_OC", "corg", float, np.uint16, lambda v: int(round(v, 2)*100)], 
        ["STU_EU_T_BD", "bulk_density", float, np.uint16, lambda v: int(v * 1000)],
        ["SMU_EU_T_TAWC", "tawc_smu", float, np.uint16, lambda v: int(round(v, 2)*100)]
    ]
    sub_files = [
        ["STU_EU_S_TEXT_CLS", "texture_class", np.uint8, np.uint8, None],
        ["STU_EU_S_TAWC", "tawc", float, np.uint16, lambda v: int(round(v, 2)*100)],
        ["STU_EU_S_SILT", "silt", np.uint8, np.uint8, None],
        ["STU_EU_S_SAND", "sand", np.uint8, np.uint8, None],
        ["STU_EU_S_CLAY", "clay", np.uint8, np.uint8, None],
        ["STU_EU_S_GRAVEL", "gravel", np.uint8, np.uint8, None],
        ["STU_EU_S_OC", "corg", float, np.uint16, lambda v: int(round(v, 2)*100)], 
        ["STU_EU_S_BD", "bulk_density", float, np.uint16, lambda v: int(v * 1000)],
        ["SMU_EU_S_TAWC", "tawc_smu", float, np.uint16, lambda v: int(round(v, 2)*100)]
    ]
    files = [
        ["STU_EU_DEPTH_ROOTS", "depth_roots", np.uint8, np.uint8, None],
        ["STU_EU_ALLOCATION", "stu_allocation", np.bool, np.bool,None]
    ]

    #all_files = []
    #for fs in [files, top_files, sub_files]:
    #    all_files.extend(fs)

    path_to_data = "C:/Users/berg.ZALF-AD/Desktop/anna-eu-bodendaten/STU_EU_Layers/"

    with h5py.File("stu_eu_layers.hdf5", "a") as hdf:
        
        #general
        general = hdf["/general"] if "general" in hdf else hdf.create_group("general")
        for f in files:
            print("general:", f)
            ds_name = f[1]
            grid = np.loadtxt(path_to_data + f[0] + ".asc", dtype=f[2], skiprows=5)
            ds = general[ds_name] if ds_name in general else general.create_dataset(ds_name, grid.shape, dtype=f[3])
            ds[:] = grid[:]

        top = hdf["/top"] if "top" in hdf else hdf.create_group("top")
        for f in top_files:
            print("top:", f)
            ds_name = f[1]
            npr = np.frompyfunc(f[4], 1, 1) if f[4] else None
            grid = np.loadtxt(path_to_data + f[0] + ".asc", dtype=f[2], skiprows=5)
            grid2 = np.empty_like(grid, f[3])
            grid2[:] = npr(grid)[:] if npr else grid
            ds = top[ds_name] if ds_name in top else top.create_dataset(ds_name, grid.shape, dtype=f[3])
            ds[:] = grid2[:]

        sub = hdf["/sub"] if "sub" in hdf else hdf.create_group("sub")
        for f in sub_files:
            print("sub:", f)
            ds_name = f[1]
            npr = np.frompyfunc(f[4], 1, 1) if f[4] else None
            grid = np.loadtxt(path_to_data + f[0] + ".asc", dtype=f[2], skiprows=5)
            grid2 = np.empty_like(grid, f[3])
            grid2[:] = npr(grid)[:] if npr else grid
            ds = sub[ds_name] if ds_name in sub else sub.create_dataset(ds_name, grid.shape, dtype=f[3])
            ds[:] = grid2[:]


def create_csv():

    wgs84 = Proj(init="epsg:4326") #proj4 -> (World Geodetic System 1984 https://epsg.io/4326)
    etrs89 = Proj(init="EPSG:3035") 
    transformer = Transformer.from_proj(etrs89, wgs84) 

    #cols = 5900
    rows = 4600
    cellsize = 1000
    xll = 1500000 
    yll = 900000
    #nodata_value = 0
    xll_center = xll + cellsize // 2
    yll_center = yll + cellsize // 2
    yul_center = yll_center + (rows - 1)*cellsize

    with h5py.File("stu_eu_layers.hdf5", "r") as hdf:
        alloc = hdf["/general/stu_allocation"]
        shape = alloc.shape
        
        path_to_csv_file = "stu_eu_layers.csv"
        if not os.path.isfile(path_to_csv_file):
            with open(path_to_csv_file, "w") as _:
                _.write("col,row,elevation,latitude,longitude,depth,OC_topsoil,OC_subsoil,BD_topsoil,BD_subsoil,Sand_topsoil,Clay_topsoil,Silt_topsoil,Sand_subsoil,Clay_subsoil,Silt_subsoil\n")

        with open(path_to_csv_file, "a", newline="") as _:
            writer = csv.writer(_, delimiter=",")

            depth_roots = hdf["/general/depth_roots"]

            corg_t = hdf["/top/corg"]
            corg_s = hdf["/sub/corg"]

            bd_t = hdf["/top/bulk_density"]
            bd_s = hdf["/sub/bulk_density"]

            sand_t = hdf["/top/sand"]
            sand_s = hdf["/sub/sand"]

            clay_t = hdf["/top/clay"]
            clay_s = hdf["/sub/clay"]

            silt_t = hdf["/top/silt"]
            silt_s = hdf["/sub/silt"]

            for row in range(83, shape[0]):
                print(row, sep=" ")
                for col in range(shape[1]-500):
                    # create a row
                    if alloc[row, col] == 1:
                        r = xll_center + col * cellsize
                        h = yul_center - row * cellsize
                        lon, lat = transformer.transform(r, h)

                        line = [
                            col,
                            row,
                            0, #elevation
                            round(lat, 4), 
                            round(lon, 4),
                            round(depth_roots[row, col] / 100, 2),
                            round(corg_t[row, col] / 100, 4),
                            round(corg_s[row, col] / 100, 4),
                            bd_t[row, col],
                            bd_s[row, col],
                            sand_t[row, col],
                            clay_t[row, col],
                            silt_t[row, col],
                            sand_s[row, col],
                            clay_s[row, col],
                            silt_s[row, col]
                        ]
                        writer.writerow(line)    
            



if __name__ == "__main__":
    #create_hdf5()
    create_csv()
