#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR: Shen Ruoque
The verification point file is used to obtain the type of the corresponding point on the identification diagram and use it as verification.
利用验证点文件，求取对应点在识别图上的类型，用作验证。
"""
#%%
import os, re
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import json
#%%
rx = re.compile(r'"coordinates": *\[ *(\d+\.\d+), (\d+\.\d+), ')
placeyrs = [
    "SouthKorea-2018",
    "SouthKorea-2021"
]
version = "0"
os.chdir("./path")
path1 = "the/way/to/the/rice/maps" #识别图路径
path2 = "the/way/to/the/validation_points" #验证点路径

for placeyr in placeyrs:
    print(placeyr)
    yr = placeyr[-4:]
    print(yr)
    with rasterio.open(os.path.join(
        path1, f"classified-.tif" #识别图数据
    )) as ds:
        classified = ds.read()# == 0
        width = ds.width
        height = ds.height
        transform = ds.transform
        crs = ds.crs

    points = gpd.read_file(f"{path2}/{placeyr}-30m.geojson") #验证点文件
    points.rename({"Name": "class"})
    xs = (points.geometry.x.to_numpy() - transform[2]) / transform[0]
    ys = (points.geometry.y.to_numpy() - transform[5]) / transform[4]
    data = np.zeros_like(xs, dtype=np.uint8)
    for i in range(len(xs)):
        y, x = int(ys[i]), int(xs[i])
        if (0 <= y < height) and (0 <= x < width): #x,y为验证点坐标
            data[i] = classified[0, y, x]
        else:
            data[i] = 0
    points["b1"] = data

    del classified

    points[["class", "b1", "geometry"]].to_csv(
        f"{placeyr}-point_validation-{version}.csv", #输出文件
        index=False, encoding="utf-8"
    )
    print("完成")

#%%