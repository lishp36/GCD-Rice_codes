#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR: Shen Ruoque
VERSION: v2023.11.23

Count the distribution of good observations in each province.
"""
#%%
import rasterio
import os
import numpy as np
import pandas as pd

#%%
def gtiffref(file):
    with rasterio.open(file) as ds:
        ref = {
            "width": ds.width,
            "height": ds.height,
            "count": ds.count,
            "dtype": ds.dtypes[0],
            "transform": ds.transform,
            "crs": ds.crs,
        }
    return ref

def readgtiff(file):
    with rasterio.open(file) as ds:
        data = ds.read()
    return data, gtiffref(file)

def writegtiff(file, data, ref, compress="DEFLATE", bigtiff=False, nthread=None):
    ref["count"] = data.shape[0]
    ref["dtype"] = data.dtype
    if compress is not None: ref["compress"] = compress; ref["tiled"] = True
    if bigtiff: ref["BIGTIFF"] = "YES"
    if nthread is not None: ref["NUM_THREADS"] = nthread
    with rasterio.open(file, "w", driver="GTiff", **ref) as ds:
        ds.write(data)
    return 0

#%%
places = [
    "SouthKorea"
]


def getdays(place):
    if place in ["SouthKorea"]:
        return "137_8_225_day", range(137, 225, 8)
    elif place in [""]:
        return "", range(, , )
    else:
        return ""

#%%
path = ""
yrs = range(1990, 2023+1)

for place in places:
    days, ranges = getdays(place)
    pixnums = np.zeros(len(ranges)+1)#csv表的列数，返回数组
    hist = pd.DataFrame()
    for yr in yrs:
        print(place, yr)
        file = f"{path}/{place}/{place}-{yr}-L57_num-{days}-WGS84.tif"
        pixnums = np.histogram(readgtiff(file)[0][0, :, :], bins=range(0, len(ranges)+2))[0]
        hist[f"{yr}_pixel_frac"] = pixnums / np.sum(pixnums)
    hist.to_csv(f"{path}/{place}-all-observation-frequency-cropland.csv", index=False)
#%%