#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR: Shen Ruoque
VERSION: v2023.08.07

Mask non-cropland pixels in resent rice map.
(Set the values of non-cropland pixels to 2).
"""
#%%
import os
import multiprocessing as mp
import rasterio
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

def writegtiff(file, data, ref, compress="DEFLATE", bigtiff=False, nodata=None, nthread=None):
    ref["count"] = data.shape[0]
    ref["dtype"] = data.dtype
    if compress is not None: ref["compress"] = compress; ref["tiled"] = True
    if bigtiff: ref["BIGTIFF"] = "YES"
    if nodata is not None: ref["nodata"] = nodata
    if nthread is not None: ref["NUM_THREADS"] = nthread
    with rasterio.open(file, "w", driver="GTiff", **ref) as ds:
        ds.write(data)
    return file

#%%
res30 = 0.000269494585235856472
res20 = 0.000179663056823904
res10 = 0.000089831528412

res = res30
vecpath = "the/way/to/the/boundary_shapefile.shp"
places = [
    "SouthKorea"
]

maskpath = f"the/way/to/the/cropland"

for place in places:
    print(place)
    inpath = f"the/way/to/the/cropland/twdtw-classified-Result-maps"
    outpath = f"the/masked/path"
    maskfile = f"{maskpath}/Cropland.tif"

    def call_mask(yr):
        print(yr)
        data, ref = readgtiff(f"{inpath}/twdtw-classified-Result-maps.tif")
        mask, _ = readgtiff(f"{maskpath}/Cropland.tif")

        data[mask != 1] = 2
        writegtiff(
            f"{outpath}/-crop-cropmask.tif",
            data, ref
        )
        del mask, data

    yrs = range(2016,2023+1)
    pool = mp.Pool(len(yrs))
    pool.map(call_mask, yrs)
    pool.close()
    pool.join()
#%%