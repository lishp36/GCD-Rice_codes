#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR: Shen Ruoque
VERSION: v2024.07.19
"""
#%%
import numpy as np
import os, re, yaml
import rasterio
import multiprocessing as mp

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

def writegtiff(file, data, ref0, compress="DEFLATE", bigtiff=False, nodata=None, nthread=None):
    ref = ref0.copy()
    if data.ndim == 2: data = data[np.newaxis, :, :]
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
place = "SouthKorea"
version = "rf7"
outpath = "path/classified-Result/"
yrs = range(1990, 2023+1)
area_file = "path/area.yaml"

#%% fill zero observation
status = os.system(
    "julia -t 20 path/fill-zero-observation.jl " +
    f"{place} {version}"
)
 
#%% recognize
def recognize(yr):
    print(yr)
    outfile = f"{outpath}/classified-{place}-{yr}-rice-v{version}f-2.tif"
    if os.path.isfile(outfile): return 0
    probfile = f"{outpath}/prob-{place}-{yr}-rice-v{version}f-2.tif"
    area = yaml.safe_load(open(area_file))[place][yr]
    status = os.system(
        "julia path/recognize-simple.jl " +
        f"{probfile} {area} WGS84 {outfile}"
    )
    return status

pools = mp.Pool(min(len(yrs), 6))
status = pools.map(recognize, yrs)

#%% frequency
ref = gtiffref(f"{outpath}/classified-{place}-1990-rice-v{version}f-2.tif")

freq = np.zeros([ref["height"], ref["width"]], dtype=np.uint8)

for yr in yrs:
    print(place, yr)
    data, _ = readgtiff(f"{outpath}/classified-{place}-{yr}-rice-v{version}f-2.tif")
    freq += data[0, :, :]

writegtiff(f"{outpath}/frequency-{place}-rice-v{version}f.tif", freq, ref, nthread=20, nodata=0)

#%% filter less year
maxnum = 5

freq, ref = readgtiff(f"{outpath}/frequency-{place}-rice-v{version}f.tif")

def filter_less(yr):
    prob, _ = readgtiff(f"{outpath}/prob-{place}-{yr}-rice-v{version}f-2.tif")
    to_edge = min(yr - 1990, 2023 - yr)#edge
    filter_num = min(maxnum, to_edge)
    print(yr, filter_num)
    mask = freq > filter_num
    prob = prob * mask
    writegtiff(f"{outpath}/prob-{place}-{yr}-rice-v{version}f-3.tif", prob, ref, nthread=20)
    return 0

pools = mp.Pool(min(len(yrs), 5))
status = pools.map(filter_less, yrs)
pools.close()
pools.join()


#%% recognize filter less year
def recognize2(yr):
    print(yr)
    outfile = f"{outpath}/classified-{place}-{yr}-rice-v{version}f-3.tif"
    probfile = f"{outpath}/prob-{place}-{yr}-rice-v{version}f-3.tif"
    area = yaml.safe_load(open(area_file))[place][yr]
    status = os.system(
        "julia path/recognize-simple.jl " +
        f"{probfile} {area} WGS84 {outfile}"
    )
    return status

pools = mp.Pool(min(len(yrs), 6))
status = pools.map(recognize2, yrs)
pools.close()
pools.join()

#%% area time series filter
status = os.system(
    "julia -t 20 path/area-time-series-filter.jl " +
    f"{place} {version}f-3"
)
#%%