#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR: Shen Ruoque
VERSION: v2022.01.22
Use rasterio.mask.mask to segment the identification map of each country,
The map of each province/city was cropped out and the rice area of each province/city was calculated for verification.
Output Excel.

使用 rasterio.mask.mask 对各国家识别图进行分割，
裁剪出每个省/市的图，计算各省/市水稻面积用于验证。
输出 Excel。
"""
#%%
import os
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import rasterio.mask
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore")
#%%
def gridarea(lonres, latres, lat):
    return (
        6371008.8 ** 2 * np.deg2rad(lonres) *
        (np.sin(np.deg2rad(lat)) - np.sin(np.deg2rad(lat + latres)))
    )
gridareavec = np.vectorize(gridarea)

def getarea(place, infile, flag, res=None):
    """
    计算省/市水稻面积

    Parameters
    ----------
    place: GeoDataFrame row, 省/市级矢量边界
    infile: Str, 省级识别结果文件
    flag: List, 水稻标识
    res: Float, 栅格分辨率

    Returns
    -------
    area: Float, 省/市水稻面积列表
    """
    geo = [place.geometry.__geo_interface__]
    with rasterio.open(infile) as ds:
        try:
            data, transform = rasterio.mask.mask(ds, geo, crop=True)
        except ValueError:
            return np.nan
    data = data[0, :, :]
    height = data.shape[0]
    pixnum = np.zeros(data.shape[0], dtype=np.uint16)
    for iflag in flag:
        pixnum += np.sum(data == iflag, axis=1, dtype=np.uint16)
    del data
    if res is None:
        lonres = transform[0]
        latres = transform[4]
        lats = np.array([transform[5] + i * latres for i in range(height)])
        pixareas = gridareavec(lonres, latres, lats)
        _grid_area = np.sum(pixnum * pixareas)
        if _grid_area > 10000:
            area = _grid_area / 10000
        else:
            area = np.nan
        del lats, pixareas
    else:
        _grid_num = np.sum(pixnum) * res * res
        if _grid_num > 10000:
            area = _grid_num / 10000
        else:
            area = np.nan
    del pixnum
    return area
#%%
version = "rf27-2"
prefix = "classified-"

regions = {
    "SouthKorea" : "SouthKorea"
}

inpath = "the/way/to/the/rice/maps/" #识别结果影像路径
vecpath = f"the/way/to/the/boundry_shapefiles.shp" #shp边界
validpath = "the/way/to/the/staticstic_area_files/" #输入.xlsx路径

for province in regions.keys():
    print(province)
    places = gpd.read_file(vecpath)
    sheet = pd.read_excel(os.path.join(
        validpath, f"{regions[province]}-area.xlsx" #输入.xlsx 名称
    ))
    for yr in range(1990, 2023+1):
        print(yr) 
        file = (
            f"{province}_{yr}.tif"
        ) #识别结果影像名称
        infile = os.path.join(inpath, file)
        print(infile)
        with rasterio.open(infile)as ds:
            crs = ds.crs
        places_infile = [
            (place, infile)
            for _, place in places.to_crs(crs).iterrows()
            if place["ID"] in list(sheet["ID"])
        ]

        sheet2 = pd.DataFrame({"ID": [i[0]["ID"] for i in places_infile]})

        def getarea_p1(place_infile):
            return getarea(*place_infile, [1]) # res=20.0 #[1]代表计算标记值为1的地方
        # def getarea_p2(place_infile):
        #     return getarea(*place_infile, [2]) # res=20.0
        # def getarea_p3(place_infile):
        #     return getarea(*place_infile, [3]) # res=20.0

        pools = mp.Pool(processes=40)
        areas = pools.map(getarea_p1, places_infile)
        sheet2[f"single_{version}_{yr}"] = pools.map(getarea_p1, places_infile)
        # sheet2[f"double_{version}_{yr}"] = pools.map(getarea_p2, places_infile)
        # sheet2[f"triple_{version}_{yr}"] = pools.map(getarea_p3, places_infile)

        pools.close()
        pools.join()
        sheet = sheet.merge(sheet2, on="ID", how="left")

sheet.to_excel(
        os.path.join(validpath, f"{province}-area-validation.xlsx"),#输出结果.xlsx 名称
        encoding="utf-8", index=False
    )
print("完成统计")

#%%