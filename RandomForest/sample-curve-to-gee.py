#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR: Shen Ruoque
VERSION: v2024.07.08

把训练曲线合成一个 csv 文件，方便传到 GEE 上
GEE 上训练点多了好像跑不起来，限制到 70000 个好像还能跑起来，不知道能不能更多
"""
#%%
import os
import pandas as pd
import numpy as np
import re
#%%
place = ""
yrs = []
satellite = "S2_L8_L9"; composite = range(1, 361+8, 8)
band = "swir1"
comp_str = f"{composite[0]}_{composite.step}_{composite[-1]}_day"
bands = [f"{i}_{band}" for i in range(len(composite))]
scale = 1000
covers = ["rice", "other_crop"]
codes = {
    "rice": 1,
    "other_crop": 0,
}

#%%
path = "path/sample-curve"
def read_curve(file, cover):
    sample_curve = pd.read_excel(os.path.join(path, file))
    sample_curve.fillna(0, inplace=True)
    sample_curve["class"] = cover
    return sample_curve

curve = pd.DataFrame()
for cover in covers:
    for edge in ["", "-edge"]:
        for yrs in ["x", "2016_2023"]:
        # for yrs in ["2016_2023"]:
            temp = read_curve(
                f"{place}-{yrs}-S2_L8_L9_swir1-1_8_361_day{edge}-{cover}.xlsx",
                codes[cover]
            )
            curve = pd.concat([curve, temp], ignore_index=True)

#%%
picknum = 70000 # 限制曲线数目
curve = curve.iloc[np.random.choice(len(curve), picknum, replace=False)]
curve.to_csv(f"{path}/{place}-S2_L8_L9_swir1-1_8_361_day-x.csv", index=False)
#%%