#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Draw a verification diagram, one containing multiple diagrams
画验证图，一张包含多个图
"""
#%%
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
from matplotlib import rcParams

#%%
rcParams["font.family"] = "Times New Roman"
os.chdir("the/way/to/the/staticstic_identified_area_files/")#识别数据.xlsx  路径

def varify(x0, y0, accuracy):
    x = []; y = []
    for i in range(len(x0)):
        if np.isnan(x0[i]) or np.isnan(y0[i]): continue
        x.append(x0[i])
        y.append(y0[i])
    rng = max(
        np.nanmax(x), np.nanmax(y)
    ) // accuracy * accuracy + accuracy
    data = pd.DataFrame({"x": x, "y": y})
    sns.regplot(
        data=data, x="x", y="y", marker=".", color="k", 
        scatter_kws={"s": 100, "lw": 0},
        line_kws={"ls": (0,(5,10)), "lw": 0.85}
    )
    plt.xlabel(""); plt.ylabel("")
    k, b = np.linalg.lstsq(np.vstack([x, np.ones(len(x))]).T, y, rcond=-1)[0]
    r, _ = stats.pearsonr(x, y)
    dtm = r ** 2
    # rmse = np.sqrt(np.mean((np.array(x) - np.array(y)) ** 2))
    mae = np.mean(np.abs(np.array(x) - np.array(y)))
    rmae = mae / np.mean(np.array(x))
    # mrae = np.mean(np.abs(np.array(x) - np.array(y)) / np.array(x))
    plt.plot([0, rng], [0, rng], color="k", lw=0.85)
    plt.xlim(0, rng); plt.ylim(0, rng)
    plt.xticks(size="x-large")
    plt.yticks(size="x-large")
    est_equ = "y = %1.2f x "%k + ("+ %.2f\n"%b if b > 0 else "− %5.2f\n"%-b)
    plt.text(
        rng - rng / 32, rng / 32,
        est_equ +
        "R$^2$ = %1.3f\n"%dtm +
        "MAE = %.2f\n"%mae +
        "RMAE = %.2f"%rmae,# +
        # "MRAE = %.2f"%mrae,
        horizontalalignment='right',
        verticalalignment='bottom',
        size="x-large"
    )
    ax = plt.gca()
    ax.set_aspect(1)
    # plt.grid()
    return rng
#%%
def varify_1yr(place, yr):
    #ricetype0 = "srice" if ricetype == "spr" else "lrice"
    fig = plt.figure(figsize=(5, 5), dpi=100, facecolor="white")
    ax = fig.add_subplot(1, 1, 1)
    xlsx = pd.read_excel(f"{place}-area-valisation.xlsx") #面积.xlsx  名称
    print(xlsx)
    x = xlsx[f"{yr}"].to_numpy() / 1000 # 1000ha -> ha   列名称
    y = xlsx[f"boro_rf27-2_{yr}"].to_numpy() / 1000   #列名称
    _ = varify(x, y, accuracy)
    # plt.suptitle(version, size=20)
    plt.savefig(f"{place}.png")
    plt.tight_layout()
    plt.close()
    return 0
#%%
def varify_3yr(place, yrs):
    # ricetype0 = "srice" if ricetype == "spr" else "lrice"
    ntot = len(yrs)
    nx = min(5, ntot)
    ny = (ntot - 1) // 5 + 1
    fig = plt.figure(figsize=(5*nx, 5*ny), dpi=100, facecolor="white")
    for iyr, yr in enumerate(yrs):
        ax = fig.add_subplot(ny, nx, iyr+1)
        xlsx = pd.read_excel(f"{place}-area-validation.xlsx") # 输入面积文件.xlsx
        print(xlsx)
        x = xlsx[f"{yr}(ha)"].to_numpy() / 1000 # 1000ha -> ha  列名称
        y = xlsx[f"aus_rf27-2_{yr}"].to_numpy() / 1000   # 列名称
        max_axis = varify(x, y, accuracy)
        plt.text(
            max_axis / 16, 0.90 * max_axis,
            "(" + chr(ord("a")+iyr) + f") {yr}",
            size="xx-large"
        )
        fig.text(
            0.5, 0.06, "Statistical area (×10³ha)",
            size="xx-large", ha="center", va="center"
        )
        fig.text(
            0.06, 0.5, "Identified area (×10³ha)",
            size="xx-large", ha="center", va="center", rotation=90
        )
    plt.suptitle(place, size=20)
    plt.savefig(f"{prefix}{place}-area-validation-compare-NESEA.png")
    plt.tight_layout()
    plt.close()
    return 0

#%%
places = {
    "SouthKorea" : list(range(1990,2023+1))
}
prefix = "classified-"
# level = "city"; accuracy = 1e2
#level = "county"; 
accuracy = 1e1
for place in places.keys():
    # for version in ["rf27-2"]:
    print(place)
    print(places)
    # varify_1yr(place,2014)
    varify_3yr(place, places[place])
    #varify_3yr(place, version, places[place], "double")
print("完成制图")

#%%