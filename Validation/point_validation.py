#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR: Shen Ruoque
VERSION: v0.1.0
"""
#%%
import os
import numpy as np
import pandas as pd
import yaml
import glob
from sklearn.metrics import confusion_matrix
from sklearn.metrics import cohen_kappa_score
#%%
places = [
    "SouthKorea-2018",
    "SouthKorea-2019"
]
version = "0"
os.chdir("path")

valid = {"acc": {}, "kappa": {}}
for place in places:
    fnames = glob.glob(f"{place}-point_validation-0.csv") #用于计算混淆矩阵的结果csv文件
    sheets = [pd.read_csv(fname, encoding="utf-8") for fname in fnames]
    validation = pd.concat(sheets, ignore_index=True)

    cm = confusion_matrix(validation["class"] != 0, validation["b1"] != 1)
    acc_p = [cm[i, i] / sum(cm[i, :]) * 100 for i in range(2)]
    acc_u = [cm[i, i] / sum(cm[:, i]) * 100 for i in range(2)]

    acc = sum([cm[i, i] for i in range(2)]) * 100 / sum(sum(cm))
    kappa = cohen_kappa_score(validation["class"] != 0, validation["b1"] != 1)
    print(
        # "%s:\n\t\tNon-Rice\tRice\n"%place +
        "%s\t%s\tRice\t%i\t%i\n"%(
            version, place, cm[0, 0], cm[0, 1]
        ) +
        "\t\t\tOther\t%i\t%i\n"%(
            cm[1, 0], cm[1, 1]
        )# +
        # "\nProductor Accuracy: Non-Rice: %2.2f%%, Rice: %2.2f%%\n"%(acc_p[0], acc_p[1]) +
        # "User Accuracy: Non-Rice: %2.2f%%, Rice: %2.2f%%\nTotal Accuracy: %2.2f%%, Kappa: %2.2f"%(
        #     acc_u[0], acc_u[1], acc, kappa
        # )
    )
    valid["acc"][place] = float(acc)
    valid["kappa"][place] = float(kappa)
    print(acc)
    print(kappa)

# with open("point_validation-10_70.yaml", "w", encoding="utf-8") as f:
#     yaml.dump(valid, f) 

#%%