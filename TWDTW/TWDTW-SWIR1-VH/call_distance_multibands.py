#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR: Shen Ruoque
VERSION: v2021.12.27

The new version of the integrated runtime script is used to identify the VH and swir1 bands.
All parameters are written in the configuration file 'method_info.yaml'.
The area per year is written in 'area.yaml'.
1. Read the '.tif' file of the remote sensing image, call 'distance.jl', and calculate the distance value.
2. Call 'recognize_multibands.jl' to identify the rice distribution.

新版整合后的运行脚本，用于结合 VH 和 swir1 两个波段进行识别。
所有参数写在配置文件 `method_info.yaml`。
每年的面积写在 `area.yaml`。
1. 读取 遥感影像的 `.tif` 文件，调用 `distance.jl`，计算距离值。
2. 调用 `recognize_multibands.jl`，识别水稻分布。
"""
#%%
import yaml, os, time, shutil, re, uuid, json
import multiprocessing as mp
from os.path import join

path = "the/way/to/current/path"
os.chdir(path)
outpath = "the/way/to/ouput/path"
yrs = [2016]
versions = {
    "Thailand" : "0",
}

def catch_error(script, tempfile, err_str):
    for i in range(3):
        status = os.system(script)
        if status != 0:
            print(f"Error calling julia {i}. Code: {status >> 8}")
            shutil.copy(tempfile, join(path, f"error-{i}-" + err_str + ".yaml"))
        else:
            os.remove(tempfile)
            break

def json_dump(config, file):
    with open(file, "w", encoding="utf-8") as f:
        json.dump(config, f, ensure_ascii=False, indent=2)

def yaml_load(file):
    with open(file, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def get_range(config):
    com = config["data"]["composite"]
    data_date = list(range(com["start"], com["end"]+1, com["step"]))
    t = config["process"]["t"]
    return {"istart": data_date.index(t[0])+1, "iend": data_date.index(t[-1])+1}

def get_multibands(config0, place, version):
    multibands = []
    for key in config0[place].keys():
        band = re.match(f"^{version}(_.+)?$", key)
        if band != None:
            if band[1] == "_VH":
                multibands.append([key, 3])
            else:
                multibands.append([key, 6])
    return multibands

def dist_config(config0, place, version, yr):
    data_info = config0[place][version]["data"]
    config = dict(
        config0[place][version]["process"],
        **get_range(config0[place][version]),
        inpath = f"{data_info['path']}/{place}/{yr}",
        outpath = outpath,
        flist = [f"{place}-{yr}-{data_info['name']}.tif"],
        version = version,
        scale = data_info["scale"],
    )
    fconfig = f"temp-distance-{place}-{yr}-{version}.yaml"
    json_dump(config, join(path, fconfig))
    return fconfig

julia = "/usr/local/julia/julia-1.9.0/bin/julia"

def call_distance(configtuple):
    i, fconfig = configtuple
    time.sleep(60 * i)
    status = os.system(julia + f" -t 15 {path}/distance.jl " + fconfig)
    return status

def call_recognize(configs):
    fconfig = configs
    status = os.system(julia + f" {path}/recognize_multibands.jl " + fconfig)
    return status

config0 = yaml_load(join(path, "method_info.yaml"))

#%% get distance using TWDTW
for place in versions.keys():
    version = versions[place]
    for version, prcs in get_multibands(config0, place, version):
        fconfigs = [dist_config(config0, place, version, yr) for yr in yrs]
        fconfigs = list(zip(range(len(fconfigs)), fconfigs))
        pools = mp.Pool(processes=1)
        print(place, time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime()))
        statuses = pools.map(call_distance, fconfigs)
        pools.close()
        pools.join()
        print(statuses)

#%% recognize
print(time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime()))
configs = []
for place in versions.keys():
    version = versions[place]
    for yr in yrs:
        config = dict(
            config0,
            place = re.sub(r"\d", "", place),
            yr = yr,
            prefix0 = "distance-",
            prefix1 = "best-",
            prefix = "classified-",
            namepattern = f"{place}-{yr}-$name-v$version.tif",
            versions = versions,
        )
        fconfig = f"temp-recognize-{place}-{yr}-{version}.yaml"
        json_dump(config, join(path, fconfig))
        
        configs.append(fconfig)

pools = mp.Pool(processes=1)
statuses = pools.map(call_recognize, configs)
pools.close()
pools.join()
print(statuses)

#%%