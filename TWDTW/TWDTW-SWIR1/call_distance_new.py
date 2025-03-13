#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR: Shen Ruoque
VERSION: v2021.12.27

The new version of the integrated runbook runs
All parameters are written in the configuration file 'method_info.yaml'.
The area per year is written in 'area.yaml'.
1. Read the '.tif' file of the remote sensing image, call 'distance.jl', and calculate the distance value.
2. Call 'recognize.jl' to identify the rice distribution.

新版整合后的运行脚本
所有参数写在配置文件 `method_info.yaml`。
每年的面积写在 `area.yaml`。
1. 读取 遥感影像的 `.tif` 文件，调用 `distance.jl`，计算距离值。
2. 调用 `recognize.jl`，识别水稻分布。
"""
#%%
import yaml, os, time, shutil, re, uuid, json
import multiprocessing as mp
from os.path import join

path = "the/way/to/current/path"
os.chdir(path)
outpath = "the/way/to/output/path"
yrs = [2023]
versions = {
    "SouthKorea" : "0"
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
    print(com)
    data_date = list(range(com["start"], com["end"]+1, com["step"]))
    print(data_date)
    t = config["process"]["t"]
    print(t)
    return {"istart": data_date.index(t[0])+1, "iend": data_date.index(t[-1])+1}

def dist_config(config0, place, version, yr):
    data_info = config0[place][version]["data"]
    name = data_info["name"]
    config = dict(
        config0[place][version]["process"],
        **get_range(config0[place][version]),
        inpath = f"{data_info['path']}/{place}/{yr}",
        outpath = outpath,
        flist = [f"{place}-{yr}-{name}.tif"],
        version = version,
        scale = data_info["scale"],
    )
    fconfig = f"temp-distance-{place}-{yr}-{version}.yaml"
    json_dump(config, join(path, fconfig))
    return fconfig

julia = "/usr/local/julia/julia-1.9.0/bin/julia"

def call_distance(configtuple):
    i, fconfig = configtuple
    time.sleep(30 * i)
    status = os.system(julia + f" -t 12 {path}/distance.jl " + fconfig)
    return status
config0 = yaml_load(join(path, "method_info.yaml"))

#%% get distance using TWDTW
for place in versions.keys():
    version = versions[place]
    print(place)
    print(version)
    print(yrs)
    print(path)
    fconfigs = [dist_config(config0, place, version, yr) for yr in yrs]
    fconfigs = list(zip(range(len(fconfigs)), fconfigs))
    pools = mp.Pool(processes=min(2, len(fconfigs)))
    print(place, time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime()))
    statuses = pools.map(call_distance, fconfigs)
    pools.close()
    pools.join()
    print(statuses)

#%% recognize
config = dict(
    config0,
    places = [re.sub(r"\d", "", place) for place in versions.keys()],
    yrs = yrs,
    prefix0 = "distance-",
    prefix = "classified-",
    namepattern = f"$place-$yr-$name-v$version.tif",
    versions = versions,
)
fconfig = f"temp-recognize-{uuid.uuid1()}.yaml"
json_dump(config, join(path, fconfig))
catch_error(
    julia +                                                                                   
    f" -p 2 {path}/recognize.jl {fconfig}", 
    join(path, fconfig), ""
)
print(time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime()))
#%%