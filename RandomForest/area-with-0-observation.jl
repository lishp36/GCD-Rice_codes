#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# AUTHOR: Shen Ruoque
# VERSION: v2023.11.23
#=
Extract pixels with o observations of each province each year.
=#

include("../gtiffio.jl")

places = [
    "SouthKorea"
]

function getdays(place)
    if place in ["SouthKorea"]
        return "137_8_225_day"
    elseif place in [""]
        return ""
    else
        return ""
    end
end

function filter_number(x)
    if x == 0
        return 0x01
    elseif x == 255
        return 0xff
    else
        return 0x00
    end
end

path = "path/to/num_30m_files"
yrs = yr1:yr2

for place = places
    days = getdays(place)
    outfile = "$path/$place/$place-Landsat_num-$days-obs_eq_0-WGS84.tif"
    ref = gtiffref("$path/$place/$place-2023-Landsat_num-$days-WGS84-cropmask.tif")
    low_obs_areas = zeros(UInt8, ref[:width], ref[:height], length(yrs))
    GC.gc()
    for yr = yr1:yr2
        println("$place $yr")
        GC.gc()
        file = "$path/$place/$place-$yr-Landsat_num-$days-WGS84-cropmask.tif"
        temp = readgtiff(file)[1][:, :, 1]
        replace!(filter_number, temp)
        low_obs_areas[:, :, yr-yrs[1]+1] = temp
        temp = nothing
    end
    temp = sum(low_obs_areas, dims=3)
    println("max year ", maximum(temp[temp .<= 27]))
    writegtiff(outfile, low_obs_areas, ref, nthread=8, nodata=0)
end
