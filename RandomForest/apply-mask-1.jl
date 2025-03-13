#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# AUTHOR: Shen Ruoque
# VERSION: v2023.12.20

#=
Mask non-cropland pixels in the observation number file.
(Set 255 to non-cropland pixels)
=#

include("../gtiffio.jl")
yrs = 2016:2023
maskyr = 2020
places = [
    "SouthKorea"
]

cd("path/num_30m_files")

function getdays(place)
    if place in ["SouthKorea"]
        return "137_8_225_day"
    elseif place in [""]
        return ""
    else
        return ""
    end
end


for place in places
    for yr = yrs
        GC.gc()
        println("$place $yr")
        days = getdays(place)
        outfile = "$place/$place-$yr-Landsat_num-$days-WGS84-cropmask.tif"
        if isfile(outfile) continue end
        mask, _ = readgtiff("path/to/LUCC_files.tif")
        data, ref = readgtiff("$place/$place-$yr-Landsat_num-$days-WGS84.tif")
        data[mask .== 0] .= 255
        writegtiff(outfile, data, ref, nthread=8, nodata=255)
        data = mask = nothing
        GC.gc()
    end
end
