#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# AUTHOR: Shen Ruoque
# VERSION: v2023.11.23
#=
Get the full study area of the province,
and rewrite the observation number file by
setting 255 to pixels out the province' range.
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

path = "path/num_30m_files"

for place = places
    days = getdays(place)
    ref = gtiffref("$path/$place/$place-2023-Landsat_num-$days-WGS84.tif")
    full_area = zeros(Bool, ref[:width], ref[:height])
    GC.gc()
    outfile = "$path/$place/$place-full_area-WGS84.tif"
    for yr in yr1:yr2
        println("$place $yr")
        GC.gc()
        file = "$path/$place/$place-$yr-Landsat_num-$days-WGS84.tif"
        data = readgtiff(file)[1][:, :, 1]
        full_area .|= data .≠ 0
    end
    writegtiff(outfile, UInt8.(full_area), ref, nthread=8, nodata=0)

    for yr in yr1:yr2
        println("$place $yr rewrite")
        GC.gc()
        file = "$path/$place/$place-$yr-Landsat_num-$days-WGS84.tif"
        data = readgtiff(file)[1][:, :, 1]
        data[.!full_area] .= 255
        writegtiff(file, UInt8.(data), ref, nthread=8, nodata=255)
    end
end
