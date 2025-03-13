#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# AUTHOR: Shen Ruoque
# VERSION: v2023.07.28

#=
Mask non-cropland pixels in downloaded images using CLCD product.
=#

include("path/gtiffio.jl")
yrs = 2019:2023
place = ""
cd("path/$place")
composite = "_day"

for yr in yrs
    GC.gc()
    println(yr)
    outfile = "$yr/$place-$yr-L5_L9_swir1-$composite-WGS84-crop.tif"
    if isfile(outfile) continue end
    mask, _ = readgtiff("path/LUCC_files.tif")
    data, ref = readgtiff("$yr/$place-$yr-L5_L9_swir1-$composite-WGS84.tif")
    data .*= mask
    writegtiff(outfile, data, ref, nthread=8)
    data = mask = nothing
    GC.gc()
end