#!/usr/bin/env julia
# -*- coding: utf-8 -*-
#=
输入 RF 预测的水稻概率，卡面积阈值。
=#
infile = ARGS[1]
crop_area = parse(Float64, ARGS[2])
proj = ARGS[3]
outpath = ARGS[4]

include("path/gtiffio.jl")
include("path/recognize_core.jl")

dist, ref = readgtiff(infile)
threshold, classified = recognize(dist[:, :, 1], ref, 10000crop_area; proj=proj, rev=true)
println("$outpath $threshold")
writegtiff(outpath, UInt8.(classified), ref; nthread=40)
