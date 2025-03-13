#!/usr/bin/env julia
# -*- coding: utf-8 -*-
#=
Identify rice growing areas：
Using two distance graphs, the sorts are added to find the new sort.
When summed, it is weighted according to the number of valid observations of the optical data.
The weight of optical data is low when the quality is poor, and high when the quality is good.

识别水稻种植区：
使用两张距离图，排序相加来求新的排序。
相加的时候，根据光学数据的有效观测数量来加权。
光学数据质量差时权重低，质量好时权重高。
=#
# using Distributed
# addprocs(2)
#=@everywhere=# using YAML, Base.Threads#, SharedArrays

fconfig = length(ARGS) == 0 ? "error-0-.yaml" : ARGS[1]
config = YAML.load(open(joinpath(@__DIR__, fconfig), "r"))

for var = [:yr, :place, :prefix0, :prefix1, :prefix, :namepattern, :versions]
    varname = string(var)
    eval(:($var = config[$varname]))
end

#=@everywhere=# include(joinpath(@__DIR__, "gtiffio.jl"))
#=@everywhere=# include(joinpath(@__DIR__, "recognize_core.jl"))
#=@everywhere=# cropareas = YAML.load(open(joinpath(@__DIR__, "area.yaml"), "r"))

#=@everywhere=# @inbounds function getorder(dist)#根据观测次数valid来定义每个像元的排名
    order = Array{Int32}(undef, size(dist))
    invalid = isnan.(dist) .| (dist .== 0)
    valid = .!invalid
    order[invalid] .= sum(valid) + 1
    invalid = nothing
    GC.gc()
    order[valid] = dist[valid] |> sortperm |> sortperm
    order
end

#=@everywhere=# function getweight!(best, validation, config_swir1)#定义每个像元的权重
    process = config_swir1["process"]
    t = process["t"]
    t_std = process["t_std"]
    tlen, stdlen = length(t), length(t_std)
    tdiff = tlen - stdlen + 1
    corresponds = [(UInt8(i), UInt8(j)) for i = 1:tdiff for j = stdlen+i-1:min(stdlen+i, tlen)]
    corresponds = [(0x0, 0x0); corresponds]

    x = 0:stdlen+1
    weights = @. floor(UInt8, 100 / (1 + exp(-2(7/5*x - 2.5)))) 

    valid = best .≠ 0
    index = findall(valid)
    nvalid = length(index)
    weight = UInt8.(.!valid)
    best = best[index]
    validation = validation[index, :]
    istart = replace(x -> corresponds[x+1][1], best)
    iend = replace(x -> corresponds[x+1][2], best)
    @inbounds @threads for j = 1:nvalid
        weight[index[j]] = weights[sum(validation[j, istart[j]:iend[j]])+1]
    end
    weight
end

#=@everywhere=# function alignment(dist2 :: Matrix{T}, ref2, ref1) where T#地理配准，两个数据集对齐。dist2：距离矩阵，ref2：第2个数据集，ref1：第1个数据集
    res = ref2[:geotransform][2]
    left1, upper1 = ref1[:geotransform][1], ref1[:geotransform][4]
    left2, upper2 = ref2[:geotransform][1], ref2[:geotransform][4]

    left = round(Int64, (left2 - left1) / res + 1)
    upper = round(Int64, (upper1 - upper2) / res + 1)

    ref3 = deepcopy(ref1)
    ref3[:geotransform][2] = ref2[:geotransform][2]
    ref3[:geotransform][6] = ref2[:geotransform][6]
    ref3[:width] = 3ref1[:width]
    ref3[:height] = 3ref1[:height]
    distout = fill(T(NaN), ref3[:width], ref3[:height])
    validwidth = left:min(ref3[:width], left+ref2[:width]-1)
    validheight = upper:min(ref3[:height], upper+ref2[:height]-1)
    distout[validwidth, validheight] = dist2[1:length(validwidth), 1:length(validheight)]
    distout, ref3
end

#=@everywhere=# function alignment0(validation :: Array{T, 3}, refx, ref1) where T#validation（验证数据），它处理的是三维数组，通常用于图像或体数据
    res = ref1[:geotransform][2]
    left1, upper1 = ref1[:geotransform][1], ref1[:geotransform][4]
    left2, upper2 = refx[:geotransform][1], refx[:geotransform][4]

    left = round(Int64, (left2 - left1) / res + 1)
    upper = round(Int64, (upper1 - upper2) / res + 1)

    ref3 = deepcopy(ref1)
    ref3[:width] = ref1[:width]
    ref3[:height] = ref1[:height]

    if (ref1[:width] == refx[:width]) && (ref1[:height] == refx[:height]) && (left == 1) && (upper == 1)
        return validation, refx
    end

    distout = fill(T(0), ref3[:width], ref3[:height], size(validation)[3])
    validwidth = left:min(ref3[:width], left+refx[:width]-1)
    validheight = upper:min(ref3[:height], upper+refx[:height]-1)
    distout[validwidth, validheight, :] = validation[1:length(validwidth), 1:length(validheight), :]
    distout, ref3
end

yrmax = cropareas[place] |> keys |> maximum
version = versions[place]
config_swir1 = config[place][version * "_swir1"]
config_VH = config[place][version * "_VH"]
composite = config_swir1["data"]["composite"]
name_VH = config_VH["data"]["name"]
name_swir1 = config_swir1["data"]["name"]
t = config_swir1["process"]["t"]

path = "the/way/to/distance--VH.tif" #the output path of last step

t0 = collect(composite["start"]:composite["step"]:composite["end"])
istart0 = findfirst(t0 .== t[1])
iend0 = findfirst(t0 .== t[end])

fname = namepattern
fname1 = replace(fname, "\$name" => name_swir1)
fname1 = replace(fname1, "\$version" => version * "_swir1")
fname2 = replace(fname, "\$name" => name_VH)
fname2 = replace(fname2, "\$version" => version * "_VH")
proj = match(r"(ALBERS|WGS84)", fname1)[1]

path_valid = "the/way/to/valid/files"
fvalid = "$place-$yr-$name_swir1.tif"
fvalid = replace(fvalid, "swir1" => "valid")
println(fvalid)

fname_out = "$place-$yr-$proj-v$(version)_2.tif"
outpath = joinpath(path, prefix * fname_out)

if !isfile(outpath)

    println(joinpath(path, prefix0 * fname1))
    dist1, ref1 = readgtiff(joinpath(path, prefix0 * fname1))
    best = readgtiff(joinpath(path, prefix1 * fname1))[1][:, :, 1]
    dist2, ref2 = readgtiff(joinpath(path, prefix0 * fname2))
    dist2, ref2 = alignment(dist2[:, :, 1], ref2, ref1)
    GC.gc()
    validation, refx = readgtiff(joinpath(path_valid, fvalid))
    validation, refx = alignment0(validation[:, :, istart0:iend0], refx, ref1)
    GC.gc()
    validnum = sum(.!isnan.(dist2))
    println("read finish")

    dist1 = dist1[:, :, 1]
    # dist1 = SharedArray(dist1)
    # dist2 = SharedArray(dist2)
    GC.gc()
    order1 = getorder(dist1)
    dist1 = nothing
    GC.gc()
    order2 = getorder(dist2)
    # order1 = fetch(order1)
    # order2 = fetch(order2)
    # dist1 = dist2 = nothing
    dist2 = nothing
    # rmprocs(workers())
    GC.gc()
    order = 9f-2repeat(order1; inner=(3, 3)) .* weight .+ 1f-2order2 .* (1f2 .- weight)
    order1 = order2 = weight = nothing
    GC.gc()
    println("order")

    croparea = 10000cropareas[place][min(yr, yrmax)]
    maxnum = croparea / gridarea(ref2[:geotransform][[2, 6, 4]]...)
    println("$maxnum, $validnum")
    threshold, classified = recognize(
        order, ref2, croparea, min(validnum - 10000, 1.2maxnum); proj=proj, mergesort=false
    )
    println(threshold)
    order = nothing
    GC.gc()
    writegtiff(outpath, UInt8.(classified), ref2)
    classified = nothing
    GC.gc()

end
