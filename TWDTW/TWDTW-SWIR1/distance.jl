#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# AUTHOR: Shen Ruoque
#=
Read the SWIR1 image file, calculate the distance using the TWDTW algorithm, 
and output to the file `distance-....tif`.

读取 SWIR1 影像文件，
使用 TWDTW 算法计算距离，
输出到文件 `distance-....tif` 中。
=#
using YAML, Dates
using Base.Threads
evaluate = eval ∘ Meta.parse
include(joinpath(@__DIR__, "gtiffio.jl"))
include("./TWDTW.jl")

foundwhere = false

fconfig = length(ARGS) == 0 ? "temp-distance.yaml" : ARGS[1]
config = YAML.load(open(joinpath(@__DIR__, fconfig), "r"))

for varname = keys(config)
    var = Symbol(varname)
    eval(:($var = config[$varname]))
end
std = Float32.(std)
scale = Float32(scale)
method = evaluate(method)
evaluate(get_dist)

for file = flist

    infile = joinpath(inpath, file)
    outfile = joinpath(outpath, "distance-" * replace(
        file, r"\.tif" => s"-v" * "$version.tif"
    ))

    if !isfile(infile) || isfile(outfile) continue end
    println(file)
    
    data, ref0 = readgtiff(infile)
    valid = dropdims(any(data .≠ 0, dims=3); dims=3)
    index = findall(valid)
    validnum = length(index)
    if validnum == 0 writegtiff(outfile, dist, ref0); continue end
    dist = .!valid .* NaN32
    data = data[index, istart:iend]

    if foundwhere == false
        @inbounds @threads for j = 1:validnum
            dist[index[j]] = getdist(data[j, :], t, std)
        end
    else
        outbest = joinpath(outpath, "best-" * replace(
            file, r"\.tif" => s"-v" * "$version.tif"
        ))
        best = UInt8.(valid)
        @inbounds @threads for j = 1:validnum
            dist0, best0 = getdist(data[j, :], t, std)
            dist[index[j]], best[index[j]] = dist0, UInt8.(best0)
        end
        writegtiff(outbest, best, ref0, nthread=12)
    end

    writegtiff(outfile, dist, ref0, nthread=12)

end
