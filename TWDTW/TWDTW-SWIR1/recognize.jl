#!/usr/bin/env julia
# -*- coding: utf-8 -*-
#=
Identify rice planting areas.
识别水稻种植区
=#
using Distributed
@everywhere using YAML

fconfig = length(ARGS) == 0 ? "temp-recognize.yaml" : ARGS[1]
config = YAML.load(open(joinpath(@__DIR__, fconfig), "r"))

@everywhere include(joinpath(@__DIR__, "gtiffio.jl"))
@everywhere include(joinpath(@__DIR__, "recognize_core.jl"))

@everywhere maize_areas = YAML.load(open(joinpath(@__DIR__, "area.yaml"), "r"))

for var = [:yrs, :places, :prefix0, :prefix, :namepattern, :versions]
    varname = string(var)
    eval(:($var = config[$varname]))
end

for place in places
    yrmax = maize_areas[place] |> keys |> maximum
    version = versions[place]
    name = config[place][version]["data"]["name"]
    @sync @distributed for yr in yrs
        path = "the/way/to/distance-.tif" # the output path in last step
        fname = namepattern
        fname = replace(fname, "\$yr" => yr)
        fname = replace(fname, "\$place" => place)
        fname = replace(fname, "\$name" => name)
        fname = replace(fname, "\$version" => version)
        proj = match(r"(ALBERS|WGS84)", fname)[1]

        infile = joinpath(path, prefix0 * fname)
        dist, ref = readgtiff(infile)
        dist = dropdims(dist; dims=3)

        maize_area = maize_areas[place][min(yr, yrmax)]
        threshold, classified = recognize(dist, ref, 10000maize_area; proj=proj)
        println("$place: $yr: $threshold")
        outpath = joinpath(path, prefix * fname)
        writegtiff(outpath, UInt8.(classified), ref)
    end
end
rmprocs(workers())