#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# AUTHOR: Shen Ruoque
# VERSION: v2023.12.21
#=
遥感影像中某些像元在该年可能一次有效观测都没有。
故而这些像元无法用 RF 识别得到概率值，
利用临近年份的概率值对这些像元进行填补。
=#

using Base.Threads
using Statistics
include("/gtiffio.jl")

place = ARGS[1]
version = ARGS[2]

numpath = ""
probpath = ""


yrs = 1990:2023
yrnum = length(yrs)

function getdays(place)
    if place in ["SouthKorea"]
        return "137_8_225_day"
    elseif place in [""]
        return ""
    else
        return ""
    end
end

days = getdays(place)

numfile = joinpath(numpath, "$place/$place-LS_num-$days-obs_eq_0-WGS84.tif")
nums, ref = readgtiff(numfile)


probfiles = [joinpath(probpath, "prob-$place-$yr-rice-v$(version)-WGS84-2.tif") for yr in yrs]
probs = stack([readgtiff(probfile)[1] for probfile in probfiles], dims=4)
probtype = eltype(probs)

println("read finish.") 

failpixs = nums .== 1
validpixs = .!failpixs .& (nums .!= 255)

pixcels = findall(any(failpixs, dims=3)[:, :, 1])

@threads for pixcel in pixcels
    fails = findall(failpixs[pixcel, :])
    valids = findall(validpixs[pixcel, :])
    for fail in fails
        formers = valids[valids .< fail]
        laters = valids[valids .> fail]
        if length(formers) != 0
            former = maximum(formers)
            if length(laters) != 0
                later = minimum(laters)
                average = mean([probs[pixcel, :, former], probs[pixcel, :, later]])
                if probtype <: Integer
                    probs[pixcel, :, fail] = round.(probtype, average)
                else
                    probs[pixcel, :, fail] = average
                end
            else
                probs[pixcel, :, fail] = probs[pixcel, :, former]
            end
        else
            if length(laters) != 0
                later = minimum(laters)
                probs[pixcel, :, fail] = probs[pixcel, :, later]
            end
        end
    end
end

println("fill finish.")

for yr in yrs
    println(yr)
    writegtiff(
        joinpath(probpath, "prob-$place-$yr-rice-v$(version)f-2.tif"),
        probs[:, :, :, yr-yrs[1]+1], ref, nthread=10
    )
end
