#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# AUTHOR: Shen Ruoque
# VERSION: v2021.09.09
#=
Extract the central coordinates of pixel points for rice and non-rice samples and write them as a CSV.
提取水稻和非水稻样本像元点的中心坐标，写为 csv
=#

using Printf, JSON
import ArchGDAL as AG
using DataFrames, CSV
cd(@__DIR__)
include("./gtiffio.jl")


function coords2geojson(fname, coordinates0, class; proj="WGS84")
    prestring = """{"type":"Feature","properties":{"class":"$class"},""" *
        """"geometry":{"type":"Point","coordinates":"""
    open(fname, "w") do f
        write(f, """{"type":"FeatureCollection","features":[\n""")
        if proj == "WGS84"
            coordinates = [@sprintf "[%4.15f, %4.15f]" i[1] i[2] for i = coordinates0]
        elseif proj == "ALBERS"
            coordinates = [@sprintf "[%24.15f, %24.15f]" i[1] i[2] for i = coordinates0]
        end
        for (i, coordinate) = enumerate(coordinates[1:end-1])
            write(f, prestring * "$coordinate}},\n")
        end
        write(f, prestring * "$(coordinates[end])}}\n]}")
    end
    fname
end

function albers2wgs84(infie, outfile)
    albers = "+proj=aea +lat_0=0 +lon_0=105 +lat_1=25 +lat_2=47" *
        " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
    run(`'C:\\Users\\DELL\\anaconda3\\Library\\bin\\ogr2ogr.exe' -s_srs $albers -t_srs 'EPSG:4296' $outfile $infie`)
    outfile
end

function geojson2csv(infile, outfile, class)
    features = open(infile, "r") do f
        JSON.parse(f)["features"]
    end
    open(outfile, "w") do f
        write(f, "class,longitude,latitude\n")
        for i in features
            lonlat = i["geometry"]["coordinates"]
            write(f, "$class,$(lonlat[1]),$(lonlat[2])\n")
        end
    end
    outfile
end

"""
    coords2list(fname0, coordinates, class)
The obtained coordinates will be converted from ALBERS to WGS84 and written in geojson format.
将得到的坐标从 ALBERS 转换成 WGS84，写成 geojson。
"""
function coords2list(fname0, coordinates0, class; proj="WGS84")
    fname = coords2geojson(fname0 * ".geojson", coordinates0, class; proj=proj)
    if proj == "ALBERS"
        fname2 = albers2wgs84(fname, fname0 * "0.geojson")
        if isfile(fname) rm(fname) end
    elseif proj == "WGS84"
        fname2 = fname
    end
    fname3 = geojson2csv(fname2, fname0 * ".csv", class)
    fname3
end

"""
    mergecsv(flist, fout)
The obtained CSV files for rice and non-rice will be merged and sorted.
将得到的水稻和非水稻的 csv 文件拼合起来，并排序。
"""
function mergecsv(flist, fout)
    s1 = open(flist[1], "r") do f
        readlines(f)
    end
    open(fout, "w") do f
        [write(f, i * "\n") for i in s1]
        for file = flist[2:end]
            s2 = open(file, "r") do f2
                readlines(f2)[2:end]
            end
            [write(f, i * "\n") for i in s2]
        end
    end
    fout
end


path = "the/way/to/kml/files"
cd(path)
placeyr = ARGS[1]
flist = []
for (class, nclass) = [("rice", 0), ("other", 1)]
    fname = "$placeyr-$class-30m"
    fout = "$placeyr-$class-30m"
    if isfile(fout * ".geojson") rm(fout * ".geojson") end

    infile = fname * ".tif"
    data, ref = readgtiff(infile)
    width = ref[:width]
    geotrans = ref[:geotransform]
    data = dropdims(data; dims=3)[:]
    rice_id = findall(data .> 90)

    picked = [[
        geotrans[1] + ((i - 1) % width + 0.5) * geotrans[2],
        geotrans[4] + ((i - 1) ÷ width + 0.5) * geotrans[6]
    ] for i = rice_id]

    push!(flist, coords2list(fout, picked, nclass; proj="WGS84"))
end
mergecsv(flist, "$placeyr-30m.csv")
for file in flist rm(file) end
