#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# AUTHOR: Shen Ruoque
# VERSION: 2022.03.10
#=
TWDTW.jl
使用 Julia 和 Fortran 重写了两个版本的 TWDTW 算法。

序列长度为 10，Float32 时：

Linux 工作站下：
Julia 版 twdtw 耗时约为 1.8 μs
Fortran 版 twdtw 耗时约为 1.2 μs
MATLAB 版 twdtw 耗时约为 6.4 μs
Julia 版 tw_dtw 耗时约为 424 ns
Fortran 版 tw_dtw 耗时约为 171 ns

Windows 台式机下：
Julia 版 twdtw 耗时约为 1.5 μs
Fortran 版 twdtw 耗时约为 4.0 μs
MATLAB 版 twdtw 耗时约为 7.2 μs
Julia 版 tw_dtw 耗时约为 414 ns
Fortran 版 tw_dtw 耗时约为 165 ns
=#
module TWDTW

export twdtw, dtw, tw, tw_dtw, tw_m_dtw

using Libdl

"""
    _tw_julia(rdate, tdate, α, β)

Using @fastmath, The time series should not contain NaNs or Infs.

version: 2022.03.10
"""
function _tw_julia end

"""
    _tw_dtw_julia(r, t, tw)

Using @fastmath, The time series should not contain NaNs or Infs.

version: 2022.03.10
"""
function _tw_dtw_julia end

"""
    _twdtw_julia(r, t, rdate, tdate, α, β)

Using @fastmath, The time series should not contain NaNs or Infs.

original MATLAB code authored by Dong Jie

convert to Julia: Shen Ruoque

version: 2022.03.10
"""
function _twdtw_julia end

"""
    _dtw_julia(r, t)

Using @fastmath, The time series should not contain NaNs or Infs.

version: 2022.03.10
"""
function _dtw_julia end

"""
    _tw_m_dtw_julia(r, t, tw)

Using @fastmath, The time series should not contain NaNs or Infs.

version: 2022.03.10
"""
function _tw_m_dtw_julia end

@fastmath @inbounds function dtw_core!(d)
    M, N = size(d)
    for m = 2:M
        d[m, 1] += d[m-1, 1]
    end
    for n = 2:N
        d[1, n] += d[1, n-1]
        for m = 2:M
            d[m, n] += min(d[m-1, n], d[m-1, n-1], d[m, n-1])
        end
    end
    d[M, N]
end

@fastmath @inbounds _dist(r, t) = @. abs(r - t')

@fastmath @inbounds function _tw_julia(rdate, tdate, α, β)
    time_diff = @. abs(rdate - tdate')
    @. 1 / (1 + exp(-α * (time_diff - β))) # -0.1, 50 % A Logistic Curve
end

_tw_julia(rdate, tdate; α=0.1, β=50) = _tw_julia(rdate, tdate, α, β)

_tw_dtw_julia(
    r :: Vector{T}, t :: Vector{T}, tw :: Matrix{T}
) where T <: AbstractFloat = dtw_core!(_dist(r, t) .+ tw)

function _tw_dtw_julia(r :: Vector{T}, t :: Vector{T}, tw) where T <: Integer
    _tw_dtw_julia(
        convert(Vector{Float64}, r), convert(Vector{Float64}, t),
        convert(Matrix{Float64}, tw)
    )
end

_tw_dtw_julia(r, t, tw; factor) = _tw_dtw_julia(r ./ factor, t, tw)

_twdtw_julia(
    r :: Vector{T}, t :: Vector{T}, rdate, tdate, α, β
) where T <: AbstractFloat = _tw_dtw_julia(
    r, t, convert(Matrix{T}, _tw_julia(rdate, tdate; α, β))
)

_twdtw_julia(
    r :: Vector{T}, t :: Vector{T}, rdate, tdate, α, β
) where T <: Integer = _tw_dtw_julia(
    convert(Vector{Float64}, r), convert(Vector{Float64}, t),
    convert(Matrix{Float64}, _tw_julia(rdate, tdate; α, β))
)

_twdtw_julia(r, t, rdate, tdate; α=0.1, β=50) = _twdtw_julia(r, t, rdate, tdate, α, β)

_dtw_julia(
    r :: Vector{T}, t :: Vector{T}
) where T <: AbstractFloat = dtw_core!(_dist(r, t))

function _dtw_julia(r :: Vector{T}, t :: Vector{T}) where T <: Integer
    _dtw_julia(convert(Vector{Float64}, r), convert(Vector{Float64}, t))
end

_dtw_julia(r, t; factor) = _dtw_julia(r ./ factor, t)

_tw_m_dtw_julia(
    r :: Vector{T}, t :: Vector{T}, tw :: Matrix{T}
) where T <: AbstractFloat = dtw_core!(_dist(r, t) .* tw)

function _tw_m_dtw_julia(r :: Vector{T}, t :: Vector{T}, tw) where T <: Integer
    _tw_m_dtw_julia(
        convert(Vector{Float64}, r), convert(Vector{Float64}, t),
        convert(Matrix{Float64}, tw)
    )
end

_tw_m_dtw_julia(r, t, tw; factor) = _tw_m_dtw_julia(r ./ factor, t, tw)

# Fortran

"""
    _twdtw_fortran(r, t, rdate, tdate, α, β)

Using Fortran, The time series should not contain NaNs or Infs.

original MATLAB code authored by Dong Jie

convert to Fortran: Shen Ruoque

version: 2022.03.10
"""
function _twdtw_fortran end

"""
    _dtw_fortran(r, t)

Using Fortran, The time series should not contain NaNs or Infs.

version: 2022.03.10
"""
function _dtw_fortran end

"""
    _tw_fortran(r, t)

Calculate Timw weight matrix using Fortran

version: 2022.03.10
"""
function _tw_fortran end

"""
    _tw_dtw_fortran(r, t, tw)

Using Fortran, The time series should not contain NaNs or Infs.

version: 2022.03.10
"""
function _tw_dtw_fortran end

"""
    _tw_m_dtw_fortran(r, t, tw)

Using Fortran, The time series should not contain NaNs or Infs.

version: 2022.03.10
"""
function _tw_m_dtw_fortran end


if !isfile(joinpath(@__DIR__, "twdtw.so")) &&
    isfile(joinpath(@__DIR__, "twdtw.f90"))
    run(`gfortran ./twdtw.f90 -Ofast -shared -fPIC -o ./twdtw.so`)
end

lib = Libdl.dlopen(joinpath(@__DIR__, "twdtw.so"))

libtwdtw = Dict(
    :Float32 => Dict(
        :Int32 => Libdl.dlsym(lib, :__twdtw_MOD_single32),
        :Int64 => Libdl.dlsym(lib, :__twdtw_MOD_single64),
    ),
    :Float64 => Dict(
        :Int32 => Libdl.dlsym(lib, :__twdtw_MOD_double32),
        :Int64 => Libdl.dlsym(lib, :__twdtw_MOD_double64),
    ),
)

for xtype in [:Float32, :Float64]
    for datetype in [:Int32, :Int64]
        eval(quote
            function _twdtw_fortran(
                r :: Vector{$xtype}, t :: Vector{$xtype},
                rdate :: Vector{$datetype}, tdate :: Vector{$datetype},
                α :: $xtype, β :: $xtype
            ) :: $xtype
                ccall(
                    $(libtwdtw[xtype][datetype]), $xtype,
                    (
                        Ref{$xtype}, Ref{$datetype}, Ref{Int32},
                        Ref{$xtype}, Ref{$datetype}, Ref{Int32},
                        Ref{$xtype}, Ref{$xtype}
                    ),
                    r, rdate, length(r), t, tdate, length(t), α, β
                )
            end
        end)
    end
end

function _twdtw_fortran(r :: Vector{T}, t :: Vector{T}, rdate, tdate; α=0.1, β=50) where T <: Integer
    _twdtw_fortran(
        convert(Vector{Float64}, r), convert(Vector{Float64}, t),
        rdate, tdate, Float64(α), Float64(β)
    )
end

function _twdtw_fortran(r :: Vector{T}, t :: Vector{T}, rdate, tdate; α=0.1, β=50) where T <: AbstractFloat
    _twdtw_fortran(r, t, rdate, tdate, T(α), T(β))
end

libdtw = Dict(
    :Float32 => Libdl.dlsym(lib, :__dtw_MOD_singlex),
    :Float64 => Libdl.dlsym(lib, :__dtw_MOD_doublex),
)

for xtype = [:Float32, :Float64]
    eval(quote
        function _dtw_fortran(
            r :: Vector{$xtype}, t :: Vector{$xtype}
        ) :: $xtype
            ccall(
                $(libdtw[xtype]), $xtype,
                (Ref{$xtype}, Ref{Int32}, Ref{$xtype}, Ref{Int32}),
                r, length(r), t, length(t)
            )
        end
    end)
end

function _dtw_fortran(r :: Vector{T}, t :: Vector{T}) where T <: Integer
    _dtw_fortran(convert(Vector{Float64}, r), convert(Vector{Float64}, t))
end

libtwdtw = Dict(
    :Float32 => Dict(
        :Int8 => Libdl.dlsym(lib, :__dtw_MOD_single8),
        :Int16 => Libdl.dlsym(lib, :__dtw_MOD_single16),
    ),
    :Float64 => Dict(
        :Int8 => Libdl.dlsym(lib, :__dtw_MOD_double8),
        :Int16 => Libdl.dlsym(lib, :__dtw_MOD_double16),
    ),
)

for ttype in [:Float32, :Float64]
    for xtype in [:Int8, :Int16]
        eval(quote
            function _dtw_fortran(
                r :: Vector{$xtype}, t :: Vector{$ttype}, factor :: $ttype
            ) :: $xtype
                ccall(
                    $(libtwdtw[ttype][xtype]), $ttype,
                    (
                        Ref{$xtype}, Ref{Int32}, Ref{$ttype}, Ref{Int32}, Ref{$ttype}
                    ),
                    r, length(r), t, length(t), factor
                )
            end
        end)
    end
end

function _dtw_fortran(r, t :: Vector{T}, factor) where T <: Integer
    _dtw_fortran(convert(Vector{Float64}, r), convert(Vector{Float64}, t), factor)
end

_dtw_fortran(r, t; factor) = _dtw_fortran(r, t, factor)

function _tw_fortran(rdate :: Vector{Int64}, tdate :: Vector{Int64}, α, β) :: Matrix{Float64}
    M = length(rdate); N = length(tdate)
    tw = Matrix{Float64}(undef, M, N)
    ccall(
        Libdl.dlsym(lib, :__dtw_core_MOD_tw_double64), Nothing,
        (
            Ref{Int64}, Ref{Int32}, Ref{Int64}, Ref{Int32},
            Ref{Float64}, Ref{Float64}, Ref{Float64}
        ),
        rdate, M, tdate, N, Float64(α), Float64(β), tw
    )
    tw
end

_tw_fortran(rdate :: Vector{Int32}, tdate :: Vector{Int32}, α, β) = _tw_fortran(
    convert(Int64, rdate), convert(Vector{Float64}, tdate), α, β
)
_tw_fortran(rdate, tdate; α=0.1, β=50) = _tw_fortran(rdate, tdate, α, β)

libdtw = Dict(
    :Float32 => Libdl.dlsym(lib, :__tw_dtw_MOD_singlex),
    :Float64 => Libdl.dlsym(lib, :__tw_dtw_MOD_doublex),
)

for xtype = [:Float32, :Float64]
    eval(quote
        function _tw_dtw_fortran(
            r :: Vector{$xtype}, t :: Vector{$xtype}, tw :: Matrix{$xtype}
        ) :: $xtype
            ccall(
                $(libdtw[xtype]), $xtype,
                (Ref{$xtype}, Ref{Int32}, Ref{$xtype}, Ref{Int32}, Ref{$xtype}),
                r, length(r), t, length(t), tw
            )
        end
    end)
end

function _tw_dtw_fortran(r :: Vector{T}, t :: Vector{T}, tw) where T <: Integer
    _tw_dtw_fortran(
        convert(Vector{Float64}, r), convert(Vector{Float64}, t),
        convert(Vector{Float64}, tw)
    )
end

libtwdtw = Dict(
    :Float32 => Dict(
        :Int8 => Libdl.dlsym(lib, :__tw_dtw_MOD_single8),
        :Int16 => Libdl.dlsym(lib, :__tw_dtw_MOD_single16),
    ),
    :Float64 => Dict(
        :Int8 => Libdl.dlsym(lib, :__tw_dtw_MOD_double8),
        :Int16 => Libdl.dlsym(lib, :__tw_dtw_MOD_double16),
    ),
)

for ttype in [:Float32, :Float64]
    for xtype in [:Int8, :Int16]
        eval(quote
            function _tw_dtw_fortran(
                r :: Vector{$xtype}, t :: Vector{$ttype},
                tw :: Matrix{$ttype}, factor :: $ttype
            ) :: $ttype
                ccall(
                    $(libtwdtw[ttype][xtype]), $ttype,
                    (
                        Ref{$xtype}, Ref{Int32},
                        Ref{$ttype}, Ref{Int32},
                        Ref{$ttype}, Ref{$ttype}
                    ),
                    r, length(r), t, length(t), tw, factor
                )
            end
        end)
    end
end

function _tw_dtw_fortran(r, t :: Vector{T}, tw, factor) where T <: Integer
    _tw_dtw_fortran(
        convert(Vector{Float64}, r), convert(Vector{Float64}, t),
        convert(Vector{Float64}, tw), factor
    )
end

_tw_dtw_fortran(r, t, tw; factor) = _tw_dtw_fortran(r, t, tw, factor)

libdtw = Dict(
    :Float32 => Libdl.dlsym(lib, :__tw_m_dtw_MOD_singlex),
    :Float64 => Libdl.dlsym(lib, :__tw_m_dtw_MOD_doublex),
)

for xtype = [:Float32, :Float64]
    eval(quote
        function _tw_m_dtw_fortran(
            r :: Vector{$xtype}, t :: Vector{$xtype}, tw :: Matrix{$xtype}
        ) :: $xtype
            ccall(
                $(libdtw[xtype]), $xtype,
                (Ref{$xtype}, Ref{Int32}, Ref{$xtype}, Ref{Int32}, Ref{$xtype}),
                r, length(r), t, length(t), tw
            )
        end
    end)
end

function _tw_m_dtw_fortran(r :: Vector{T}, t :: Vector{T}, tw) where T <: Integer
    _tw_m_dtw_fortran(
        convert(Vector{Float64}, r), convert(Vector{Float64}, t),
        convert(Vector{Float64}, tw)
    )
end

libtwdtw = Dict(
    :Float32 => Dict(
        :Int8 => Libdl.dlsym(lib, :__tw_m_dtw_MOD_single8),
        :Int16 => Libdl.dlsym(lib, :__tw_m_dtw_MOD_single16),
    ),
    :Float64 => Dict(
        :Int8 => Libdl.dlsym(lib, :__tw_m_dtw_MOD_double8),
        :Int16 => Libdl.dlsym(lib, :__tw_m_dtw_MOD_double16),
    ),
)

for ttype in [:Float32, :Float64]
    for xtype in [:Int8, :Int16]
        eval(quote
            function _tw_m_dtw_fortran(
                r :: Vector{$xtype}, t :: Vector{$ttype},
                tw :: Matrix{$ttype}, factor :: $ttype
            ) :: $ttype
                ccall(
                    $(libtwdtw[ttype][xtype]), $ttype,
                    (
                        Ref{$xtype}, Ref{Int32},
                        Ref{$ttype}, Ref{Int32},
                        Ref{$ttype}, Ref{$ttype}
                    ),
                    r, length(r), t, length(t), tw, factor
                )
            end
        end)
    end
end

function _tw_m_dtw_fortran(r, t :: Vector{T}, tw, factor) where T <: Integer
    _tw_m_dtw_fortran(
        convert(Vector{Float64}, r), convert(Vector{Float64}, t),
        convert(Vector{Float64}, tw), factor
    )
end

_tw_m_dtw_fortran(r, t, tw; factor) = _tw_m_dtw_fortran(r, t, tw, factor)


# 序列长度在 10 左右时：
# Linux 工作站下，Julia 版速度慢于 Fortran。耗时约为 Fortran版的 150%。
# Windows 台式机下，Julia 版速度略快于 Fortran。（2.5 倍）
backend = Sys.iswindows() ? :_julia : :_fortran
for funcname = [:_twdtw, :_dtw, :_tw, :_tw_dtw, :_tw_m_dtw]
    funcname0 =  Symbol(funcname, backend)
    eval(:($funcname = $funcname0))
end

"""
    twdtw(r, t, rdate, tdate; α=0.1, β=50)

Return the distance of two series calculted by TWDTW

The time series should not contain NaNs or Infs.

# Arguments
- `r`: real series
- `t`: typical series
- `rdate`: dates of real series (day)
- `tdate`: dates of typical series (day)
- `α`: how rapid the Logistic time weight were given, default 0.1
- `β`: (day) when the time weight reaches 0.5, default 50 days

# Examples
```julia-repl
julia> twdtw(rand(10), rand(10), collect(1:10), collect(1:10))
2.4634249103374004
```
"""
twdtw= _twdtw

"""

    dtw(r, t, rdate, tdate; α=0.1, β=50)

Return the distance of two series calculted by DTW

The time series should not contain NaNs or Infs.

# Arguments
- `r`: real series
- `t`: typical series
- `rdate`: dates of real series (day), Not needed
- `tdate`: dates of typical series (day), Not needed
- `α`: how rapid the Logistic time weight were given, default 0.1, Not needed
- `β`: (day) when the time weight reaches 0.5, default 50 days, Not needed

# Examples
```julia-repl
julia> dtw(rand(10), rand(10), collect(1:10), collect(1:10))
3.0264608748300112

julia> dtw(rand(10), rand(10))
2.426256327862013
```
"""
dtw(r, t, rdate, tdate; α=0.1, β=50, factor) = _dtw(r, t; factor)
dtw(r, t; α=0.1, β=50, factor) = _dtw(r, t; factor)
dtw(r, t, rdate, tdate; α=0.1, β=50) = _dtw(r, t)
dtw(r, t; α=0.1, β=50) = _dtw(r, t)

"""

    tw(rdate, tdate; α=0.1, β=50)

Return the time weight matrix of two date series

The time series should not contain NaNs or Infs.

# Arguments
- `rdate`: dates of real series (day)
- `tdate`: dates of typical series (day)
- `α`: how rapid the Logistic time weight were given, default 0.1
- `β`: (day) when the time weight reaches 0.5, default 50 days

# Examples
```julia-repl
julia> tw(collect(1:4), collect(2:5))
4×4 Matrix{Float64}:
 0.00739154  0.00816257  0.0090133   0.0099518
 0.00669285  0.00739154  0.00816257  0.0090133
 0.00739154  0.00669285  0.00739154  0.00816257
 0.00816257  0.00739154  0.00669285  0.00739154
```
"""
tw = _tw

"""

    tw_dtw(r, t, tw; factor)

Return the distance of two series calculted by TWDTW

The time series should not contain NaNs or Infs.

# Arguments
- `r`: real series
- `t`: typical series
- `tw`: time weight matrix
- `factor`: factor for real series, optional

# Examples
```julia-repl
julia> tweight = tw(collect(1:10), collect(2:11));

julia> tw_dtw(rand(10), rand(10), tweight)
2.689593066906015

julia> tw_dtw(floor.(Int16, 1000rand(10)), rand(10), tweight, 1000e0)
2.689593066906015
```
"""
tw_dtw = _tw_dtw

"""

    tw_m_dtw(r, t, tw; factor)

Return the distance of two series calculted by TWDTW, 
using multipy instead of plus for time weights

The time series should not contain NaNs or Infs.

# Arguments
- `r`: real series
- `t`: typical series
- `tw`: time weight matrix
- `factor`: factor for real series, optional

# Examples
```julia-repl
julia> tweight = tw(collect(1:10), collect(2:11));

julia> tw_m_dtw(rand(10), rand(10), tweight)
0.01988789922291201

julia> tw_m_dtw(floor.(Int16, 1000rand(10)), rand(10), tweight, 1000e0)
0.019811554880000866
```
"""
tw_m_dtw = _tw_m_dtw

end

using .TWDTW
