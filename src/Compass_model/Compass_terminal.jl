#using Plots
using CairoMakie
using RollingFunctions
using ProgressMeter
using LinearAlgebra
using DelaunayTriangulation
using DelimitedFiles
include("Compass.jl")

px = []
py = []
pz = []

for z in range(0,1,step=0.05)
    for x in range(0,1-z,step=0.05)
        y = round(1 - x - z, digits=2)
        push!(px, x)
        push!(py, y)
        push!(pz, z)
    end
end

p = collect(zip(px,py,pz))

Sx = []
Sy = []
Tx = []
Ty = []

@showprogress for i in 1:length(p)
    sx, sy, tx, ty = main(p[i], 16)
    push!(Sx, sx)
    push!(Sy, sy)
    push!(Tx, tx)
    push!(Ty, ty)
end

open("data/data_L16_2_SxSyTxTy.txt", "w") do io
    writedlm(io, [Sx Sy Tx Ty])
end