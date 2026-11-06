#!/usr/bin/env julia
using LinearAlgebra
const ATLAS_DIR=normpath(joinpath(@__DIR__,".."))
include(normpath(joinpath(ATLAS_DIR,"..","..","Julia","polytope_symmetric_N6_uniform_attainment.jl")))
rows=[parse.(Float64,split(line,',')) for line in readlines(joinpath(ATLAS_DIR,"data","n6_boundary_inputs.csv"))]
output=zeros(length(rows),6)
for (j,row) in enumerate(rows)
    x,y,z=row
    lower=compatibility_bottom(x,max(y,5.));upper=compatibility_upper(x,max(y,5.))
    projection=min(x,y-5,15-y,(y+3)/3-x)
    residual=min(abs(z-lower),abs(z-upper),abs(x),abs(y-5))
    inside=min(z-lower,upper-z,projection)>=-1e-7
    output[j,:].=[j,residual,z-lower,upper-z,projection,inside && residual<1e-7 ? 1 : 0]
end
open(joinpath(ATLAS_DIR,"data","n6_boundary_audit_rows.csv"),"w") do io
    for row in eachrow(output);println(io,join(row,','));end
end
println("Current outer-graph boundary audit: ",sum(output[:,6]),"/",length(rows)," on boundary; max residual ",maximum(output[:,2]),
    "; min floor gap ",minimum(output[:,3]),"; min roof gap ",minimum(output[:,4]),"; min projection gap ",minimum(output[:,5]))
