#!/usr/bin/env julia
# Exact marker inventory from the manuscript renderer; no figure rendering.
using LinearAlgebra
const ATLAS_DIR=normpath(joinpath(@__DIR__,".."))
include(normpath(joinpath(ATLAS_DIR,"..","..","Julia","polytope_symmetric_N6_uniform_attainment.jl")))
function write_json(io,x)
    if x isa AbstractString || x isa Symbol
        s=replace(string(x),"\\"=>"\\\\","\""=>"\\\"","\n"=>"\\n");print(io,'"',s,'"')
    elseif x isa Bool;print(io,x ? "true" : "false")
    elseif x isa Real;isfinite(x) || error("Nonfinite");print(io,x)
    elseif x isa AbstractDict
        print(io,'{');first=true
        for (k,v) in x;first || print(io,',');first=false;write_json(io,k);print(io,':');write_json(io,v);end
        print(io,'}')
    else
        print(io,'[');first=true
        for v in x;first || print(io,',');first=false;write_json(io,v);end
        print(io,']')
    end
end

# These are verbatim the node inventory expressions in
# polytope_symmetric_N6_uniform_attainment.jl, lines 535–561.
landmarks=vcat(KNOWN_VERTEX_STATES,KNOWN_BOUNDARY_STATES,[ANTICOHERENT_EQUALITY_STATE])
labeled_nodes=vcat([(collect(s.target),s.label) for s in landmarks],
    [(collect(NEIGHBORHOOD_ATTAINED_ENDPOINT),L"|\psi_0\rangle"),
     (collect(ACHIEVABLE_SEGMENT_ENDPOINT),L"|\Xi_*\rangle")],
    uniform_c5_landmarks(),uniform_sharp_c2_landmarks())
attained_points=[p for (p,_) in labeled_nodes]
nodes=Any[]
for point in compatibility_isolated_corners()
    any(norm(point-p)<1e-8 for p in attained_points) && continue
    push!(nodes,Dict("point"=>Float64.(point),"attained"=>false,"status"=>"outer-only",
        "manuscriptLabel"=>"","source"=>"compatibility_isolated_corners()"))
end
for (point,label) in labeled_nodes
    caption=uniform_landmark_caption(point,label)[1]
    push!(nodes,Dict("point"=>Float64.(point),"attained"=>true,"status"=>"attained",
        "manuscriptLabel"=>string(caption),"source"=>"labeled_nodes in manuscript uniform renderer"))
end
length(nodes)==27 || error("Current manuscript marker inventory is not 27")
length(labeled_nodes)==16 || error("Expected 16 violet state markers")
for i in eachindex(nodes),j in 1:i-1
    norm(nodes[i]["point"]-nodes[j]["point"])>1e-8 || error("Duplicate node")
end
callouts=[
    Dict("point"=>c3_spectral_point(1/18),"label"=>"ψᵣ family","manuscriptLabel"=>"|\\psi_r\\rangle","family"=>"C3 spectral arc"),
    Dict("point"=>c4_attained_point(7/16),"label"=>"Arc: Ω_f, β ≤ f ≤ 1/2","manuscriptLabel"=>"\\mathrm{arc:}\\ \\{|\\Omega_f\\rangle\\}_{\\beta\\le f\\le1/2}","family"=>"C4"),
    Dict("point"=>c5_attained_point(.87),"label"=>"Rim: χ_p, 1/2 ≤ p ≤ 1","manuscriptLabel"=>"\\mathrm{rim:}\\ \\{|\\chi_p\\rangle\\}_{1/2\\le p\\le1}","family"=>"C5"),
    Dict("point"=>n6rank4_segment_point(N6RANK4_QMAX/2),"label"=>"Segment: Ψ_t","manuscriptLabel"=>"\\mathrm{segment:}\\ |\\Psi_t\\rangle","family"=>"Rank-four equality segment"),
]
open(joinpath(ATLAS_DIR,"data","n6-manuscript-nodes-base.json"),"w") do io
    write_json(io,Dict("nodes"=>nodes,"callouts"=>callouts))
end
println("Exact current manuscript marker inventory: ",length(nodes)," nodes, ",length(labeled_nodes)," violet, ",length(nodes)-length(labeled_nodes)," gray; four non-node family leader attachments.")
