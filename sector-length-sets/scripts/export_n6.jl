#!/usr/bin/env julia
# Export existing proved geometry and explicit pure-state witnesses only.
# Does not render or modify any existing manuscript/figure/library source.
using LinearAlgebra, Random, Printf
const ATLAS_DIR=normpath(joinpath(@__DIR__,".."))
const PROJECT_DIR=normpath(joinpath(ATLAS_DIR,"..",".."))
include(joinpath(PROJECT_DIR,"Julia","polytope_symmetric_N6_uniform_attainment.jl"))

function atlas_json(io,x)
    if x===nothing; print(io,"null")
    elseif x isa Bool; print(io,x ? "true" : "false")
    elseif x isa AbstractString || x isa Symbol
        s=replace(string(x),"\\"=>"\\\\","\""=>"\\\"","\n"=>"\\n","\r"=>"\\r","\t"=>"\\t")
        print(io,'"',s,'"')
    elseif x isa Real
        isfinite(x) || error("Nonfinite JSON value"); print(io,x)
    elseif x isa AbstractDict || x isa NamedTuple
        print(io,'{');first=true
        for (k,v) in pairs(x)
            first || print(io,',');first=false;atlas_json(io,string(k));print(io,':');atlas_json(io,v)
        end; print(io,'}')
    else
        print(io,'[');first=true
        for v in x;first || print(io,',');first=false;atlas_json(io,v);end
        print(io,']')
    end
end

"""All sectors from independent Dicke-Schmidt marginal purities."""
function atlas_sectors(psi)
    n=length(psi)-1;s=[1.0]
    for t in 1:n
        schmidt=zeros(ComplexF64,t+1,n-t+1)
        for a in 0:t,r in 0:n-t
            schmidt[a+1,r+1]=psi[a+r+1]*sqrt(binomial(t,a)*binomial(n-t,r)/binomial(n,a+r))
        end
        rho=schmidt*schmidt';purity=real(tr(rho*rho))
        value=binomial(n,t)*(2.0^t*purity-sum(binomial(t,k)/binomial(n,k)*s[k+1] for k in 0:t-1))
        push!(s,abs(value)<1e-11 ? 0.0 : value)
    end
    s[2:end]
end

const RECORDS=Any[];const NAMED=Any[];const SEEN=Set{Tuple}()
function add_state(psi,family;label="",id="",parameters=Dict(),definition="",origin="",deduplicate=true)
    isnothing(psi) && return false
    abs(norm(psi)-1)<2e-8 || error("Unnormalized $family $(norm(psi))")
    psi=ComplexF64.(psi)/norm(psi)
    s=atlas_sectors(psi);p=s[1:3]
    maximum(abs.(p-cn_marginal_sector_point(psi)))<1e-10 || error("Marginal conversion mismatch")
    all(s .>= -1e-8) || error("Negative sector")
    abs(sum(s)-63)<1e-8 || error("Pure-state total mismatch")
    key=Tuple(round.(p;digits=7))
    deduplicate && isempty(label) && key in SEEN && return false
    push!(SEEN,key)
    lower=compatibility_bottom(p[1],max(p[2],5.));upper=compatibility_upper(p[1],max(p[2],5.))
    boundary_residual=min(abs(p[3]-upper),abs(p[3]-lower),abs(p[1]),abs(p[2]-5))
    projection=min(p[1],p[2]-5,15-p[2],(p[2]+3)/3-p[1])
    isboundary=min(p[3]-lower,upper-p[3],projection)>=-1e-7 && boundary_residual<1e-7
    rec=Dict("id"=>(isempty(id) ? "n6-$(length(RECORDS)+1)" : id),"label"=>label,
        "N"=>6,"family"=>family,"point"=>p,"sectors"=>s,
        "coefficients"=>[[abs(real(c))<1e-14 ? 0. : real(c),abs(imag(c))<1e-14 ? 0. : imag(c)] for c in psi],
        "parameters"=>parameters,"definition"=>definition,"proofOrigin"=>origin,
        "status"=>"attained","onOuterBoundary"=>isboundary,"boundaryResidual"=>boundary_residual)
    push!(RECORDS,rec);isempty(label) || push!(NAMED,rec)
    true
end
two(k,l,p)=_cn_two_level(k,l,p)
function safe_constructor(f)
    try f() catch e
        e isa ArgumentError || e isa DomainError || rethrow();nothing
    end
end

println("Exporting named states");flush(stdout)
add_state(two(0,6,1.),"C6";label="Product",id="n6-product",definition="|0⟩⊗6",origin="Spin-coherent product state")
add_state(two(0,6,.5),"C6";label="GHZ₆",id="n6-ghz",definition="(D₆⁽⁰⁾+D₆⁽⁶⁾)/√2",origin="Exact Dicke coefficients")
add_state(two(1,5,1.),"C5";label="W₆",id="n6-w",definition="D₆⁽¹⁾",origin="Exact Dicke state")
add_state(two(1,5,.5),"Octahedral";label="Octa",id="n6-octa",definition="(D₆⁽¹⁾+D₆⁽⁵⁾)/√2",origin="Exact octahedral representative")
for k in 2:3
    c=zeros(ComplexF64,7);c[k+1]=1
    add_state(c,"Dicke";label="D₆⁽$(k)⁾",id="n6-dicke-$k",definition="D₆⁽$(k)⁾",origin="Exact Dicke state")
end
add_state(two(0,4,.25),"C4";label="C₄ (pyramid)",id="n6-c4",definition="(D₆⁽⁰⁾+√3D₆⁽⁴⁾)/2",origin="Complete spectral equality classification")
add_state(two(0,4,3/8),"C4";label="Ω₄",id="n6-omega4",definition="√(3/8)D₆⁽⁰⁾+√(5/8)D₆⁽⁴⁾",origin="Attained C4 arc")
add_state(two(2,5,20/21),"C3";label="ψ₆",id="n6-psi6",definition="(√20D₆⁽²⁾+D₆⁽⁵⁾)/√21",origin="Spectral S4 equality state")
add_state(two(0,5,4/7),"C5";label="χ₆",id="n6-chi6",definition="(2D₆⁽⁰⁾+√3D₆⁽⁵⁾)/√7",origin="Spectral S4 equality state; same sectors as ψ₆, inequivalent orbit")
add_state(two(0,5,.7),"C5";label="χ₇/₁₀ (C₅ junction)",id="n6-chi7-10",definition="√(7/10)D₆⁽⁰⁾+√(3/10)D₆⁽⁵⁾",origin="Exact attained C5 rim")
add_state(two(0,5,.5),"C5";label="χ₁/₂",id="n6-chi1-2",definition="(D₆⁽⁰⁾+D₆⁽⁵⁾)/√2",origin="Exact attained C5 rim")
for (id,label,p) in (("n6-ac2","AC₂",[0.,5.,640/81]),("n6-psi0","ψ₀",[0.,5.1,8.]))
    add_state(c3_state_at_value(p...),"C3";label=label,id=id,definition="C3 three-level representative; exact coefficients shown",origin="Complete C3 spectral equality family")
end
add_state(n6gray_C2_null_witness(0.,F0_NULL_S2_MAX),"C2";label="Ξ*",id="n6-xi",definition="Sharp zero-dipole null-roof endpoint; C2 coefficients shown",origin="Proved sharp fixed S1=0,S3=8 bound")
add_state(n6research_omega_state(),"C2";label="ω",id="n6-omega",definition="Odd C2 radical representative; coefficients shown",origin="Complete rank-four supporting-plane equality classification")
add_state(n6rank4_estar_state(),"C2";label="E*",id="n6-estar",definition="Exceptional even C2 equality orbit; coefficients shown",origin="Complete rank-four supporting-plane equality classification")
add_state(n6pairedrim_endpoint_state(),"C2";label="η₆ (sharp C₂ endpoint)",id="n6-eta6",definition="(5D₀+√35D₂+√15D₄−√21D₆)/√96",origin="Sharp paired-null threshold S1=2/27")

println("Sampling explicit boundary families");flush(stdout)
for p in range(0,1;length=101), (k,l,family) in ((0,6,"C6"),(0,5,"C5"),(0,4,"C4"),(1,5,"Octa–W"),(0,1,"W–product"))
    add_state(two(k,l,p),family;parameters=Dict("p"=>p),definition="√p D₆⁽$(k)⁾+√(1−p) D₆⁽$(l)⁾",origin="Exact two-Dicke-level family")
end
for r in range(0,1/9;length=81)
    p=c3_spectral_point(r)
    add_state(c3_state_at_value(p...),"C3 spectral arc";parameters=Dict("r"=>r),origin="Exact spectral equality family")
end
for t in range(1,N6RANK4_TMAX;length=81)
    add_state(n6rank4_equality_state(t),"Rank-four equality segment";parameters=Dict("t"=>t),origin="Explicit radical family; complete equality classification")
end
for h in range(0,1/3;length=81)
    add_state(n6pairedrim_state(h),"Paired null rim";parameters=Dict("h"=>h),origin="Explicit radical family; sharp threshold 2/27")
end
for f in range(0,.5;length=41)
    bounds=n6gray_C2_null_bounds(f);isnothing(bounds) && continue
    for y in range(bounds[1],bounds[2];length=17)
        add_state(safe_constructor(()->n6gray_C2_null_witness(f,y)),"C2 null-roof band";
            parameters=Dict("f"=>f,"S2"=>y),origin="Exact C2 monotone inverse; attainment only outside certified sharp strips")
    end
end
for r in range(.005,.995;length=25),phi in range(0,pi/2;length=17)
    for (constructor,family) in ((n6op_state,"Octa–product roof"),(n6gp_state,"GHZ–product roof"))
        add_state(constructor(r,phi),family;parameters=Dict("r"=>r,"phi"=>phi),origin="Exact Householder state family")
    end
end
for k in vcat(0.,10.0.^range(-3,1.5;length=29)),c in range(1/3,1;length=18)
    add_state(n6fo_state(k,c),"Filtered octahedral lower plane";
        parameters=Dict("k"=>k,"c"=>c),origin="Complete SL(2,C) filtered-octahedral lower-plane image")
end
for theta in range(0,pi;length=25),phi in range(0,2pi;length=25)
    add_state(n6seven_state(cos(theta),sin(theta)*cos(phi),sin(theta)*sin(phi)),"Seven-level null conic";
        parameters=Dict("theta"=>theta,"phi"=>phi),origin="Exact normalized null-conic pure-state family")
end
for t in range(.05,1;length=23),k in range(-1,1;length=21)
    add_state(n6smallfour_state(t,k),"Four-level null roof";parameters=Dict("t"=>t,"k"=>k),origin="Exact four-level null pure-state family")
end
# Explicit x=0 C3 side states, without filling a convex hull.
for s in range(0,1;length=31),phase in range(0,pi/2;length=21)
    add_state(c3_three_level_state(s,0.,phase),"C3 zero-dipole side";
        parameters=Dict("s"=>s,"phase"=>phase),origin="Exact three-level family")
end
for a in range(-1.7,1.7;length=17),b in range(0,1.7;length=13)
    add_state(n6ww_state(a,b),"Rotated W–W roof";parameters=Dict("a"=>a,"b"=>b),origin="Exact superposition of two rotated W states")
end

# Targeted witnesses on the complete side and paired-purity faces.
for x in range(0,8/3;length=31)
    for z in range(5x,compatibility_upper(x,5.);length=17)
        add_state(safe_constructor(()->cn_c2_even_y5_witness(x,z)),"C2 quadrupole-free side";
            parameters=Dict("S1"=>x,"S3"=>z),origin="Exact two-triangle C2 reconstruction on S2=5")
    end
end
for x in range(0,3/8;length=31)
    qlo,qhi=n6research_paired_q_bounds(x)
    for q in range(qlo,qhi;length=19)
        y=7+3x+8q;z=8+2x-8q
        add_state(safe_constructor(()->cn_boundary_witness([x,y,z])),"C2 paired-purity face";
            parameters=Dict("S1"=>x,"q"=>q),origin="Explicit C2 reconstruction; complete sharp paired-face image")
    end
end
for r in range(0,1;length=35)
    lo,hi=n6gray_F0_attained_bounds(r)
    for y in range(lo,hi;length=17)
        certificate=safe_constructor(()->n6gray_F0_certificate(r,y))
        isnothing(certificate) && continue
        add_state(certificate.state,"Zero-dipole attained side";
            parameters=Dict("r"=>r,"S2"=>y),origin="Audited continuous C2/C3/halfturn state constructions")
    end
end

println("Sampling current exact outer graphs");flush(stdout)
A,b=violet_view_linear_constraints()
patches=exact_curved_region_meshes(A,b;u_samples=17,y_samples=29)
exported=Any[];edges=Any[]
function export_patch(kind,points,faces,surface)
    flags=Bool[]
    for face in faces
        p=sum(points[i] for i in face)/3
        push!(flags,kind==:fourcopy_lower ? false : uniform_cn_mask(p,surface))
    end
    push!(exported,Dict("kind"=>string(kind),"surface"=>string(surface),"positions"=>points,
        "triangles"=>[[i-1 for i in f] for f in faces],"attained"=>flags))
end
for patch in patches
    surface=uniform_cn_surface(patch)
    # Rebuild two side faces with vertical resolution to display attained patches.
    if surface in (:x0,:y5)
        pts=Vector{Float64}[];fs=Vector{Int}[];steps=40;levels=20
        for a in range(0,1;length=steps+1)
            x,y=surface==:x0 ? (0.,5+10a) : (8a/3,5.)
            lo=compatibility_bottom(x,y);hi=compatibility_upper(x,y)
            for v in range(0,1;length=levels+1);push!(pts,[x,y,lo+v*(hi-lo)]);end
        end
        for j in 0:steps-1,i in 1:levels
            aa=j*(levels+1)+i;bb=aa+1;cc=aa+levels+1;dd=cc+1
            push!(fs,[aa,bb,cc],[bb,dd,cc])
        end
        export_patch(patch.kind,pts,fs,surface)
    else
        coords=GeometryBasics.coordinates(patch.mesh)
        used=sort!(unique([convert(Int,i) for f in GeometryBasics.faces(patch.mesh) for i in f]));ids=Dict(v=>i for (i,v) in enumerate(used))
        pts=[Float64.(collect(coords[i])) for i in used]
        fs=[[ids[convert(Int,i)] for i in f] for f in GeometryBasics.faces(patch.mesh)]
        export_patch(patch.kind,pts,fs,surface)
    end
    for rim in patch.rims
        length(rim)>1 || continue
        for part in uniform_cn_rim_parts(rim)
            pp=part.points;indices=sort!(unique(vcat(collect(1:5:length(pp)),length(pp))))
            push!(edges,Dict("points"=>[pp[i] for i in indices],"attained"=>part.attained,"kind"=>string(patch.kind)))
        end
    end
    println("Exported patch ",patch.kind," / ",surface);flush(stdout)
end
for path in uniform_body_edge_paths()
    pp=path.points;idx=sort!(unique(vcat(collect(1:8:length(pp)),length(pp))))
    push!(edges,Dict("kind"=>string(path.name),"attained"=>true,"points"=>[pp[i] for i in idx]))
end
knownpoints=[s["point"] for s in NAMED]
graynodes=[Dict("point"=>Float64.(p),"status"=>"outer-only","label"=>"Outer-region junction") for p in compatibility_isolated_corners()
    if all(norm(p-q)>1e-7 for q in knownpoints)]
data=Dict("N"=>6,"axes"=>["S1","S2","S3"],"dimension"=>3,
    "description"=>"Current proved necessary outer region, explicit attained boundary patches, and actual pure-state witnesses. No convex hull is used as an achievability claim.",
    "normalization"=>"S_k is the squared Pauli-correlation norm summed over all weight-k Pauli strings; pure-state total sum S_k=63 (k=1,…,6).",
    "outer"=>Dict("patches"=>exported,"edges"=>edges,"nodes"=>graynodes,
        "proofOrigin"=>"N6_proven_outer_region_20261015.tex, current 64/81 four-copy refinement",
        "meshNote"=>"Display triangulation samples exact analytical graphs; triangle colors are conservative center classifications, not a claim that their chords are pure-state images."),
    "states"=>NAMED,"witnesses"=>RECORDS)
mkpath(joinpath(ATLAS_DIR,"data"))
jsonpath=joinpath(ATLAS_DIR,"data","n6.json")
open(jsonpath,"w") do io;atlas_json(io,data);end
open(joinpath(ATLAS_DIR,"data","n6.js"),"w") do io
    print(io,"globalThis.SECTOR_N6_DATA = ");atlas_json(io,data);println(io,";")
end
println("N6 atlas exported: ",length(NAMED)," named states; ",length(RECORDS)," total witnesses; ",
    sum(length(p["triangles"]) for p in exported)," triangles; ",length(edges)," edge runs.")
