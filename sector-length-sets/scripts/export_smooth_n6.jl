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

# Keep full double-precision analytical graph coordinates in this isolated export.
function compatibility_mesh(points,faces)
    good=TriangleFace{Int}[]
    for f in faces
        a,b,c=(Float64.(collect(points[i])) for i in f)
        norm(cross(b-a,c-a))>1e-12 || continue
        push!(good,TriangleFace{Int}(Int(f[1]),Int(f[2]),Int(f[3])))
    end
    isempty(good) && return nothing
    used=sort!(unique([Int(i) for f in good for i in f])); mapping=Dict(i=>j for (j,i) in enumerate(used))
    ps=GeometryBasics.Point3{Float64}[GeometryBasics.Point3{Float64}(points[i]...) for i in used]
    fs=[TriangleFace{Int}(mapping[Int(f[1])],mapping[Int(f[2])],mapping[Int(f[3])]) for f in good]
    ns=[zeros(3) for _ in ps]
    for f in fs
        a,b,c=(ps[Int(i)] for i in f)
        n=cross(b-a,c-a)
        for i in f;ns[Int(i)].+=n;end
    end
    normals=GeometryBasics.Vec3{Float64}[
        GeometryBasics.Vec3{Float64}((norm(n)>1e-14 ? n/norm(n) : [0.,0.,1.])...) for n in ns]
    GeometryBasics.Mesh(ps,fs;normal=normals)
end

const STARTED=time()
const EXPORTED=Any[]
const EDGES=Any[]
const STATS=Any[]
const MASK_CACHE=Dict{Tuple{Symbol,Float64,Float64,Float64},Bool}()
function projected(p,surface)
    q=Float64.(p)
    surface==:upper && (q[3]=compatibility_upper(q[1],max(q[2],5.)))
    surface==:lower && (q[3]=compatibility_bottom(q[1],max(q[2],5.)))
    q
end
function flag(p,surface,kind)
    kind in (:fourcopy_lower,:compatibility_lower) && return false
    q=projected(p,surface)
    key=(surface,round(q[1];digits=11),round(q[2];digits=11),round(q[3];digits=11))
    get!(MASK_CACHE,key) do
        uniform_cn_mask(q,surface)
    end
end
normalworld(n)=begin
    v=Float64[-4.5n[1],-11.5n[3],-8.75n[2]]
    norm(v)>1e-14 ? v/norm(v) : Float64[0,1,0]
end

# Refine only constructive color transitions. Split shared edges on both
# incident triangles, so the refinement introduces no cracks or T-junctions.
# Coordinates are always evaluated on the existing analytical graph.
function refine_color_transition!(pts,ns,fs,surface,kind;levels=3)
    kind in (:fourcopy_lower,:compatibility_lower) && return fs
    edge(a,b)=a<b ? (a,b) : (b,a)
    for level in 1:levels
        flags=[flag(p,surface,kind) for p in pts]
        marked=Set{Tuple{Int,Int}}()
        for f in fs
            a,b,c=f
            values=(flags[a],flags[b],flags[c])
            homogeneous=values[1]==values[2]==values[3]
            center=flag((pts[a]+pts[b]+pts[c])/3,surface,kind)
            if !homogeneous || center!=values[1]
                push!(marked,edge(a,b),edge(b,c),edge(c,a))
            end
        end
        isempty(marked) && break
        midpoints=Dict{Tuple{Int,Int},Int}()
        for (a,b) in sort!(collect(marked))
            p=projected((pts[a]+pts[b])/2,surface)
            n=ns[a]+ns[b];n=norm(n)>1e-14 ? n/norm(n) : ns[a]
            push!(pts,p);push!(ns,n);midpoints[(a,b)]=length(pts)
        end
        refined=Vector{Int}[]
        add(a,b,c)=push!(refined,[a,b,c])
        for f in fs
            a,b,c=f
            ab=get(midpoints,edge(a,b),0);bc=get(midpoints,edge(b,c),0);ca=get(midpoints,edge(c,a),0)
            if ab!=0 && bc!=0 && ca!=0
                add(a,ab,ca);add(ab,b,bc);add(ca,bc,c);add(ab,bc,ca)
            elseif ab!=0 && bc!=0
                add(a,ab,c);add(ab,bc,c);add(ab,b,bc)
            elseif ab!=0 && ca!=0
                add(a,ab,ca);add(ab,b,c);add(ab,c,ca)
            elseif bc!=0 && ca!=0
                add(a,b,ca);add(b,bc,ca);add(ca,bc,c)
            elseif ab!=0
                add(a,ab,c);add(ab,b,c)
            elseif bc!=0
                add(a,b,bc);add(a,bc,c)
            elseif ca!=0
                add(a,b,ca);add(ca,b,c)
            else
                push!(refined,f)
            end
        end
        fs=refined
    end
    fs
end

function clipped_patch(patch)
    surface=uniform_cn_surface(patch);kind=patch.kind
    pts=[projected(Float64.(collect(p)),surface) for p in GeometryBasics.coordinates(patch.mesh)]
    ns=[Float64.(collect(n)) for n in patch.mesh.normal]
    fs=[[convert(Int,i) for i in f] for f in GeometryBasics.faces(patch.mesh)]
    fs=refine_color_transition!(pts,ns,fs,surface,kind;levels=3)
    positions=Vector{Float64}[];normals=Vector{Float64}[];triangles=Vector{Int}[];attained=Bool[]
    ids=Dict{Tuple,Int}()
    function vertex(p,n)
        key=Tuple(vcat(p,round.(n;digits=8)))
        get!(ids,key) do
            push!(positions,copy(p));push!(normals,round.(normalworld(n);digits=8))
            length(positions)
        end
    end
    function emit(vs,yes)
        length(vs)>=3 || return
        for j in 2:length(vs)-1
            tri=(vs[1],vs[j],vs[j+1])
            norm(cross(tri[2].p-tri[1].p,tri[3].p-tri[1].p))>1e-13 || continue
            # As in the established renderer, every violet cell is guarded by
            # its exact graph centroid; failed membership means uncertified.
            accepted=yes && flag(sum(v.p for v in tri)/3,surface,kind)
            indices=[vertex(v.p,v.n)-1 for v in tri]
            length(unique(indices))==3 || continue
            push!(triangles,indices);push!(attained,accepted)
        end
    end
    function interpolate(a,b,t)
        p=projected((1-t)*a.p+t*b.p,surface)
        n=(1-t)*a.n+t*b.n;n=norm(n)>1e-14 ? n/norm(n) : a.n
        (p=p,n=n,attained=flag(p,surface,kind))
    end
    function crossing(a,b)
        lo=0.;hi=1.
        for _ in 1:20
            mid=(lo+hi)/2
            flag((1-mid)*a.p+mid*b.p,surface,kind)==a.attained ? (lo=mid) : (hi=mid)
        end
        interpolate(a,b,(lo+hi)/2)
    end
    function split(vs,depth=0)
        flags=[v.attained for v in vs]
        if all(==(flags[1]),flags)
            center_flag=flag(sum(v.p for v in vs)/3,surface,kind)
            if center_flag!=flags[1] && depth<2
                n=sum(v.n for v in vs);n=n/norm(n)
                mid=(p=projected(sum(v.p for v in vs)/3,surface),n=n,attained=center_flag)
                for j in 1:3;split([vs[j],vs[mod1(j+1,3)],mid],depth+1);end
            else
                emit(vs,all(flags)&&center_flag)
            end
            return
        end
        crossings=Dict{Int,Any}()
        for j in 1:3
            k=mod1(j+1,3);flags[j]!=flags[k] && (crossings[j]=crossing(vs[j],vs[k]))
        end
        for yes in (false,true)
            polygon=Any[]
            for j in 1:3
                flags[j]==yes && push!(polygon,vs[j])
                haskey(crossings,j) && push!(polygon,crossings[j])
            end
            emit(polygon,yes)
        end
    end
    # Fast, rigorous exclusions avoid asking a family inverse on known excluded
    # open lower sheets. Endpoints retain their separate manuscript markers.
    if kind in (:fourcopy_lower,:compatibility_lower)
        for f in fs
            vs=[(p=projected(pts[i],surface),n=ns[i],attained=false) for i in f];emit(vs,false)
        end
    else
        vertices=[(p=projected(pts[i],surface),n=ns[i],attained=flag(pts[i],surface,kind)) for i in eachindex(pts)]
        for (i,f) in enumerate(fs)
            split([vertices[j] for j in f])
            i%5000==0 && (println("  ",kind," classified ",i," / ",length(fs));flush(stdout))
        end
    end
    # Rebuild a shared, area-weighted normal field from the FINAL indexed
    # world geometry, including graph-projected clipping vertices. Computing
    # crosses after the axis swap incorporates its negative determinant.
    wp=[[2(p[1]+1)/9-1,2(p[3]+1)/23-1,2(p[2]+1.5)/17.5-1] for p in positions]
    accum=[zeros(3) for _ in positions]
    for f in triangles
        a,b,c=(wp[i+1] for i in f);n=cross(b-a,c-a)
        for i in f;accum[i+1].+=n;end
    end
    for i in eachindex(normals)
        norm(accum[i])>1e-16 && (normals[i]=round.(accum[i]/norm(accum[i]);digits=8))
    end
    push!(STATS,Dict("kind"=>string(kind),"surface"=>string(surface),"triangles"=>length(triangles),"violet"=>count(attained),"positions"=>length(positions)))
    push!(EXPORTED,Dict("kind"=>string(kind),"surface"=>string(surface),"positions"=>positions,
        "triangles"=>triangles,"worldNormals"=>normals,"attained"=>attained))
end

println("Building high-resolution analytical graph meshes");flush(stdout)
A,b=violet_view_linear_constraints()
patches=exact_curved_region_meshes(A,b;u_samples=65,y_samples=101)
println("Clipping constructive color masks to transition curves");flush(stdout)
for patch in patches
    surface=uniform_cn_surface(patch)
    if surface in (:x0,:y5)
        ps=Vector{Float64}[];fs=TriangleFace{Int}[];steps=100;levels=48
        for a in range(0,1;length=steps+1)
            x,y=surface==:x0 ? (0.,5+10a) : (8a/3,5.)
            lo=compatibility_bottom(x,y);hi=compatibility_upper(x,y)
            for v in range(0,1;length=levels+1);push!(ps,[x,y,lo+v*(hi-lo)]);end
        end
        for j in 0:steps-1,i in 1:levels
            aa=j*(levels+1)+i;bb=aa+1;cc=aa+levels+1;dd=cc+1
            push!(fs,TriangleFace{Int}(aa,bb,cc),TriangleFace{Int}(bb,dd,cc))
        end
        mesh=compatibility_mesh(ps,fs);clipped_patch((kind=patch.kind,mesh=mesh,rims=patch.rims))
    else
        clipped_patch(patch)
    end
    println("Finished patch ",patch.kind," / ",surface," at ",round(time()-STARTED;digits=1)," s");flush(stdout)
    for rim in patch.rims
        length(rim)>1 || continue
        for part in uniform_cn_rim_parts(rim)
            push!(EDGES,Dict("points"=>[round.(Float64.(p);digits=10) for p in part.points],"attained"=>part.attained,"kind"=>string(patch.kind)))
        end
    end
end
for path in uniform_body_edge_paths()
    push!(EDGES,Dict("kind"=>string(path.name),"attained"=>true,"points"=>[round.(Float64.(p);digits=10) for p in path.points]))
end
max_graph_error=0.;max_linear_violation=0.;max_side_violation=0.
for patch in EXPORTED
    surface=Symbol(patch["surface"])
    for p in patch["positions"]
        global max_linear_violation=max(max_linear_violation,maximum(A*p-b))
        if surface==:upper
            global max_graph_error=max(max_graph_error,abs(p[3]-compatibility_upper(p[1],max(p[2],5.))))
        elseif surface==:lower
            global max_graph_error=max(max_graph_error,abs(p[3]-compatibility_bottom(p[1],max(p[2],5.))))
        else
            global max_side_violation=max(max_side_violation,compatibility_bottom(p[1],p[2])-p[3],p[3]-compatibility_upper(p[1],p[2]))
        end
    end
end
println("AUDIT max graph error=",max_graph_error," linear violation=",max_linear_violation," side violation=",max_side_violation);flush(stdout)
output=Dict("patches"=>EXPORTED,"edges"=>EDGES,"metadata"=>Dict(
    "origin"=>"Current analytical graph functions, with clipped constructive color transitions",
    "coordinates"=>"Sector coordinates (S1,S2,S3), sampled directly from unchanged bounds",
    "normals"=>"Area-weighted worldNormals computed from final indexed mesh under mapping [2(x+1)/9-1,2(z+1)/23-1,2(y+1.5)/17.5-1], retaining patch creases",
    "triangles"=>sum(s["triangles"] for s in STATS),"patchStats"=>STATS,"elapsedSeconds"=>time()-STARTED,
    "nodeChanges"=>0,"inequalityChanges"=>0,"maxGraphError"=>max_graph_error,"maxLinearViolation"=>max_linear_violation,"maxSideViolation"=>max_side_violation,"colorTransitions"=>"Three conforming local-refinement levels plus twenty bisections at mixed membership edges; no mesh-grid edge is drawn"))
# Publish complete assets atomically so a live browser never sees a partial file.
json_path=joinpath(ATLAS_DIR,"data","n6-smooth-geometry.json")
js_path=joinpath(ATLAS_DIR,"data","n6-smooth-geometry.js")
open(json_path*".tmp","w") do io;atlas_json(io,output);end
open(js_path*".tmp","w") do io
    print(io,"globalThis.N6_SMOOTH_GEOMETRY=");atlas_json(io,output);println(io,';')
end
mv(json_path*".tmp",json_path;force=true)
mv(js_path*".tmp",js_path;force=true)
println("SMOOTH EXPORT COMPLETE: ",sum(s["triangles"] for s in STATS)," triangles in ",time()-STARTED," s");flush(stdout)
