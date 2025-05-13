#=
*
*      MSIAC: a Julia module designed to post-process data with SIAC filters
*                                Version 1.0
*
*      Copyright 2022, Licensed under The GNU Lesser General Public License, version 2.1
*      See http://www.opensource.org/licenses/lgpl-2.1.php
*      Written by: Julia Docampo Sanchez
*      Email: julia.docampo@bsc.es
*
=#
__precompile__()

module MSIAC
#using Revise
#include VTK functions:
const useVTK = false

using LinearAlgebra, Einsum,DelimitedFiles, Distributions, Printf,CircularArrays,
Jacobi, DistributedArrays, Distributed, Triangle

import BSplines, FastGaussQuadrature, Jacobi, Polynomials

import Gmsh: gmsh
const JP   = Jacobi
const GQ   = FastGaussQuadrature

if useVTK
	using  WriteVTK, VTKDataIO
end

using PyCall, PyPlot
const fTOL = 1.e-11
const dbl  = Float64
const vecD = Vector{Float64}
const vecI = Vector{Int}
const vecU = Union{Vector{UInt8},Vector{UInt64}}
const uint = Union{UInt8,UInt64}
const matD = Matrix{Float64}
const matI = Matrix{Int}


ele_ca(type::Int) = repeat([1:type;],2)

function tensor_index(n1::Int,n2::Int; n3 = 0, reverse = false::Bool)

	if n3 == 0
		x = vcat([fill(j,n2) for j = 1:n1]...)
		y = repeat([1:n2;], n1)
		if reverse
			return y,x
		else return x,y
		end
	else
		z = repeat([1:n3;], n1 * n2)
		y = repeat( vcat([fill(j,n3) for j = 1:n2]...),n1)
		x = vcat([fill(j,n3 * n2) for j = 1:n1]...)
		if reverse return z, y, x
		else       return x, y, z
		end
	end
end

"""
Structure for quadrature
------------------------
nodes :
	quadrature points
weights:

eQT:
	type: Legendre, Lobatto, etc
"""
struct GZW
   type    :: String
   Q       :: Int
   nodes   :: vecD
   weights :: vecD
   function GZW(tQ::String, Q::Int)
	   if tQ == "Legendre" || tQ == "legendre" || tQ == "LEGENDRE"
	       z,w = GQ.gausslegendre(Q)
	   elseif tQ == "Lobatto" || tQ == "lobatto" || tQ == "LOBATTO"
	       z,w = GQ.gausslobatto(Q)
	   elseif tQ == "Radau" || tQ == "radau" || tQ == "RADAU"
		   z,w = GQ.gaussradau(Q)
	   else
	       print(" Quadrature type ", tQ," not supported. Add it ! \n")
	       return nothing
	   end
	   return new(tQ, Q, z, w)
   end
end

const eTagDict = Dict(1 => 2, 2 => 3, 3 => 4, 4 => 4, 5=>8)

struct MESH
	dim        :: Int
	N          :: Int
	structured :: Bool
	type       :: Int
	toler      :: dbl
	FMAP       :: Vector{Vector{Int}}
	ELEMENTS   :: Dict{Int,vecI}
	NODES      :: Dict{Int,Tuple{vecD,vecI}}
	EDGES      :: Dict{vecI,vecI}
	FACETS     :: Dict{vecI,vecI}
	BNODES     :: Dict{Int,vecI}
	BEDGES     :: Dict{vecI,vecI}
	BFACES     :: Dict{vecI,vecI}

	function MESH(dim::Int,N::Int,FMAP::Vector{Vector{Int}},ELEMENTS::Dict{Int,vecI},NODES::Dict{Int,Tuple{vecD,vecI}},
		          EDGES::Dict{vecI,vecI},FACETS::Dict{vecI,vecI}; BNODES=Dict{Int,vecI}()::Dict{Int,vecI},
				  BEDGES=Dict{vecI,vecI}()::Dict{vecI,vecI},BFACES=Dict{vecI,vecI}()::Dict{vecI,vecI},
				  structured=false::Bool, type=4::Int,toler=1.e-14::dbl)
		 return new(dim,N,structured,type,toler,FMAP,ELEMENTS,NODES,EDGES,FACETS,BNODES, BEDGES, BFACES)
	end

end

"""
holds footprint information: element reference coordinates (z)
kernel coordinates (t) and sorts by elements (e)
"""
mutable struct BP
	z :: vecD
	t :: vecD
	e :: vecI
	v :: Int
	f :: vecI
	BP(z::vecD,t::vecD,e::vecI,v :: Int, f =[-1,-1]::vecI) = new(deepcopy(z), deepcopy(t),deepcopy(e), v, sort(f))
end

function find_data(name::String, data::Matrix{Any})
	i1 = findfirst(x->x == name, data)[1]
	i2 = findfirst(x->x == "\$End"*strip(name,['\$']), data)[1]
	if i1 == nothing || i2 == nothing
		@error " !!!!! MISSING KEYWORD $name IN DATA MATRIX \n"
		return nothing
	end
	info     = data[i1 + 1,:]
	dType    = zeros(Int,info[1])
	datatype = typeof(data[i1+2,3])
	dVals    = Vector{datatype}[[] for m = 1:info[1]]
	for j = i1 + 2:i2-1
		j_st = findfirst(x->x == "",data[j,:])
		dl   =  j_st == nothing ? length(data[j,:]) : j_st - 1 #- length(x)
		dVals[data[j,1]] = datatype.(data[j,3:dl])  # data not necessarily in order
		dType[data[j,1]] = data[j,2]     # data not necessarily in order
	end
	info, dVals, dType, i2
end

include("kernel.jl")
include("mesh.jl")
include("geometry.jl")
include("fields2D.jl")
include("fields3D.jl")
include("convolution2D.jl")
include("convolution3D.jl")





const EXPANSION = Union{EXP_2D, EXP_3D}

mutable struct DATAFIELD
	fields :: Union{Vector{FIELD_2D},Vector{FIELD_3D}}
	fdict  :: Dict{String,Int}
	mesh   :: MESH
	function DATAFIELD(f::Union{Vector{FIELD_2D},Vector{FIELD_3D}}, m::MESH)
		fdict = Dict("" => 1)
		[fdict[f[k].name] = k for k = 1: length(f)]
		return new(f,fdict,m)
	end
end



polynomial_degree(data::DATAFIELD) = data.fields[data.fdict[""]].basis.degree

"""
Modal reconstruction based on expansion (basis) and f sampled at zw12(3) points
Parameters:
===========
basis: a EXPANSION structure (basis, mass matrix,etc)
f: the projection function at a particular element size length(zw1) * length(zw2)(*length(zw3))
zw1: points sampled along 1st dir
zw2: points sampled along 2nd dir
(zw3: points sampled along 3rd dir)

Returns:
========
modes matrix with dimension [data.degree + 1, #elements]
"""
modal_values(basis::EXPANSION, f::matD, zw::Vector{GZW};
			 reverse=false::Bool) = typeof(basis) == EXP_2D ? modal_2Dvalues(basis,f,zw,reverse=reverse) :
                                                              modal_3Dvalues(basis,f,zw,reverse=reverse)

phys_values(basis::EXPANSION, modes::matD, zw::Vector{GZW};
			reverse=false::Bool) = typeof(basis) == EXP_2D ? phys_2Dvalues(basis,modes,zw, reverse=reverse) :
					                                         phys_3Dvalues(basis,modes,zw, reverse=reverse)

get_exact_data(data::DATAFIELD, type::String, zw::Vector{Vector{dbl}};
			   time = 0.0::dbl) = data.mesh.dim == 2 ? get_exact_2Ddata(data.mesh, type, zw,time = time) :
				                                       get_exact_3Ddata(data.mesh, type, zw,time = time)


function create_analytic_field(q::Vector{Int}, eQT::Vector{String}, mesh::Union{MESH,String}, f_name::String;
						       type = "sine"::String, time = 0.0::dbl)

	zw = [GZW(eQT[k],q[k]) for k = 1 : length(eQT)]

	typeof(mesh) == String && ( mesh = load_mesh(mesh) )
	@info " CREATE FIELD $q"
    open(f_name, "w") do io
		println(io,"\$Fields")
		if mesh.dim  == 2
			println(io,mesh.N," ",q[1]," ",eQT[1]," ",q[2]," ",eQT[2])
		else
			println(io,mesh.N," ",q[1]," ",eQT[1]," ",q[2]," ",eQT[2]," ",q[3]," ",eQT[3])
		end
		for j = 1:mesh.N
			j_nodes = mesh.ELEMENTS[j]
			print(io,j," ",length(j_nodes)," ")
			verts   = element_coords(mesh,j)
			t,dumb  = local_region(getfield.(zw,:nodes), verts)
			f = mesh.dim == 3 ? eval_3Danalytic_field(t,type,time=time) :
			                    eval_2Danalytic_field(t,type,time=time)
			[ print(io,dbl(f[l])," ") for l = 1:length(f)]
			print(io,"\n")
		end
		println(io,"\$EndFields")
		end
	return f_name
end
## Geometry functions

"""
Given a mesh & field file (ONE file if .vtk .vtu), creates a DATAFIELD structure storing the mesh and reconstructing a modal approximation from field points

Parameters:
===========
mesh file
field file
** one single .vtu (.vtk) file

Optional:
=========
structured: if mesh is structured, indicate true (faster computation)
modalExp:  expansion type (default: legendre for quads and hierarchy for tris (modified legendre))
degree: polynomial degree for the expansion
reverse: by default (false), y runs first and x second. Reverse = true to switch it
Returns: a DATAFIELD  structure
"""
function load_data(f1::String, f2=nothing::Union{String,Nothing}; structured=false,
	               degree = -1, modalExp = "legendre",reverse= false, structrued = false)
	ext = split(f1, ".")[end]
	if ext == "vtu" || ext == "vtk"
		return load_vtu(f1)
	end
	mesh  = load_mesh(f1; structured = structured)
	field = mesh.dim == 2  ? load_2Dfield(f2, degree = degree, modalExp = modalExp, reverse=reverse ) :
	                         load_3Dfield(f2, degree = degree, modalExp = modalExp, reverse=reverse )
	return DATAFIELD([field], mesh)
end


"""
Updates the field(s) in data evaluated at the quadrature points (zw)
Parameters:
===========
data: a DATAFIELD / FIELD structure
zw: a 2D vector of GZW points

Optional:
=========
idx: if "all" updates all the fields in data, otherwise specify a particular field (p,u,v,etc)
"""
function update_phys!(data::DATAFIELD,q::Vector{String}, n::vecI;idx = "all"::String)
	piv = idx != "all" ? [idx] : getfield.(data.fields,:name)
	for j in piv
		data.mesh.dim == 2 ? update_2Dphys!(data.fields[data.fdict[j]],q,n) :
		                     update_3Dphys!(data.fields[data.fdict[j]],q,n)
	end
end


"""
Given a field name (u,v,w etc) returns the data stored for that variable
"""
get_field(data::DATAFIELD, idx=""::String) = data.fields[data.fdict[idx]].f

function add_break!(breaks::Vector{BP}, pt::BP)
	#@info " +++++++++++ ADDING BREAK $pt"
	isempty(breaks) && (breaks = Vector{BP}([BP(pt.z,pt.t,pt.e,pt.v)]))
	piv = findfirst(x -> norm(x .- pt.z) < fTOL, getfield.(breaks,:z) )
	if piv == nothing
		push!(breaks, BP(pt.z,pt.t,pt.e,pt.v))
	else
		breaks[piv].e = unique!(vcat(breaks[piv].e,pt.e))
	end
	return breaks
end

function wall_distance!(data::MESH, breaks::Vector{BP}, pC::BP, cN::vecD,tDir::vecD)
	verts = element_coords(data,pC.e[1])
	cPN   = (pC.z,cN)
	qc    = ele_ca(length(verts))
	for f in data.FMAP
		face = data.ELEMENTS[pC.e[1]][f]
		length(data.FACETS[sort(face)]) > 1 && continue
		if data.dim == 2
			edge = (coords(data,face[1]),coords(data,face[2]))
			abs(ndot(edge[2] .- edge[1], tDir) ) > 0.9 && continue
			res = line_intersection(cPN,edge,data.toler)
		else
			res = plane_intersection([coords(data,v) for v in face],cPN,data.toler)
		end
		isempty(res) && continue
		v =  data.dim == 2 ? is_vertex(res, [coords(data,j) for j in face], data.toler) : data.dim == 3 && res[1] == 2 ? res[2] : nothing
		v != nothing && ( pC.v = v)
		pC.z .= data.dim == 2 ? res : res[3]
		tACC  = norm(tDir) < data.toler ? 0.0 : norm(cPN[1] .- pC.z) ./ norm(tDir)
		pC.t += norm(cPN[1] .- pC.z) .* tDir
		return tACC, add_break!(breaks,pC)
	end
	return 0.0, breaks
end

function pull_vertex!(data::MESH, pC::BP, dir::vecD)

	# get all the vertex ids. It can have more than one if point is periodic
	vper  = haskey(data.BNODES,pC.v) ? vcat(pC.v,data.BNODES[pC.v]) : [pC.v]
	veles = vcat( [data.NODES[k][2] for k in vper]...)
	vmain = Int[]

	[ (pC.v in data.ELEMENTS[e]) && push!(vmain,e) for e in veles]
	adjs = filter(x -> x != -1, [adjacent(data,pC.e[1],j) for j in data.FMAP])
	aux  = filter(x-> !isempty(intersect(data.ELEMENTS[pC.e[1]],data.ELEMENTS[x])), adjs)
	o1   = filter(x->x in aux, vmain)
	o2   = setdiff(adjs,aux)
	eles = unique(vcat(pC.e[1],o1,vmain,o2,veles))
	# We are going to pull from vertex in eles list so that periodic elements are consider last
	for j in eles
		aux = findfirst(x->x ∈ vper, data.ELEMENTS[j])
		aux == nothing && continue

		vs  = element_coords(data,j)
		cPN = (vs[aux],vs[aux].+ dir)
		for k in data.FMAP
			foundV = false
			if any( k .== aux )
				vl = filter(x -> norm(x - vs[aux]) > data.toler, vs[k])
				for j in vl
					in_segment(j, cPN, data.toler) && (foundV = true ; break)
				end
			elseif (data.dim == 2 && !isempty(line_intersection(cPN,(vs[k[1]],vs[k[2]]),data.toler))) ||
				   (data.dim == 3 && !isempty(plane_intersection(vs[k],cPN, data.toler)))
				   	foundV = true
			end
			!foundV && continue
			pC.v = data.ELEMENTS[j][aux] ; pC.e[1] = j ; pC.z .= vs[aux]
			return
		end
	end
	pC.e[1] = -1
end

function periodic_coords_through_edge(data::MESH, pC::BP)
	per = data.BEDGES[pC.f]
	adj = data.EDGES[per][1]

	aux = findfirst( x-> sort(data.ELEMENTS[pC.e[1]][x]) == pC.f, data.FMAP)
	e1  = data.ELEMENTS[pC.e[1]][data.FMAP[aux]]
	v1  = [ coords(data, j) for j in e1 ]
	t   = global_coords(pC.z,v1)

	i1  = filter( x -> x in data.ELEMENTS[adj], data.BNODES[e1[1]])[1]
	i2  = filter( x -> x in data.ELEMENTS[adj], data.BNODES[e1[2]])[1]

	v2  = [ coords(data, i1), coords(data, i2) ]
	return BP(v2[1] + t .* (v2[2] - v2[1]), pC.t, [adj], 0, per )
end

structCopy(x::T) where T = T([deepcopy(getfield(x, k)) for k ∈ fieldnames(T)]...)

function pull_edge!(data::MESH, pC::BP; dir=dbl[]::vecD, plotFig=false::Bool)

	if isempty(dir)
		# we come from an intersection, just pick the right adjoint
		if data.dim == 2
			haskey(data.BEDGES,pC.f) ? (pC = periodic_coords_through_edge(data,pC)) :
			                           (pC.e[1] = data.EDGES[pC.f][1])
		end
		return pC
	end

	periodic_edge = haskey(data.BEDGES,pC.f)

	fper  = periodic_edge ? [pC.f,data.BEDGES[pC.f]] : [pC.f]
	feles = unique(vcat( pC.e[1], vcat([data.EDGES[k] for k in fper]...)))
	if !in_segment(pC.z,[coords(data,pC.f[1]),coords(data,pC.f[2])], data.toler)
		throw(@info " We don't need to pull anything. Point is not in edge ")
		return pC
	end

	aux = is_vertex(pC.z, [coords(data,pC.f[1]),coords(data,pC.f[2])], data.toler)

	if aux != nothing
		throw(@error " this is a vertex, better to pull $(pC.v) from coordinates ")
		return pC
	end

	if plotFig
		pygui(true)
		fig = plt.figure(10)
		ax  = data.dim == 3 ? fig.add_subplot(projection="3d") : fig.add_subplot(111)
		vs  = hcat(element_coords(data,pC.e[1])...)
		for k in data.FMAP
			fl   = k[vcat( [1:length(k);],1) ]
			ndir = norm(dir)

			if data.dim == 3
				ax.plot(vs[1,fl], vs[2,fl], vs[3,fl])
				ax.scatter(pC.z[1],pC.z[2],pC.z[3], color = "red")
				ax.plot([pC.z[1],pC.z[1] + dir[1]/ndir],[pC.z[2],pC.z[2]+dir[2]/ndir],[pC.z[3],pC.z[3]+dir[3]/ndir], color = "green")
			else
				ax.plot(vs[1,fl], vs[2,fl])
				ax.scatter(pC.z[1],pC.z[2], color = "red")
				ax.plot([pC.z[1],pC.z[1] + dir[1]/ndir],[pC.z[2],pC.z[2]+dir[2]/ndir], color = "green")
			end

		end
		for a in feles
			a == -1 && continue
			va = hcat(element_coords(data,a)...)
			nv = data.ELEMENTS[a]
			for l in data.FMAP
				fl = l[vcat([1:length(l);],1)]
				if data.dim == 3
					ax.plot(va[1,fl], va[2,fl], va[3,fl])
					[ax.text(va[1,fl[i]],va[2,fl[i]],va[3,fl[i]],string(nv[fl[i]])) for i = 1 : length(l) ]
				else
					ax.plot(va[1,fl], va[2,fl])
					[ax.text(va[1,fl[i]],va[2,fl[i]],string(nv[fl[i]])) for i = 1 : length(l) ]
				end
			end
			ca = element_centre(data,a)
			data.dim == 3 ? ax.text(ca[1],ca[2],ca[3], string(a), color="indigo") :
			                ax.text(ca[1],ca[2], string(a), color="indigo")
		end
	end

	back = structCopy(pC)

	stopSearch = -1
	for j in feles
		pC = periodic_edge && j ∉ data.EDGES[back.f] ? periodic_coords_through_edge(data,back) : pC = structCopy(back)
		cPN = (pC.z,pC.z.+ dir)

		vs  = element_coords(data,j)
		ne  = data.ELEMENTS[j]
		for k in data.FMAP
			if data.dim == 2
				any( x->x == sort([ne[k[1]], ne[k[2]]]), fper ) && continue
				!isempty(line_intersection(cPN, (vs[k[1]], vs[k[2]]), data.toler)) && ( stopSearch = j )
			else
				if issubset(pC.f,ne[k])
					!in_plane(vs[k], cPN, data.toler) && continue
					for i = 1:length(k)
						edge = [ne[k[i]], ne[k[mod1(i+1,length(ne))]] ]
						sort(edge) == sort(pC.f) && continue
						res = line_intersection(cPN, (vs[edge[1]], vs[edge[2]]), data.toler)
						if !isempty(res) && res != cPN[1]
							stopSearch = j
							break
						end
					end
				else
					res = plane_intersection(vs[k],cPN, data.toler)
					if !isempty(res) && res[3] != cPN[1]
						stopSearch = j
						break
					end
				end
			end
			stopSearch != -1 && break
		end
		stopSearch != -1 && break
	end
	pC.e[1] = stopSearch
	return pC
end


function pull_face!(data::MESH, pC::BP, dir::vecD, oriFace::vecI; plotFig=false::Bool)

	if in_face([coords(data,j) for j in oriFace], pC.z, data.toler) == [-1,0]
		throw(@info " We don't need to pull anything. Point is not in face ")
		return pC
	end

	face = sort(oriFace)
	periodic_face = haskey(data.BFACES,face)

	fper  = periodic_face ? [face, data.BFACES[face]] : [face]
	feles = unique(vcat( pC.e[1], vcat([data.FACETS[k] for k in fper]...)))


	if plotFig
		pygui(true)
		fig = plt.figure(10)
		ax  = data.dim == 3 ? fig.add_subplot(projection="3d") : fig.add_subplot(111)
		vs  = hcat(element_coords(data,pC.e[1])...)
		for k in data.FMAP
			fl   = k[vcat( [1:length(k);],1) ]
			ndir = norm(dir)

			if data.dim == 3
				ax.plot(vs[1,fl], vs[2,fl], vs[3,fl])
				ax.scatter(pC.z[1],pC.z[2],pC.z[3], color = "red")
				ax.plot([pC.z[1],pC.z[1] + dir[1]/ndir],[pC.z[2],pC.z[2]+dir[2]/ndir],[pC.z[3],pC.z[3]+dir[3]/ndir], color = "green")
			else
				ax.plot(vs[1,fl], vs[2,fl])
				ax.scatter(pC.z[1],pC.z[2], color = "red")
				ax.plot([pC.z[1],pC.z[1] + dir[1]/ndir],[pC.z[2],pC.z[2]+dir[2]/ndir], color = "green")
			end

		end
		for a in feles
			a == -1 && continue
			va = hcat(element_coords(data,a)...)
			nv = data.ELEMENTS[a]
			for l in data.FMAP
				fl = l[vcat([1:length(l);],1)]
				if data.dim == 3
					ax.plot(va[1,fl], va[2,fl], va[3,fl])
					[ax.text(va[1,fl[i]],va[2,fl[i]],va[3,fl[i]],string(nv[fl[i]])) for i = 1 : length(l) ]
				else
					ax.plot(va[1,fl], va[2,fl])
					[ax.text(va[1,fl[i]],va[2,fl[i]],string(nv[fl[i]])) for i = 1 : length(l) ]
				end
			end
			ca = element_centre(data,a)
			data.dim == 3 ? ax.text(ca[1],ca[2],ca[3], string(a), color="indigo") :
			                ax.text(ca[1],ca[2], string(a), color="indigo")
		end
	end

	back = structCopy(pC)

	stopSearch = -1
	for j in feles
		pC  = periodic_face && j ∉ data.FACETS[face] ? periodic_coords_through_face(data,back) : pC = structCopy(back)
		cPN = (pC.z,pC.z.+ dir)

		vs  = element_coords(data,j)
		ne  = data.ELEMENTS[j]
		for k in data.FMAP
			if sort(face) == sort(ne[k])
				!in_plane(vs[k], cPN, data.toler) && continue
				stopSearch = j
				break
			else
				res = plane_intersection(vs[k],cPN, data.toler)
				if !isempty(res) && res[3] != cPN[1]
					stopSearch = j
					break
				end
			end
			stopSearch != -1 && break
		end
		stopSearch != -1 && break
	end
	pC.e[1] = stopSearch
	return pC
end


"""
Given a mesh & field data structure, filters the entire field
Parameters:
===========
data: mesh & field data structrue (FIELD)
kType: kernel type (line / tensor)

Optional:
=========
rx: number of splines in the 1st direction / line direction

lx: order of splines in the 1st direction / line direction

ry: number of splines in the 2nd direction (only useful for 2D tensor filter)

ly: order of splines in the 2nd direction (only useful for 2D tensor filter)

dumpFile: if true, writes an output file with the filtered data
"""
function filter_data(msh::MESH, fld::Union{FIELD_2D,FIELD_3D}, kType::String; rx =(-1)::Int, lx=(-1)::Int, ry=(-1)::Int, ly=(-1)::Int,
	                 dumpFile="nofile"::String, theta="auto"::Union{String,dbl,Vector{dbl}},
					 scaling="auto"::Union{String,dbl,Vector{dbl}},inexact=false::Bool, parallel=true::Bool, track_footprint=false::Bool)

	(rx == -1) && (rx = 2 * fld.basis.degree + 1)
	(ry == -1) && (ry = 2 * fld.basis.degree + 1)
	(lx == -1) && (lx =     fld.basis.degree + 1)
	(ly == -1) && (ly =     fld.basis.degree + 1)

	kx = [rx,lx,rx,lx,1] ; ky =[ry,ly,ry,ly,1]

	idx   = tensor_index(fld.zw[1].Q,fld.zw[2].Q,n3 = msh.dim == 2 ? 0 : fld.zw[3].Q, reverse=fld.reverse)
	Qt    = length(idx[1])
	filt  = zeros(Qt,msh.N)

	if parallel
		iC  = [ CircularArray(fld.zw[j].nodes[ idx[j] ]) for j = 1:msh.dim]
		eC  = vcat([fill(e,Qt) for e in 1:msh.N]...)
		aux = zeros(Qt*msh.N)
		Threads.@threads for e in 1:msh.N * Qt
			aux[e] = filter_convolution(msh,fld,eC[e], msh.dim == 2 ? [iC[1][e],iC[2][e]] : [iC[1][e],iC[2][e],iC[3][e]],kType,kx = kx, ky = ky,theta = theta, scaling = scaling, inexact=inexact)
		end
		filt = reshape(aux, Qt, msh.N)
	else
		for e in 1:msh.N
			for i in 1:Qt
				i = 1 ; e = 102
				zeta = [ fld.zw[j].nodes[ idx[j][i] ] for j = 1:msh.dim]
				filt[i,e] = filter_convolution(msh,fld,e, zeta, kType,kx =kx, ky = ky,theta = theta, scaling = scaling, inexact=inexact, track_footprint=track_footprint)
				@info " RUNNING NO PARALLLEL $zeta "
				return filt
			end
		end
	end

	dumpFile == "nofile" && (return filt)

	ext = split(dumpFile, ".")[end]
	ext == dumpFile && (dumpFile *=".txt")
	open(dumpFile, "w") do io
		println(io,"\$Fields")
		println(io,msh.N," ", fld.zw[1].Q," ", fld.zw[1].type," ", fld.zw[2].Q," ",fld.zw[2].type)
		for j = 1:msh.N
			print(io,j," ",1," ")
			for k = 1:fld.zw[1].Q * fld.zw[2].Q
				print(io,filt[k,j]," ")
			end
			print(io,"\n")
		end
		println(io,"\$EndFields")
	end
	print(" \n -------- WRITTEN ON FILE ", dumpFile,"\n")
	return filt
end

filter_data(data::DATAFIELD, kType::String; f=""::String,rx =(-1)::Int, lx=(-1)::Int, ry=(-1)::Int, ly=(-1)::Int,
	                 dumpFile="nofile"::String, theta="auto"::Union{String,dbl,Vector{dbl}},
					 scaling="auto"::Union{String,dbl,Vector{dbl}}, parallel=true::Bool, inexact=false::Bool, track_footprint=false::Bool) =
					 filter_data(data.mesh,data.fields[data.fdict[f]],kType, rx = rx, lx = lx, ry = ry, ly=ly,
					 dumpFile = dumpFile, theta=theta, scaling=scaling, parallel=parallel, inexact=inexact, track_footprint=track_footprint)


function parallel_conv!(res, eC, i1, i2, i3, msh, fld,kType, kx,ky,theta,scaling,inexact)
	for e in 1:length(eC)
		res[e] = filter_convolution(msh,fld,eC[e],msh.dim == 2 ? [i1[e],i2[e]] : [i1[e],i2[e],i3[e]],
		                            kType,kx=kx,ky=ky,theta=theta,scaling=scaling,inexact=inexact)
	end
end

function filter_data_hpc(msh::MESH, fld::Union{FIELD_2D, FIELD_3D}, kType::String;
 	                     rx =(-1)::Int, lx=(-1)::Int, ry=(-1)::Int, ly=(-1)::Int,
 	                     dumpFile="nofile"::String, theta="auto"::Union{String,Float64,Vector{Float64}},
						 scaling="auto"::Union{String,Float64,Vector{Float64}},inexact=false::Bool)

 	(rx == -1) && (rx = 2 * fld.basis.degree + 1)
 	(ry == -1) && (ry = 2 * fld.basis.degree + 1)
 	(lx == -1) && (lx =     fld.basis.degree + 1)
 	(ly == -1) && (ly =     fld.basis.degree + 1)

	idx   = tensor_index(fld.zw[1].Q,fld.zw[2].Q,n3 = msh.dim == 2 ? 0 : fld.zw[3].Q, reverse=fld.reverse)
	Qt    = length(idx[1])

 	filt = zeros(Qt,msh.N)
 	i1C  = vcat(repeat(fld.zw[1].nodes[idx[1]], msh.N)...)
 	i2C  = vcat(repeat(fld.zw[2].nodes[idx[2]], msh.N)...)
	i3C  =  msh.dim == 3 ? vcat(repeat(fld.zw[3].nodes[idx[3]], msh.N)...) : Int[]

 	eC   = vcat([fill(e,Qt) for e in 1:msh.N]...)

 	aux   = dzeros((Qt*msh.N,), workers(), nworkers())
 	d_i1C = distribute(i1C, procs = workers())
 	d_i2C = distribute(i2C, procs = workers())
	d_i3C = !(isempty(i3C)) ? distribute(i3C, procs = workers()) : Int[]

 	d_eC  = distribute(eC , procs = workers())
	@sync begin
		@info " Distributed computing using $(workers())"
	 	for p in workers()
 			@spawnat p parallel_conv!(localpart(aux),localpart(d_eC),localpart(d_i1C),localpart(d_i2C),
									  isempty(i3C) ? Int[] : localpart(d_i3C),
 				                      msh,fld,kType,[rx,lx,rx,lx,1],[ry,ly,ry,ly,1],theta,scaling,inexact)
 		end
 	end
 	filt = reshape(convert(Array, aux), Qt, msh.N)
 	if dumpFile != "nofile"
 		ext = split(dumpFile, ".")[end]
 		ext == dumpFile && (dumpFile *=".txt")
 		open(dumpFile, "w") do io
 			println(io,"\$Fields")
 			println(io,msh.N," ", fld.zw[1].Q," ", fld.zw[1].type," ", fld.zw[2].Q," ",fld.zw[2].type)
 			for j = 1:msh.N
 				print(io,j," ",1," ")
 				for k = 1:fld.zw[1].Q * fld.zw[2].Q
 					print(io,filt[k,j]," ")
 				end
 				print(io,"\n")
 			end
 			println(io,"\$EndFields")
 		end
 		print(" \n -------- WRITTEN ON FILE ", dumpFile,"\n")
 	end
 	#showWarnings && @info " Finished postprocessing field with $kType kernel. Total points: $(data.N * data.zw[1].Q * data.zw[2].Q) "
 	#@everywhere GC.gc()
 	return filt
end

filter_data_hpc(data::DATAFIELD, kType::String; f=""::String,rx =(-1)::Int, lx=(-1)::Int, ry=(-1)::Int, ly=(-1)::Int,
               dumpFile="nofile"::String, theta="auto"::Union{String,Float64,Vector{Float64}},
   			   scaling="auto"::Union{String,Float64,Vector{Float64}}, parallel=true::Bool, inexact=false::Bool) =
 			   filter_data_hpc(data.mesh,data.fields[data.fdict[f]],kType, rx = rx, lx = lx, ry = ry, ly=ly,
 			   dumpFile = dumpFile, theta=theta, scaling=scaling,inexact=inexact)



function get_element_list!(ker::Vector{KERNEL}, data::MESH, pC::BP, kdir::vecD;
	                       kT = [eSYM, ePOS]::vecI,secondAxis = false::Bool, skip_knots = false)

   	eqK = length(kT) == 2 && ker[1].rC == ker[2].rC && ker[1].l == ker[2].l
	for t in kT
		update_knot_matrix!(ker[t])
		kL,kR = skip_knots ? ([minimum(ker[t].T)],[maximum(ker[t].T)]) : split_knots(ker[t].T)
		b0 = secondAxis ? add_break!(BP[],pC) : BP[]
		bl1, dL = data.dim == 2 ? support2D_fits!(data,structCopy(pC), kdir, kL,breaks = b0) :
		                          support3D_fits!(data,structCopy(pC), kdir, kL,breaks = b0)

		(t == eSYM && dL > 0.0 && !eqK) && continue
		bline, dR = data.dim == 2 ? support2D_fits!(data,structCopy(pC), kdir, kR, breaks = isempty(bl1) ? b0 : reverse!(bl1)) :
		                            support3D_fits!(data,structCopy(pC), kdir, kR, breaks = isempty(bl1) ? b0 : reverse!(bl1))

		(t == eSYM && dR < 0.0 && !eqK) && continue
		shift = dL > 0.0 ? dL : dR < 0.0 ? dR : 0.0
		if shift != 0.0
			length(kT) == 2 && (t = ePOS)
			update_knot_matrix!(ker[t]; lX = shift)
			kL,kR = skip_knots ? ([minimum(ker[t].T)],[maximum(ker[t].T)]) : split_knots(ker[t].T)
			b0 = secondAxis ? add_break!(BP[],pC) : BP[]
			bl1, dL = data.dim == 2 ? support2D_fits!(data,structCopy(pC), kdir, kL, breaks = b0) :
			                          support3D_fits!(data,structCopy(pC), kdir, kL, breaks = b0)

			bline, dR = data.dim == 2 ? support2D_fits!(data,structCopy(pC), kdir, kR, breaks= isempty(bl1)  ? b0 :  reverse!(bl1)) :
									    support3D_fits!(data,structCopy(pC), kdir, kR, breaks= isempty(bl1)  ? b0 :  reverse!(bl1))
		end
		!iszero([dL,dR]) && break
		return t, bline
	end
	return 0,  BP[]
end



"""
Convolution at a particular point. Returns de filtered point
Parameters
==========
data: mesh & field data structure (FIELD)
eID : local element
zeta: relative position (in reference coordinates) at eID
type: either line or tensor

Optional
========
kx: 1st direction: list of number & order of bsplines for symmetric and shifted (default: 2p+1 splines of order p + 1)
ky: 2nd direction (only for tensor filter): list of number & order of bsplines for symmetric and shifted (default: 2p+1 splines of order p + 1)
theta: rotation in the kernel axis
scaling: kernel scaling (sized, ie, scaling = 1 means size 1, not size 1 * element_size )
"""
function filter_convolution(msh::MESH, fld::Union{FIELD_2D, FIELD_3D}, eID::Int, zeta::vecD,type::String;
	 						kx=zeros(Int,0)::vecI, ky=zeros(Int,0)::vecI, theta = "auto"::Union{String,dbl,vecD},
							scaling = "auto"::Uniont{String,dbl,vecD}, inexact = false::Bool, track_footprint=false::Bool)
	p = fld.basis.degree
	isempty(kx) && (kx = [2+p+1,p+1,2+p+1,p+1,1])
	if msh.dim == 3
		return line_3Dconvolution(msh,fld,varKernel(kx),eID, zeta, theta,scaling, track_footprint)
	end
	if type == "tensor"
		isempty(kx) && (kx = [2*p+1,p+1,2*p+1,p+1,1])
		isempty(ky) && (ky = [2*p+1,p+1,2*p+1,p+1,1])
		return tensor_convolution(msh,fld,[varKernel(kx) varKernel(ky)],eID, zeta, theta, scaling, inexact,track_footprint)
	else
		return line_convolution(msh,fld,varKernel(kx),eID, zeta, theta,scaling, track_footprint)
	end
end



"""
Test function: filters a particular data point
Parameters:
===========
rx: number of splines in the 1st direction / line direction
lx: order of splines in the 1st direction / line direction
ry: number of splines in the 2nd direction (for line filter, write 0)
ly: order of splines in the 2nd direction (for line filter, write 0)
data: mesh & field data structure (twoD_mesh_field)
e: element id
zeta: reference coordinates at element
kType: kernel type (line / tensor)

Optional:
========
theta = kernel orientation (default Cartesian for tensor and diagonal for line)
scaling = kernel scaling (default element size (hx,hy) and element diagonal (for line))
"""
test_filter(rx::Int, lx::Int, ry::Int, ly::Int, msh::MESH,fld::Union{FIELD_2D, FIELD_3D}, e::Int,zeta::vecD,kType::String;
	                 theta = "auto"::Union{vecD,dbl,String}, scaling="auto"::Union{vecD,dbl,String}, track_footprint=false::Bool) =
					 msh.dim == 2 ? filter_convolution(msh,fld, e, zeta , kType, track_footprint=track_footprint,
					                                    kx =[rx,lx,rx,lx,1], ky =[ry,ly,ry,ly,1], theta = theta, scaling=scaling) :
									line_3Dconvolution(msh,fld, varKernel([rx,lx,ry,ly,1]), e, zeta ,theta, scaling, track_footprint)

test_filter(rx::Int, lx::Int, ry::Int, ly::Int, data::DATAFIELD, e::Int,zeta::vecD,kType::String;
	                 theta = "auto"::Union{vecD,dbl,String}, scaling="auto"::Union{vecD,dbl,String}, track_footprint=false::Bool) =
                  filter_convolution(data.mesh,data.fields[data.fdict[""]], e, zeta, kType, track_footprint=track_footprint,
					 kx =[rx,lx,rx,lx,1], ky =[ry,ly,ry,ly,1], theta = theta, scaling=scaling)


plot_field(data::DATAFIELD;
           idx=""::String,f=matD(undef,0,0)::matD, showMesh=false::Bool, block=false) =
		    data.mesh.dim == 2 ? plot_2Dfield(data.mesh, data.fields[data.fdict[idx]],f=f,showMesh = showMesh, block=block) :
			                     plot_3Dfield(data.mesh, data.fields[data.fdict[idx]],f=f,showMesh = showMesh, block=block)

"""
Physical (field) values given an element an its local coordinates (in reference system)
Parameters:
===========
DATAFIELD structure or FIELD structure
e: the element ID
z: points coordinates in reference system (chi: [-1,1]^2)

Returns:
========
Field values at point
"""
function evaluate_point(data::DATAFIELD,e::Int,z::vecD, idx="all"::String)
	piv = idx != "all" ? [idx] : getfield.(data.fields,:name)
	val = zeros(length(piv))
	for j = 1 : length(piv)
		val[data.fdict[piv[j]]] = data.mesh.type <= 4 ? evaluate_2Dpoint(data.fields[j],e,z) :
		                                                evaluate_3Dpoint(data.fields[j],e,z)
	end
	return val
end

get_evaluation_weights(data::DATAFIELD,idx=""::String) = get_evaluation_weights(data.fields[data.fdict[idx]])
get_evaluation_points(data::DATAFIELD, idx=""::String) = get_evaluation_points(data.fields[data.fdict[idx]])

end
