#=
*
*      MSIAC: a Julia module designed to post-process data with SIAC filters
*
*      Copyright 2022, Licensed under The GNU Lesser General Public License, version 2.1
*      See http://www.opensource.org/licenses/lgpl-2.1.php
*      Written by: Julia Docampo Sanchez
*      Email: julia.docampo@bsc.es
*
=#

#using Revise


struct EXP_2D
	type     :: String
	degree   :: Int
	J        :: Vector{Polynomials.Polynomial}
	mass     :: matD
	collapse :: Bool
	m1       :: vecI
	m2       :: vecI
	function EXP_2D(type::String,p::Int, collapse::Bool)
		if type == "Pk" || type == "pk" || type == "PK"
			m1 = vcat([fill(i,p+1-i) for i = 0:p]...)
			m2 = vcat([[0:p-i;] for i = 0:p]...)
			return new(type,p,Vector{Polynomials.Polynomial}(undef,1),ones(1,1),collapse,m1,m2)
		else
			m1,m2 = tensor_index(p+1,p+1)
			if type == "legendre" || type == "LEGENDRE" || type == "Legendre"
				J  = [JP.poly_legendre(i) for i = 0:p]
			else
				if p == 0
					J  = [JP.poly_legendre(0)]
				else
					PJ = [JP.poly_jacobi(i,1.0,1.0) for i = 0:p-2]
					J   = vcat(Polynomials.Polynomial([0.5,.-0.5]),
						  [Polynomials.Polynomial([0.25,0.0,.-0.25]) * PJ[i] for i = 1:p-1],
						       Polynomials.Polynomial([0.5,0.5]))
				end
		    end
			# Get mass matrix
			q1    = GZW("legendre", p + 1)
			q2    = GZW("legendre", p + 1)
			t1,t2 = tensor_index(q1.Q , q2.Q)
			pl2   = zeros((p+1)^2,(p+1)^2)
			n1    = q1.nodes[t1] ;	n2    = q2.nodes[t2]
			we    = q1.weights[t1] .* q2.weights[t2]
			[pl2[i,j] = sum(J[m1[i]].(n1) .* J[m2[i]].(n2) .* J[m1[j]].(n1) .* J[m2[j]].(n2) .* we) for i = 1 : length(m1), j = 1:length(m1) ]
			mij  = inv(pl2)
			return new(type,p,J,mij,collapse,m1,m2)
		end
	end
end


mutable struct FIELD_2D
	zw      :: Vector{GZW}
	reverse :: Bool
	f       :: matD
	modes   :: matD
	name    :: String
	basis   :: EXP_2D
	function FIELD_2D(type::String, p::Int, collapse::Bool; reverse=false::Bool,zw = GZW[]::Vector{GZW},
		              f = zeros(0,0)::matD, m = zeros(0,0)::matD,name = ""::String)

				   		basis = EXP_2D(type,p,collapse)
				   		qt = ["legendre","legendre"] ; Q = p + 2
			   		   #collapse && (qt[2] = "Radau")
			   		   isempty(zw) && (zw = [GZW(qt[1],Q),GZW(qt[2],Q)])
					   ff = !isempty(f) ? f : !isempty(m) ? phys_2Dvalues(basis,m,zw,reverse)  : zeros(0,0)
					   mm = !isempty(m) ? m : !isempty(f) ? modal_2Dvalues(basis,f,zw,reverse) : zeros(0,0)
					   return new(zw,reverse,deepcopy(ff),deepcopy(mm),name,basis)
	end
end



"""
Computes the L2 errors for the entire field
Parameters:
===========
h: reference size (scalar ( h1 * h2))
f1: first field
f2: second field
w1,w2: Gauss weights in each direction
"""
function get_eL2(f1::matD, f2::matD, w1::vecD, w2::vecD)
	i1,i2 = tensor_index(length(w1), length(w2))
	if size(f1,1) != length(i1)
		@info "get_eL2 !!! miss-match between data and quadrature weigths!! f1 has $(size(f1,1)) points (nodes) and we have $(length(i1)) weights"
		return 0.0
	end
	return sqrt(0.25 * sum([sum((f1[:,e] .- f2[:,e]).^2 .* w1[i1] .* w2[i2]) for e = 1:n]))
end

get_exact_2Ddata(mesh::MESH, type::String, zw::Vector{Vector{dbl}}; time = 0.0::dbl) =
	reshape(eval_2Danalytic_field(mesh_global_to_local(mesh,zw),type, time = time),length(zw[1]) * length(zw[2]) , mesh.N)

function eval_2Danalytic_field(t::Vector{Vector{dbl}},type::String;time=0.0::dbl)
	if type == "sincos"
		return sinpi.(t[1]) .* cospi.(t[2])
	elseif type == "sinex"
		return sinpi.(t[1])
	elseif type == "ctt"
		return fill(1.0, length(t[1]))
	elseif type == "poly_xy"
		return t[2]#.* t[1]
	elseif type == "poly_x+y"
		return t[2].+ t[1]
	elseif type == "poly_x2"
		return t[1].^2
	elseif type == "poly_xy2"
		return t[1] .* t[2].^2
	elseif type == "poly_x2y2"
		return t[1].^2 .* t[2].^2
	elseif type == "poly_x3y3"
		return t[1].^2 .* t[2].^3
	elseif type == "burgers_sinXpY"
		return sinpi.(t[1] .+ t[2] .- 2.0 * time)
	elseif type == "exp"
		return exp.(.-t[1].^2 .-t[2].^2)
	else
		type != "sine_xpy" && @info " Non given. I will assume analytic solution is sin(x + y)"
		return sinpi.(t[1].+ t[2])
	end
end

"""
Given a field file, reconstructs a field solution in modal form stored in a FIELD structure

Parameters:
===========
mesh file
field file

Optional:
=========
modalExp:  expansion type (default: legendre for quads and hierarchy for tris (modified legendre))
degree: polynomial degree for the expansion
reverse: by default (false), y runs first and x second. Reverse = true to switch it
Returns: a FIELD data structure
"""
function load_2Dfield(fileF::String; degree = -1, modalExp = "legendre",reverse= false)
    ## START LOADING MESH INFO
    #############################
	ext = split(fileF, ".")[end]
	if ext == "vtu" || ext == "vtk"
		return load_vtu(fileF)
	end
    data = readdlm(fileF)
	F, fInfo, fType = find_data("\$Fields", data)
	(F == nothing) && (error(" No field data !! ") && return nothing)
	n    = F[1]
	Q    = F[2]
	type = string(F[3])
	zwS1 = GZW(type, Q)
	Q    = F[4]
	type = string(F[5])
	zwS2 = GZW(type, Q)
	p    = degree == -1 ? Q - 1 : degree
	 # ASSUMINE ONLY ONE FIELD (SOLUTION) IS GIVEN
	bt   = fType[1] == 3 && modalExp == "legendre" ? "hierarchy" : modalExp
	return FIELD_2D(bt,p, fType[1] == 3, zw = [zwS1,zwS2], f = hcat(fInfo...),reverse=reverse)
end

load_modes(q::Int, modes::matD, btype::String, collapse::Bool,name=""::String) = FIELD_2D(btype,q-1, collapse,modes, name)

"""
Given a modes file, creates a FIELD data structure

Parameters:
===========
p = degree
m = total modes (x-dir + y-dir)
n = elements
btype = basis type corresponding to modes
collapse = true if triangular mesh
f = file name
Optional:
=========
name: name the field variable (pressure, u, etc)
Returns: a FIELD data structure
"""
function load_modes_file(p::Int, m::Int, n::Int, btype::String, collapse::Bool,f::String;name=""::String)
	data  = readdlm(f)
	modes = dbl.(reshape(data, m , n))
	 # ASSUMINE ONLY ONE FIELD (SOLUTION) IS GIVEN
	return FIELD_2D(btype,p,collapse,m = modes)
end

function load_vtu(fname::String)
    ext = split(fname, ".")[end]
    if ext != "vtu" && ext != "vtk"
        @error " Invalid data file  $fname needs extension .vtu or .vtk"
        return nothing
    end

	NODES    = Dict{Int, Tuple{vecD,vecI}}()
	EDGES    = Dict{vecI,vecI}()
	FACETS   = Dict{vecI,vecI}()
	ELEMENTS = Dict{Int, vecI}()
	BFACES   = Dict{vecI,Int}()
	BEDGES   = Dict{vecI,Int}()
	EMAP     = Dict{Int,Int}()

    gridreader = read_vtk(fname)
	coords = gridreader.point_coords
	conne  = gridreader.cell_connectivity
	# Elements map: adjacent nodes + coordinates
	nE      = length(conne)
	nV      = length(coords[1,:])
	nm      = [1:nV;]
	pdata   = gridreader.point_data
	fnames  = collect(keys(pdata))
	nf      = length(fnames)
	zw      = GZW("Lobatto",2)
	p       = 1
	nv      = length(conne[1])
	fields  = [FIELD_2D("hierarchy",p,nv == 3,zw=[zw,zw],n = fnames[j],f = zeros(nv,nE), m = zeros((p+1)^2,nE), reverse=true) for j = 1:nf]

	ca  = length(conne[1]) == 4 ? [1,2,4,3,1] : [1,2,3,1]
	for j = 1:nE
		ele = conne[j][ca]
		ELEMENT[j] = ele
		for f = 1 : nf
			phi = pdata[fnames[f]]
			fields[f].f[:,j] = phi[ele]  # x runs first
		end
		for v in ele
			if haskey(NODES,v)
				aux = NODES[v]
				NODES[v] = (aux[1],vcat(aux[2],tBE[j]))
			else
				NODES[v] = (xyz[:,v], [tBE[j]])
			end
		end
		for f in FMAP
			face = vcat(f,f[1])
			for e in 1:length(face)-1
				edge = sort(face[e:e+1])
				if haskey(EDGES,edge)
					push!(EDGES[edge], tBE[j])
				else
					EDGES[edge]=[tBE[j]]
				end
			end
			face = sort(ele[f])
			if haskey(FACETS,face)
				push!(FACETS[face],tBE[j])
			else
				FACETS[face] = [tBE[j]]
			end
		end
	end
	mesh = MESH(2,NODES,EDGES,FACETS,ELEMENTS,BEDGES,BFACES,FMAP,EMAP)
	return fields, mesh
end

function eval_basis(basis::EXP_2D, z1::vecD, z2::vecD)
	P     = zeros(length(z1),length(basis.m1))
	if basis.type == "pk" || basis.type == "Pk" || basis.type == "PK"
		for j = 1:length(basis.m1)
	    	c       = (0.5)^(basis.m1[j] + 1 / 2) * sqrt((2 * basis.m1[j] + 1) * (basis.m1[j] + basis.m2[j] + 1))
	        JYvals  = poly_jacobi(basis.m2[j],2 * basis.m1[j] + 1, 0.0)
			JXvals  = poly_legendre(basis.m1[j])
			P[:, j] = c .* JXvals.(z1) .* JYvals.(z2) .* ((1 .- z2).^(basis.m1[j]))
	    end
	else
		[P[:,j] = basis.J[basis.m1[j]].(z1) .* basis.J[basis.m2[j]].(z2) for j = 1:length(basis.m1)]
	end
	return P
end


function phys_at_element(data::FIELD_2D,z1::vecD,z2::vecD,ele::Int)
	phi = eval_basis(data.basis, z1, z2)
	return [sum(data.modes[:,ele] .* phi[j,:]) for j = 1 : length(z1)*length(z2)]
end

function update_element_modes!(field::FIELD_2D, zw::Vector{GZW}, e::Int)
	# set mass matrix
	f = field.f[:,e]
	Q = length(f)
	i1,i2 = tensor_index(zw[1].Q,zw[2].Q, reverse=field.reverse)
	if (Q != length(i1))
		@error " Inside modes_f :: MISMATCH BETWEEN WEIGHTS AND FIELD $Q != $(zw[1].Q * zw[2].Q)"
		return nothing
	end
	# compute (f, phi_j) for j = 1 : order
	w12  = zw[1].weights[i1] .* zw[2].weights[i2]
	phi  = eval_basis(field.basis,zw[1].nodes[i1],zw[2].nodes[i2])
	if field.basis.type == "Pk" || field.basis.type == "pk" || field.basis.type == "Pk"
		field.modes[:,e] = [sum(f .* w12 .* phi[:, j] .* 0.5 .* (1 .- zw[2].nodes[i2])) for j = 1:length(data.m1)]
	else
		field.modes[:,e] = field.basis.mass * [sum(f .* phi[:,j] .* w12) for j = 1 : length(phi[1,:])]
	end
end




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
function evaluate_2Dpoint(data::FIELD_2D,e::Int,z::vecD)
	P = eval_basis(data.basis,[z[1]],[z[2]])
	return sum(data.modes[:,e] .* P[1,:])
end


"""
Returns the evaluation points (nodes) stored in data
Parameters:
===========
a FIELD structure
"""
get_evaluation_points(data::FIELD_2D) = data.zw[1].nodes, data.zw[2].nodes


"""
Returns the evaluation weights stored in data
Parameters:
===========
a FIELD structure
"""
get_evaluation_weights(data::FIELD_2D) = data.zw[1].weights, data.zw[2].weights

"""
Modal reconstruction based on expansion (basis) and f sampled at zw12 points
Parameters:
===========
basis: a EXP_2D structure (basis, mass matrix,etc)
f: the projection function at a particular element size length(zw1) * length(zw2)
zw1: points sampled along 1st dir
zw2: points sampled along 2nd dir

Returns:
========
modes matrix with dimension [data.degree + 1, #elements]
"""
function modal_2Dvalues(basis::EXP_2D, f::matD, zw::Vector{GZW}, reverse=false::Bool)
	# set mass matrix
	Q,n = size(f)
	i1,i2 = tensor_index(zw[1].Q , zw[2].Q, reverse=reverse)
	if (Q != length(i1))
		@error " Inside modes_f :: MISMATCH BETWEEN WEIGHTS AND FIELD $Q != $(zw[1].Q * zw[2].Q)"
		return nothing
	end
	w12  = zw[1].weights[i1] .* zw[2].weights[i2]
	phi  = eval_basis(basis,zw[1].nodes[i1],zw[2].nodes[i2])
	if basis.type == "Pk" || basis.type == "pk" || basis.type == "Pk"
		return [sum(f[:, e] .* w12 .* phi[:, j] .* 0.5 .* (1 .- zw[2].nodes[i2])) for j = 1:length(basis.m1), e = 1:n]
	else
		return hcat([basis.mass * [sum(f[:,e] .* phi[:,j] .* w12) for j = 1 : length(phi[1,:])] for e = 1:n]...)
	end
end

function phys_2Dvalues(basis::EXP_2D, modes::matD, zw::Vector{GZW}, reverse=false::Bool)
	i1, i2 = tensor_index(zw[1].Q, zw[2].Q,reverse=reverse)
	P      = eval_basis(basis,zw[1].nodes[i1],zw[2].nodes[i2])
	proj   = zeros(length(i1),size(modes,2))
	@einsum proj[k,e] = modes[j,e] * P[k,j]
	return proj
end

function update_2Dphys!(data::FIELD_2D,q::Vector{String}, n::vecI)
	data.zw[1] = GZW(q[1],n[1]); data.zw[2] = GZW(q[2],n[2])
	data.f = phys_2Dvalues(data.basis, data.modes, data.zw, data.reverse)
end

function plot_2Dfield(msh::MESH,fld::FIELD_2D;f=matD(undef,0,0)::matD, showMesh=false::Bool, block=false)
	pygui(true)
	ff = isempty(f) ? fld.f : f
	fig  = plt.figure(1)
	ax   = fig.add_subplot(111,projection="3d")
	t    = mesh_global_to_local(msh,[fld.zw[1].nodes,fld.zw[2].nodes])
	surf = hcat(ff...)
	if showMesh
		for e = 1 : msh.N
			verts = element_coords(msh,e)
			plot(vcat(verts[1,:],verts[1,1]),vcat(verts[2,:],verts[2,1]), c="indigo")
		end
	end
	ax.scatter(t[1,:],t[2,:], surf, cmap="viridis", c = surf, marker = ".")
	ax.view_init(azim=0, elev=90)
	show(block=block)
	sleep(2)
	close(fig)
end
