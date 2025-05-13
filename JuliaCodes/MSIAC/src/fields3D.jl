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

struct EXP_3D
	type     :: String
	degree   :: Int
	J        :: Vector{Polynomials.Polynomial}
	mass     :: matD
	collapse :: Bool
	m1       :: vecI
	m2       :: vecI
	m3       :: vecI
	function EXP_3D(type::String,p::Int; collapse=false::Bool)
		m1,m2,m3 = tensor_index(p+1,p+1,n3=p+1)
		if type == "legendre" || type == "LEGENDRE" || type == "Legendre"
			J  = [JP.poly_legendre(i) for i = 0:p]
		else
			if p == 0
				J  = [JP.poly_legendre(0)]
			else
				PJ = [JP.poly_jacobi(i,1.0,1.0) for i = 0:p-2]
				J  = vcat(Polynomials.Polynomial([0.5,.-0.5]),
						 [Polynomials.Polynomial([0.25,0.0,.-0.25]) * PJ[i] for i = 1:p-1],
						  Polynomials.Polynomial([0.5,0.5]))
			end
		end
		# Get mass matrix
		q        = [GZW("legendre", p + 1), GZW("legendre", p + 1), GZW("legendre", p + 1)]
		t1,t2,t3 = tensor_index(p+1,p+1,n3=p+1)
		pl2      = zeros((p+1)^3,(p+1)^3)
		n1       = q[1].nodes[m1] ;	n2 = q[2].nodes[m2] ; n3 = q[3].nodes[m3]
		we       = q[1].weights[m1] .* q[2].weights[m2] .* q[3].weights[m3]
		[pl2[i,j] = sum(J[m1[i]].(n1) .* J[m2[i]].(n2) .* J[m3[i]].(n3) .* J[m1[j]].(n1) .* J[m2[j]].(n2) .* J[m3[j]].(n3) .* we) for i = 1 : length(m1), j = 1:length(m1) ]
		mij  = inv(pl2)
		return new(type,p,J,mij,collapse,m1,m2,m3)
	end
end

mutable struct FIELD_3D
	zw      :: Vector{GZW}
	reverse :: Bool
	f       :: matD
	name    :: String
	modes   :: matD
	basis   :: EXP_3D
	function FIELD_3D(type::String,p::Int ; reverse=false::Bool,zw = GZW[]::Vector{GZW},
		              f = zeros(0,0)::matD, m = zeros(0,0)::matD,n = ""::String)
				      	basis = EXP_3D(type, p)
				      	qt = ["legendre","legendre", "legendre"] ; Q = p + 2
		   		      	isempty(zw) && (zw = [GZW(qt[1],Q) GZW(qt[2],p) GZW(qt[3],Q)])
				   		ff = !isempty(f) ? f : !isempty(m) ? phys_3Dvalues(basis,m,zw,reverse)  : zeros(0,0)
				   		mm = !isempty(m) ? m : !isempty(f) ? modal_3Dvalues(basis,f,zw,reverse) : zeros(0,0)
				   		return new(zw,reverse,ff,"u",mm,basis)
	end
end


get_exact_3Ddata(mesh::MESH, type::String, zw::Vector{Vector{dbl}}; time = 0.0::dbl) =
	reshape(eval_3Danalytic_field(mesh_global_to_local(mesh,zw),type, time = time), length(zw[1]) * length(zw[2]) * length(zw[3]),mesh.N)

function eval_3Danalytic_field(t::Vector{Vector{dbl}},type::String;time=0.0::dbl)
	if type == "sine_x+y+z"
		return sinpi.(t[1] + t[2] + t[3])
	elseif type == "ctt"
		return fill(1.0, length(t[1]))
	elseif type == "poly_xyz"
		return t[3].*t[2].* t[1]
	elseif type == "poly_xyz2"
		return (t[3].*t[2].* t[1]).^2
	elseif type == "poly_x+y+z"
		return t[1].+t[2].+ t[3]
	else
		@info " Non given. I will return zeros"
		return zeros(length(t[1]))
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
function load_3Dfield(fileF::String; degree = -1, modalExp = "legendre", reverse= false)

	## START LOADING MESH INFO
    data = readdlm(fileF)
	F, fInfo, fType = find_data("\$Fields", data)
	(F == nothing) && (error(" No field data !! ") && return nothing)
	n    = F[1]
	zwS  = [GZW(string(F[3]), F[2]),GZW(string(F[5]), F[4]),GZW(string(F[7]), F[6]) ]
	p    = degree == -1 ? zwS[1].Q - 1 : degree
	 # ASSUMINE ONLY ONE FIELD (SOLUTION) IS GIVEN
	return FIELD_3D(modalExp,p,zw=zwS, f = hcat(fInfo...),reverse=reverse)
end

load_modes(q::Int, modes::matD, type::Int,name=""::String) = FIELD_3D(type,q-1,modes,name)

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
function load_modes_file(p::Int, m::Int, n::Int, type::Int,f::String;name=""::String)
	data  = readdlm(f)
	modes = dbl.(reshape(data, m , n))
	 # ASSUMINE ONLY ONE FIELD (SOLUTION) IS GIVEN
	return FIELD_3D(type,p,m = modes)
end

function eval_basis(basis::EXP_3D, z1::vecD, z2::vecD, z3::vecD)
	P       = zeros(length(z1),length(basis.m1))
	[P[:,j] = basis.J[basis.m1[j]].(z1) .* basis.J[basis.m2[j]].(z2) .* basis.J[basis.m3[j]].(z3) for j = 1:length(basis.m1)]
	return P
end


function phys_values(basis::EXP_3D, modes::matD, zw::Vector{GZW})
	i1, i2, i3 = tensor_index(zw[1].Q, zw[2].Q,n3=zw[3].Q)
	P      = eval_basis(basis,zw[1].nodes[i1],zw[2].nodes[i2],zw[3].nodes[i3])
	proj   = zeros(length(i1),size(modes,2))
	@einsum proj[k,e] = modes[j,e] * P[k,j]
	return proj
end

function update_3Dphys!(data::FIELD_3D,q::Vector{String}, n::vecI)
	data.zw[1] = GZW(q[1],n[1]); data.zw[2] = GZW(q[2],n[2]) ; data.zw[3] = GZW(q[3],n[3])
	data.f = phys_values(data.basis, data.modes, data.zw)
end

function phys_at_element(data::FIELD_3D,z1::vecD,z2::vecD,z3::vecD, ele::Int)
	phi = eval_basis(data.basis, z1, z2, z3)
	return [sum(data.modes[:,ele] .* phi[j,:]) for j = 1 : length(z1)*length(z2)*length(z3)]
end

function update_element_modes!(field::FIELD_3D, zw::Vector{GZW}, e::Int)
	# set mass matrix
	f = field.f[:,e]
	Q = length(f)
	i1,i2,i3 = tensor_index(zw[1].Q,zw[2].Q,n3 = zw[3].Q, field.reverse)
	if (Q != length(i1))
		@error " Inside modes_f :: MISMATCH BETWEEN WEIGHTS AND FIELD "
		return nothing
	end
	# compute (f, phi_j) for j = 1 : order
	w123 = zw[1].weights[i1] .* zw[2].weights[i2].* zw[3].weights[i3]
	phi  = eval_basis(field.basis,zw[1].nodes[i1],zw[2].nodes[i2], zw[3].nodes[i3])
	field.modes[:,e] = field.basis.mass * [sum(f .* phi[:,j] .* w123) for j = 1 : length(phi[1,:])]
end

"""
Updates the mode in data from the field values (also in data)
"""
update_modes!(data::FIELD_3D) = (data.modes = modal_values(data.basis, data.f, data.zw, data.reverse))


"""
Physical (field) values given an element an its local coordinates (in reference system)
Parameters:
===========
FIELD_3D structure
e: the element ID
z: points coordinates in reference system (chi: [-1,1]^2)

Returns:
========
Field values at point
"""
function evaluate_3Dpoint(data::FIELD_3D,e::Int,z::vecD)
	modes  = data.modes[:,e]
	P      = eval_basis(data.basis,[z[1]],[z[2]],[z[3]])
	return sum(modes .* P[1,:])
end


"""
Returns the evaluation points (nodes) stored in data
Parameters:
===========
a FIELD structure
"""
get_evaluation_points(data::FIELD_3D) = data.zw[1].nodes, data.zw[2].nodes, data.zw[3].nodes

"""
Returns the evaluation weights stored in data
Parameters:
===========
a FIELD structure
"""
get_evaluation_weights(data::FIELD_3D) = data.zw[1].weights, data.zw[2].weights, data.zw[3].weights



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
function modal_3Dvalues(basis::EXP_3D, f::matD, zw::Vector{GZW}, reverse=false::Bool)
	# set mass matrix

	Q,n = size(f)
	i1,i2,i3 = tensor_index(zw[1].Q , zw[2].Q, n3=zw[3].Q, reverse=reverse)
	if (Q != length(i1))
		@error " Inside modes_f :: MISMATCH BETWEEN WEIGHTS AND FIELD $Q != $(zw[1].Q * zw[2].Q* zw[3].Q) )"
		return nothing
	end
	w123 = zw[1].weights[i1] .* zw[2].weights[i2] .* zw[3].weights[i3]
	phi  = eval_basis(basis,zw[1].nodes[i1],zw[2].nodes[i2],zw[3].nodes[i3])
	return hcat([basis.mass * [sum(f[:,e] .* phi[:,j] .* w123) for j = 1 : length(phi[1,:])] for e = 1:n]...)
end

function plot_3Dfield(msh::MESH,fld::FIELD_3D;f=matD(undef,0,0)::matD, showMesh=false::Bool, block=false)
	@info "PLOT 3D "
	pygui(true)
	ff   = isempty(f) ? fld.f : f
	fig  = plt.figure(figsize=(10, 10))
	ax   = fig.add_subplot(111,projection="3d")

	t    = mesh_global_to_local(msh,[fld.zw[1].nodes, fld.zw[2].nodes, fld.zw[3].nodes] )
	surf = hcat(ff...)
	if showMesh
		for e = 1 : msh.N
			verts = element_coords(msh,e)
			vf = hcat(verts...)
			for f in msh.FMAP
				fl = vcat(f,f[1])
				plot(vf[1,fl], vf[2,fl], vf[3,fl], c="indigo")
			end
		end
	end
	#ax.scatter(t[1,:],t[2,:], t[3,:],  cmap="viridis", c = surf, marker = ".", s = 10)
	show(block=block)

end
