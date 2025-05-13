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

"""
Creates a 2D mesh file
Parameters:
===========
type: 3 = triangles, 4 = quads
nX: number of elements in the x direction
nY: number of elements in the y direction
a: domain boundary for x
b: domain boundary for y
periodic: if true, tags periodic vertices and assigns the map to its period
name: name of the mesh file

Optional:
=========
yFirst: indeces go vertically (default, we run x first)
pert: perturbate the uniform grid
structured: true or false
"""
function create_2D_mesh(type::Int, nX::Int, nY::Int, a::vecD, b::vecD, periodic::Bool, name::String ;
	                    yFirst = false , pert = [0.0,0.0], structured = true)

						if type==3 && structured == true && iszero(pert)
							return create_2D_tri_mesh(nX, nY, yFirst, a, b, periodic, name)
						end

	points = zeros((nX+1)*(nY+1),2)
	dx     = (a[2]-a[1])/nX
	dy     = (b[2]-b[1])/nY


	t      = yFirst ? [b[1] .+ dy*(j-1) for j = 1:nY+1] :
	                  [a[1] .+ dx*(j-1) for j = 1:nX+1]

	n      = yFirst ? [nX+1,nY+1] : [nY+1,nX+1]

	epsX = pert[1] == 0 ? zeros(n[1]) : rand(Uniform(-dx * pert[1],dx * pert[1]),n[1])
	epsY = pert[2] == 0 ? zeros(n[1]) : rand(Uniform(-dy * pert[2],dy * pert[2]),n[1])

	id  = 1
	row = 0
	for j = 1:n[1]
		bp = [1:n[2];]
		ip = bp[2:end-1]
		if !structured
			(pert[1] > 0.0 ) && (epsX = rand(Uniform(-dx * pert[1],dx * pert[1]),n[2]))
			(pert[2] > 0.0 ) && (epsY = rand(Uniform(-dy * pert[2],dy * pert[2]),n[2]))
		end
		if yFirst
			points[row.+bp,1]  .= a[1] .+ dx*(j-1)
			points[row.+bp,2]  .= t
			points[row.+ip,2] .+= epsY[ip]

			if j != 1 && j != n[1]
				if structured
					points[row.+bp,1] .+= epsX[id]
					id += 1
				else
					points[row.+bp,1] += epsX
				end
			end
		else
			points[row.+bp,1] .= t
			points[row.+bp,2] .= b[1] .+ dy*(j-1)
			points[row.+ip,1] .+= epsX[ip]

			if j != 1 && j != n[1]
				if structured
					points[row.+bp,2] .+= epsY[id]
					id += 1
				else
					points[row.+bp,2] += epsY
				end
			end
		end
		row += n[2]
	end
	open(name, "w") do io
		println(io,"\$Nodes")
	    println(io,(nX+1) * (nY+1))
		ji = 0
		for i = 1:nY + 1
			for j = 1:nX+1
				ti = 1
				ji += 1
				if i == 1 || i == nY + 1
					ti = (j == 1 || j == nX + 1) ? -1 : 0
				elseif j == 1 || j == nX + 1
					ti = 0
				end
				println(io,ji," ", ti," ",points[ji,1]," ",points[ji,2]," ",0.0)
			end
		end
		println(io,"\$EndNodes")
		## ELEMENT MAP
		println(io,"\$Elements")
		if type == 3
			tess = basic_triangulation(points, [1:length(points[:,1]);])
			sm   = size(tess,1)
			println(io,sm)
			[ println(io,i," ", 3," ", tess[i][1]," ", tess[i][2]," ", tess[i][3]) for i = 1:sm]

			#tess = delaunay(points)
			#sm = size(tess.simplices)
			#println(io,sm[1])
			#for i = 1:sm[1]
			#	tri = tess.simplices[i,:]
			#	println(io,i," ", 3," ",tri[1]," ",tri[2]," ",tri[3])
			#end
		else
			if yFirst
		        println(io,nX * nY)
				q  = [0 nY+1 nY+2 1]
				for i = 1:nX
					[println(io,(i-1) * nY +j," ", 4," ",j+q[1]," ",j+q[2]," ",j+q[3]," ",j+q[4]," ") for j = 1:nY]
					q .+= nY+1
		        end
			else
				println(io,nX * nY)
				q  = [0 1 nX+2 nX+1]
				for i = 1:nY
					[println(io,(i-1) * nX +j," ", 4," ",j+q[1]," ",j+q[2]," ",j+q[3]," ",j+q[4]," ") for j = 1:nX]
					q .+= nX+1
		        end
			end
		end
		println(io,"\$EndElements")
		println(io,"\$Boundaries")
		println(io,4)
		print(io,"1 1 ")
		n = yFirst ? nY : nX
        [print(io, i," ") for i =1:n+1]
		print(io,"\n")

		print(io,"2 1 ")
        [print(io,i," ") for i =(nX + 1) * (nY + 1) - n:(nX + 1) * (nY + 1)]
		print(io,"\n")

		print(io,"3 1 ")
        [print(io, i," ") for i =1:n+1:(nX +1) * (nY+1)]
		print(io,"\n")

		print(io,"4 1 ")
        [print(io, i," ") for i =n+1:n+1:(nX +1) * (nY+1)]
		print(io,"\n")
		println(io,"\$EndBoundaries")
		println(io,"\$BCs")
		if periodic
			println(io,"2")
			println(io,"1 -1 1 2")
			println(io,"2 -1 3 4")
		else
			println(io,"4")
			println(io,"1 1 ")
			println(io,"2 1 ")
			println(io,"3 1 ")
			println(io,"4 1 ")
		end
		println(io,"\$EndBCs")
    end
	@info " Written in $name"
    return name
end

function create_2D_tri_mesh(nX::Int, nY::Int, orient::Bool, a::vecD, b::vecD, periodic::Bool, name::String)
    open(name, "w") do io
        println(io, "\$Nodes")
        println(io, (nX + 1) * (nY + 1))
        hx = [a[1] + j * (a[2] - a[1]) / nX for j = 0:nX]
        hy = [b[1] + j * (b[2] - b[1]) / nY for j = 0:nY]
        for i = 1:nY+1
            for j = 1:nX+1
                ti = 1
                if i == 1 || i == nY + 1
                    ti = (j == 1 || j == nX + 1) ? -1 : 0
                elseif j == 1 || j == nX + 1
                    ti = 0
                end
                println(io, (i - 1) * (nX + 1) + j, " ", ti, " ", hx[j], " ", hy[i], " ", 0.0)
            end
        end
        println(io, "\$EndNodes")
        println(io, "\$Elements")
        println(io, 2 * nX * nY)
        #Decide if diagonals have angle pi/4 or 3pi/4
        if orient
            q1 = [0 1 nX + 2]
            q2 = [0 nX + 2 nX + 1]
        else
            q2 = [1 nX + 1 0]
            q1 = [1 nX + 2 nX + 1]
        end
        for i = 1:nY
            for j = 1:nX
                println(io, 2 * (i - 1) * nX + 2 * j - 1, " ", 3, " ", j + q2[3], " ", j + q2[1], " ", j + q2[2], " ")
                println(io, 2 * (i - 1) * nX + 2 * j, " ", 3, " ", j + q1[1], " ", j + q1[2], " ", j + q1[3], " ")
            end
            q1 .+= nX + 1
            q2 .+= nX + 1
        end
        println(io, "\$EndElements")
        println(io, "\$Boundaries")

        #South
        println(io, 4)
        print(io, "1 1 ")
        [print(io, i, " ") for i = 1:nX+1]
        print(io, "\n")
        #North
        print(io, "2 1 ")
        [print(io, i, " ") for i = (nX+1)*(nY+1)-nX:(nX+1)*(nY+1)]
        print(io, "\n")
        #West
        print(io, "3 1 ")
        [print(io, i, " ") for i = 1:nX+1:(nX+1)*(nY+1)]
        print(io, "\n")
        #East
        print(io, "4 1 ")
        [print(io, i, " ") for i = nX+1:nX+1:(nX+1)*(nY+1)]
        print(io, "\n")
        println(io, "\$EndBoundaries")
        println(io, "\$BCs")
        if periodic
            println(io, "2")
            println(io, "1 -1 1 2")
            println(io, "2 -1 3 4")
        else
            println(io, "4")
            println(io, "1 1 ")
            println(io, "2 1 ")
            println(io, "3 1 ")
            println(io, "4 1 ")
        end
        println(io, "\$EndBCs")
    end
    @info "Written in $name"
    return name
end


"""
Creates a structured 3D mesh file using GMSH
Parameters:
===========
mt: type hex or tet
nX: number of elements in the x direction
nY: number of elements in the y direction
nZ: number of elements in the z direction
dX: domain boundary for x
dY: domain boundary for y
dZ: domain boundary for y

Return:
=======
mesh structure
"""
function gmsh_3Dmesh(mt::String, nX::Int, nY::Int, nZ::Int, dX::vecD, dY::vecD, dZ::vecD, fname::String,showmesh=false)

	gmsh.initialize()
	factory = gmsh.model.geo
	gmsh.model.add(mt)

	lc = 1
	# points
	factory.addPoint(dX[1], dY[1], dZ[1],lc, 1)
	factory.addPoint(dX[2], dY[1], dZ[1],lc, 2)
	factory.addPoint(dX[2], dY[2], dZ[1],lc, 3)
	factory.addPoint(dX[1], dY[2], dZ[1],lc, 4)

	#lines
	l1 = factory.addLine(1, 2, 1)
	l2 = factory.addLine(2, 3, 2)
	l4 = factory.addLine(4, 1, 3)
	l3 = factory.addLine(3, 4, 4)

	#loop
	factory.addCurveLoop([l1,l2,l3,l4], 1)

	#surface
	plane = factory.addPlaneSurface([1], 1)

	factory.extrude([(2, plane)], 0., 0., dZ[2],[nZ], Cdouble[], mt == "hex")
	factory.synchronize()
	meshFact = gmsh.model.mesh

	if mt == "hex"
		meshFact.setTransfiniteCurve(l1, nX+1)
		meshFact.setTransfiniteCurve(l2, nY+1)
		meshFact.setTransfiniteCurve(l3, nX+1)
		meshFact.setTransfiniteCurve(l4, nY+1)

		meshFact.setTransfiniteSurface(plane)
		meshFact.generate(2)
		#gmsh.model.geo.mesh.setRecombine(2, 1)
		mt == "hex" && meshFact.recombine()
	end
	meshFact.generate(3)
	showmesh && gmsh.fltk.run()
	file = fname*".msh"
	gmsh.write(file)
	gmsh.finalize()
	return file
end

"""
Creates a structured 2D mesh file using GMSH
Parameters:
===========
mt: type hex or tet
nX: number of elements in the x direction
nY: number of elements in the y direction
dX: domain boundary for x
dY: domain boundary for y

Return:
=======
mesh structure
"""
function gmsh_2Dmesh(mt::String, nX::Int, nY::Int, dX::vecD, dY::vecD,fname::String, showmesh=false)

	gmsh.initialize()
	factory = gmsh.model.geo
	gmsh.model.add(mt)

	lc = 1e-1
	# points
	factory.addPoint(dX[1], dY[1], 0.0,lc, 1)
	factory.addPoint(dX[2], dY[1], 0.0,lc, 2)
	factory.addPoint(dX[2], dY[2], 0.0,lc, 3)
	factory.addPoint(dX[1], dY[2], 0.0,lc, 4)

	#lines
	l1 = factory.addLine(1, 2, 1)
	l2 = factory.addLine(2, 3, 2)
	l4 = factory.addLine(4, 1, 3)
	l3 = factory.addLine(3, 4, 4)

	#loop
	factory.addCurveLoop([l1,l2,l3,l4], 1)

	#surface
	plane = factory.addPlaneSurface([1], 1)
	factory.synchronize()


	meshFact = gmsh.model.mesh
	meshFact.setTransfiniteCurve(l1, nX+1)
	meshFact.setTransfiniteCurve(l2, nY+1)
	meshFact.setTransfiniteCurve(l3, nX+1)
	meshFact.setTransfiniteCurve(l4, nY+1)

	meshFact.setTransfiniteSurface(plane)
	meshFact.generate(2)
	mt == "quad" && meshFact.recombine()
	showmesh && gmsh.fltk.run()
	file = fname*".msh"
	gmsh.write(file)
	gmsh.finalize()
	return file
end

num_vertex(m::MESH, j::Int)  = length(m.ELEMENTS[j])
eFace(m::MESH,e::Int,j::Int) = m.ELEMENTS[e][m.FMAP[j]]
coords(m::MESH, j::Int)      = m.NODES[j][1]

function adjacent(m::MESH, e::Int, f::Union{vecI,Int})
	adj = Int[]
	fg  = typeof(f) == Int ? sort(m.ELEMENTS[e][m.FMAP[f]]) : sort(m.ELEMENTS[e][f])
	if length(fg) == 2
		edge = haskey(m.BEDGES, fg) ? m.BEDGES[fg] : fg
		adj  = unique!(filter(x-> x != e, m.EDGES[edge]))
	else
		face = haskey(m.BFACES, fg) ? m.BFACES[fg] : fg
		adj  = filter(x-> x != e, m.FACETS[face])
	end
	isempty(adj) ? (return -1) : length(adj) == 1 ? (return adj[1]) : (return adj)
end

periodic_edge(m::MESH, e::Int,j::Int)  = haskey(m.BEDGES,sort( m.ELEMENTS[e][[j,mod1(j+1,length(m.ELEMENTS[e]))]] ))
periodic_edge(m::MESH, e::Int,j::vecI) = haskey(m.BEDGES,sort(j))



"""
Returns the coordinates of the element vertices (ordered)
"""
element_coords(m::MESH,e::Int) = [m.NODES[k][1] for k in m.ELEMENTS[e]]

"""
Returns the physical coordinates of the element average centre
Parameters:
===========
data: mesh & field data structure (twoD_mesh_field)
eID: element id
"""
element_centre(mesh::MESH, eID::Int) = (sum(element_coords(mesh,eID)) / length(mesh.ELEMENTS[eID]))
element_centre(v::Vector{Vector{dbl}}) = sum(v) / length(v)

function max_length(mesh::MESH, eID::Union{uint, Int})
	xyz = hcat(element_coords(mesh,eID)...)
	return norm([maximum(xyz[j,:])-minimum(xyz[j,:]) for j = 1 : 3])
end

function element_lengths(mesh::MESH, eID::Union{uint, Int})
	xyz   = element_coords(mesh,eID) ; n = mesh.type
	sides = n <= 4 ? [1:n;1] : vcat([1:4;1],[5:8;5],[1,5,8,4,1],[2,6,7,8,2])
	return [norm(xyz[sides[j+1]] - xyz[sides[j]]) for j = 1:length(sides)-1]
end

get_frame(mesh::MESH) = mesh.dim == 2 ? [[1.0,0.0,0.0] [0.0,1.0,0.0]] : [[1.0,0.0,0.0] [0.0,1.0,0.0] [0.0,0.0,1.0]] #WARNING being ignored. returning cartesian axis

function load_mesh(fileM::String; ordered = true::Bool, structured=false::Bool)
	ext = split(fileM, ".")[end]
	ext == "msh" && return load_gmsh(fileM)
	data = readdlm(fileM)
	idN  = 1
	N,nInfo,nType,aux = find_data("\$Nodes", data)
	ordered && (idN += aux)
	(N == nothing) && (error("Nodes not found in mesh file\n") && return nothing)
	E,eInfo, eType,aux = find_data("\$Elements", data[idN:end,:])
	ordered && (idN += aux)
	(E == nothing) && (error("Elements not found in mesh file\n") && return nothing)
	L,lInfo, lType,aux = find_data("\$Boundaries",data[idN:end,:])
	ordered && (idN += aux)
	(L == nothing) && (error("Boundaries not found in mesh file\n") && return nothing)
	B,bInfo, bType,aux = find_data("\$BCs",data[idN:end,:])
	ordered && (idN += aux)
	(B == nothing) && (error("BCs not found in mesh file\n") && return nothing)

	NODES    = Dict{Int, Tuple{vecD,vecI}}()
	EDGES    = Dict{vecI,vecI}()
	FACETS   = Dict{vecI,vecI}()
	ELEMENTS = Dict{Int,vecI}()
	BNODES   = Dict{Int,vecI}()
	BEDGES   = Dict{vecI,vecI}()
	BFACES   = Dict{vecI,vecI}()

	ca       = vcat([1:eType[1];],1)
	FMAP     = [ca[j:j+1] for j = 1 : length(ca)-1]

	# create a boundary map -> periodic vertex index
	periodic = false
	for b = 1:B[1]
		bType[b] != -1 && continue
		periodic = true
		b1 = lInfo[bInfo[b][1]] ; b2 = lInfo[bInfo[b][2]]
		for j = 1:length(b1)
			if haskey(BNODES,b1[j])
				push!(BNODES[b1[j]], b2[j])
			else
				BNODES[b1[j]] = [b2[j]]
			end
			if haskey(BNODES,b2[j])
				push!(BNODES[b2[j]], b1[j])
			else
				BNODES[b2[j]] = [b1[j]]
			end
			j == length(b1) && continue
			BEDGES[sort(b1[j:j+1])] = sort(b2[j:j+1])
			BEDGES[sort(b2[j:j+1])] = sort(b1[j:j+1])
		end
	end
	for (b,k) in BNODES, kj in k
		haskey(BNODES,kj) && ( BNODES[b] = unique!(vcat(BNODES[b], filter(x -> x != b, BNODES[kj]))) )
	end
	for j = 1:E[1]
		face        = eInfo[j]
		ELEMENTS[j] = face
		for f in FMAP
			vp = face[f]
			if haskey(NODES,vp[1])
				aux = NODES[vp[1]]
				NODES[vp[1]] = (aux[1],vcat(aux[2],j))
			else
				NODES[vp[1]] = (nInfo[vp[1]],[j])
			end
			edge = sort(vp)
			if haskey(FACETS,edge)
				push!(FACETS[edge], j)
			else
				FACETS[edge]=[j]
			end
		end
	end
	return MESH(2,length(ELEMENTS),FMAP,ELEMENTS, NODES,FACETS,FACETS,BNODES=BNODES, BEDGES=BEDGES, structured = structured, type = eType[1], toler=1.e-13 )
end

function load_gmsh(fname::String, periodic=false::Bool)

	ext = split(fname, ".")[end]
    if ext != "msh"
        @error " Invalid data file $fname needs extension .msh for GMSH data"
        return nothing
    end

	NODES    = Dict{Int, Tuple{vecD,vecI}}()
	EDGES    = Dict{vecI,vecI}()
	FACETS   = Dict{vecI,vecI}()
	ELEMENTS = Dict{Int,vecI}()
	BFACES   = Dict{vecI,vecI}()

	BEDGES   = Dict{vecI,vecI}()

	gmsh.initialize()
	gmsh.open(fname)
	gmsh.model.setFileName(fname)
	gm = gmsh.model.mesh

	# the nodes
	aux = gm.getNodes()
	nT  = Int.(aux[1])
	nV  = length(nT)
	xyz = reshape(aux[2],3,nV)

	# mesh 1D boundaries
	aux = gm.getElements(1)
	nBE = eTagDict[aux[1][1]] # number of vertices per line
	tBE = Int.(aux[2][1])   # tags
	lBE = reshape(Int.(aux[3][1]), nBE, length(tBE)) # vertex list per tBE
	[BEDGES[sort(lBE[:,j])] = sort(lBE[:,j]) for j = 1: length(tBE)]

	# mesh 2D boundaries
	dim = 3

	if !isempty(gm.getElements(3)[1])
		aux = gm.getElements(2)
		nBF = eTagDict[aux[1][1]]# number of vertices per facet
		tBF = Int.(aux[2][1])   # tags
		lBF = reshape(Int.(aux[3][1]), nBF, length(tBF)) # vertex list per tBF

		[BFACES[sort(lBF[:,j])] = sort(lBF[:,j]) for j = 1: length(tBF)]
		aux = gm.getElements(3)
	else
		dim = 2
		aux = gm.getElements(2)
	end

	eShape = aux[1][1]
	# mesh elements
	nBE = eTagDict[eShape] # vertices per element
	bias = Int(aux[2][1][1]) -1
	tBE  = Int.(aux[2][1]) .- bias   # tags (element IDs)
	nE   = length(tBE)
	lBE  = reshape(Int.(aux[3][1]), nBE, nE) # vertex list per tBE
	FMAP = Int[]
	if eShape <= 3
		ca   = vcat([1:nBE;],1)
		FMAP = [ca[j:j+1] for j = 1 : length(ca)-1]
	elseif eShape == 5
		FMAP = [[1,2,3,4],[5,6,7,8],[1,5,8,4],[2,6,7,3],[1,5,6,2],[8,7,3,4]]
	else
		FMAP = [[1,2,3],[2,3,4],[3,1,4],[1,4,2]]
	end

	for j = 1:nE
		ele = lBE[:,j]
		ELEMENTS[tBE[j]] = ele
		for v in ele
			if haskey(NODES,v)
				aux = NODES[v]
				NODES[v] = (aux[1],vcat(aux[2],tBE[j]))
			else
				NODES[v] = (xyz[:,v], [tBE[j]])
			end
		end
		for f in FMAP
			face = ELEMENTS[tBE[j]][vcat(f,f[1])]
			for e in 1:length(face)-1
				edge = sort(face[e:e+1])
				if eShape <= 3
					if haskey(FACETS,edge)
						push!(FACETS[edge], tBE[j])
					else
						FACETS[edge]=[tBE[j]]
					end
				else
					if haskey(EDGES,edge)
						tBE[j] ∉ EDGES[edge] && push!(EDGES[edge], tBE[j])
					else
						EDGES[edge]=[tBE[j]]
					end
				end
			end
			eShape <= 3 && continue
			sort!(unique!(face))
			if haskey(FACETS,face)
				push!(FACETS[face],tBE[j])
			else
				FACETS[face] = [tBE[j]]
			end
		end
	end
	gmsh.finalize()
	return MESH(dim,length(ELEMENTS),FMAP,ELEMENTS,NODES,EDGES,FACETS,BFACES=BFACES,BEDGES=BEDGES,structured=true, type=nBE, toler = 1.e-11)
end

function load_trixi_mesh(msh::Any)

	NODES    = Dict{Int, Tuple{vecD,vecI}}()
	EDGES    = Dict{vecI,vecI}()
	FACETS   = Dict{vecI,vecI}()
	ELEMENTS = Dict{Int,vecI}()
	BFACES   = Dict{vecI,vecI}()
	BEDGES   = Dict{vecI,vecI}()

	# the nodes
	nV  = msh.n_corners
	xyz = vcat(msh.corners, zeros(nV)')

	# dim for meshes means 2 = 2D elements with points in 3D (z = 0.0) and 3 = 3D elements
	dim = length(msh.corners[:,1])
	vE  = length(msh.element_node_ids[:,1]) #WARNING assuming all elements same type

	FMAP   = Int[] #ccwise edges

	eShape = dim == 2 ? vE - 1 : vE  # GMSH: triangles = 2, quads = 3, tets = 4, hexes = 5

	if eShape <= 3
		ca   = vcat([1:vE;],1)
		FMAP = [ca[j:j+1] for j = 1 : length(ca)-1]
	else
		throw(@error " 3D not implemented yet....")
	end

	nE = msh.n_elements

	for j = 1:nE
		ele         = msh.element_node_ids[:,j]
		ELEMENTS[j] = ele

		for v in ele
			if haskey(NODES,v)
				aux = NODES[v]
				NODES[v] = (aux[1],vcat(aux[2],j))
			else
				NODES[v] = (xyz[:,v], [j])
			end
		end

		for f in FMAP
			face = ELEMENTS[j][vcat(f,f[1])]
			for e in 1:length(face)-1
				edge = sort(face[e:e+1])
				if eShape <= 3
					if haskey(FACETS,edge)
						push!(FACETS[edge], j)
					else
						FACETS[edge]=[j]
					end
				else
					throw(@error " 3D not implemented yet ")
					#=if haskey(EDGES,edge)
						j ∉ EDGES[edge] && push!(EDGES[edge],j)
					else
						EDGES[edge]=[j]
					end=#
				end
			end
			eShape <= 3 && continue
			sort!(unique!(face))
			if haskey(FACETS,face)
				push!(FACETS[face],j)
			else
				FACETS[face] = [j]
			end
		end
	end

	#BOUNDARY STUFF !!!
	return MESH(dim,length(ELEMENTS),FMAP,ELEMENTS,NODES,FACETS,FACETS,BFACES=BFACES,BEDGES=BEDGES,structured=false, type=vE, toler = 1.e-13)
end

function plot_mesh(m::MESH ; labels=false::Bool, block=false::Bool)
	pygui(true)
	fig = plt.figure()
	ax  = m.dim == 3 ? fig.add_subplot(projection="3d") : fig.gca()
	ax.axis("off") ; ax.grid(b=nothing)
	for (e,v) in m.ELEMENTS
		if labels
			cent = element_centre(m,e)
			m.dim == 3 ? ax.text(cent[1],cent[2],cent[3],string(Int(e)),fontsize=16) :
				         ax.text(cent[1],cent[2],string(Int(e)))
		end
		for f in m.FMAP
			poly = m.ELEMENTS[e][vcat(f,f[1])]
			vs   = hcat([coords(m,k) for k in poly]...)
			m.dim == 3 ? ax.plot3D(vs[1,:],vs[2,:],vs[3,:]) :
				         ax.plot(vs[1,:],vs[2,:])
		end
		!labels && continue
		for j in v
			xyz = coords(m,j)
			m.dim == 3 ? ax.text(xyz[1],xyz[2],xyz[3],string(j), c = "r",fontsize=16) :
					     ax.text(xyz[1],xyz[2],string(j), c = "r")
		end
	end
	show(block=block)

	sleep(2)
	close(fig)
end

function plot_boundary(m::MESH ; labels=false::Bool, block=false::Bool)
	pygui(true)
	fig = plt.figure()
	ax  = m.dim == 3 ? fig.add_subplot(projection="3d") : fig.gca()
	ax.axis("off") ; ax.grid(b=nothing)
	for (e,v) in m.BEDGES
		vs   = hcat([coords(m,k) for k in e]...)
		m.dim == 3 ? ax.plot3D(vs[1,:],vs[2,:],vs[3,:]) :  ax.plot(vs[1,:],vs[2,:])
		!labels && continue
		for j = 1:length(e)
			m.dim == 2 ? ax.text(vs[1,j],vs[2,j],        string(e[j]), c = "r") :
			  	     	 ax.text(vs[1,j],vs[2,j],vs[3,j],string(e[j]), c = "r")
		end
	end
	if m.dim == 3
		for (e,v) in m.BFACES
			ele = m.ELEMENTS[m.FACETS[v][1]]
			aux = findfirst(x->issubset(x,v), [ele[j] for j in m.FMAP])
			face = ele[m.FMAP[aux]]
			vs  = hcat([coords(m,j) for j in face]...)
			if labels
				for j = 1:length(face)
					ax.text(vs[1,j],vs[2,j],vs[3,j],string(face[j]), c = "r")
				end
			end
			id = [1:length(face);]
			tri = [[0,1,2],[0,2,3]]
			ax.plot_trisurf(vs[1,id],vs[2,id],vs[3,id], triangles=tri, color="red",edgecolor="none", antialiased=false, linewidth=0,alpha=0.15, shade=false, zorder = 1)
			id = vcat([1:length(face);],1)
			ax.plot3D(vs[1,id],vs[2,id],vs[3,id],color="indigo",lw = 1,zorder = 6)
			!labels && continue
			ce = [sum(vs[j,:]) / length(face) for j = 1 : 3]
			ax.text(ce[1],ce[2],ce[3],string(m.FACETS[e][1]))

		end
	end
	show(block=block)
	sleep(2)
	close(fig)
end

function tet_to_hex(eta)
	x = zeros(3)
	ftol = 1.e-12
	if abs(z[3]-1.0) < ftol # Very top point of the tetrahedron
	  	x[1] = -1.0; x[2] = -1.0; x[3] = z[3];
	else
		if abs(z[2]-1.0) < ftol
			x[1] = -1.0;
		elseif abs(z[2] + z[3]) < ftol
			x[1] = -1.0;
		else
			x[1] = 2.0 * (1.0 + z[1]) / (-z[2]-z[3]) - 1.0;
		end
  		x[2] = 2.0 * (1.0 + z[2]) / (1.0-z[3]) - 1.0;
  		x[3] = z[3];
	end
	return x
end


function hex_to_tet(eta)
	te1 = (1.0 .+ eta[1]) .* (1.0 .- eta[3]) .* 0.5 .- 1.0;
	xi1 = (1.0 .+ te1)    .* (1.0 .- eta[2]) .* 0.5 .- 1.0;
	xi2 = (1.0 .+ eta[2]) .* (1.0 .- eta[3]) .* 0.5 .- 1.0;
	return [xi1,xi2,eta[3]]
end


"""
Given a set of quadrature points and a quadrilateral, computes the
2D quadrature region
Params:
-------
z: quadrature points
v: element vertices
Returns:
--------
Physical coordinates and Jacobian
"""
function local_region(z::Vector{Vector{dbl}}, v::Vector{Vector{dbl}})
	d = length(v[1])
	n = length(v)
	if size(z,1) == 2
		i1,i2 = tensor_index(length(z[1]),length(z[2]))
		mz1   = 0.5 .*(1.0 .- z[1][i1]) ; pz1 = 0.5 .*(1.0 .+ z[1][i1])
		mz2   = 0.5 .*(1.0 .- z[2][i2]) ; pz2 = 0.5 .*(1.0 .+ z[2][i2])
		if n == 4
			xy = [v[1][j] .* mz1 .* mz2 .+ v[2][j] .* pz1 .* mz2 .+
				  v[4][j] .* mz1 .* pz2 .+ v[3][j] .* pz1 .* pz2 for j = 1:d]
			# JACOBIAN
			dz1 = [mz2 .* (v[2][j]-v[1][j]) .+ pz2 .* (v[3][j]-v[4][j]) for j = 1:d]
			dz2 = [mz1 .* (v[4][j]-v[1][j]) .+ pz1 .* (v[3][j]-v[2][j]) for j = 1:d]
			return xy, 0.25 .*(dz1[1] .* dz2[2] .- dz2[1] .* dz1[2])
		else
			xy = [v[1][j] .*(mz1 .* mz2).+v[2][j].*(pz1 .* mz2).+ v[3][j].* pz2 for j = 1:d]
			j1 = [-0.5 .* v[1][j] .* mz2 .+ 0.5 .* v[2][j] .* mz2 for j = 1:2]
			j2 = [-0.5 .* v[1][j] .* mz1 .- 0.5 .* v[2][j] .* pz1 .+ 0.5 .* v[3][j] for j = 1:2]
			return xy, (j1[1] .* j2[2] .- j2[1] .* j1[2])

			# [-1,1]^2 to triangular region
			xi1 = 0.5 .*(1.0 .+ z[1][i1]).*(1.0 .- z[2][i2]).-1.0
			xi2 = z[2][i2]
			#xi1, xi1 represent the collapsed axis
			xy  = 0.5 .* [v[1][j].*(.-xi2 .-xi1).+v[2][j].*(1.0.+xi1).+ v[3][j].*(1.0.+xi2) for j = 1:d]
			j1  = 0.5 .* (v[2]-v[1]) ; j2  = 0.5 .* (v[3]-v[1]) ;
			return xy, 0.5 .* (1.0 .-z[2][i2]).*(j1[1] .* j2[2] .- j2[1] .* j1[2])
		end
	else
		i1,i2,i3 = tensor_index(length(z[1]),length(z[2]), n3=length(z[3]))
		if n == 8
			mz1 = 0.5 .*(1.0 .- z[1]) ; pz1 = 0.5 .*(1.0 .+ z[1])
			mz2 = 0.5 .*(1.0 .- z[2]) ; pz2 = 0.5 .*(1.0 .+ z[2])
			mz3 = 0.5 .*(1.0 .- z[3]) ; pz3 = 0.5 .*(1.0 .+ z[3])
			xyz = [v[1][j].*mz1[i1].*mz2[i2].*mz3[i3] .+ v[2][j].*pz1[i1].*mz2[i2].*mz3[i3] .+
			       v[3][j].*pz1[i1].*pz2[i2].*mz3[i3] .+ v[4][j].*mz1[i1].*pz2[i2].*mz3[i3] .+
				   v[5][j].*mz1[i1].*mz2[i2].*pz3[i3] .+ v[6][j].*pz1[i1].*mz2[i2].*pz3[i3] .+
				   v[7][j].*pz1[i1].*pz2[i2].*pz3[i3] .+ v[8][j].*mz1[i1].*pz2[i2].*pz3[i3] for j = 1 : d]
		else
			# collapse coordinates
			xi = hex_to_tet([z[1][i1],z[2][i2],z[3][i3]])
			# map coordinates [-1,1] ---> [0,1]
			x01 = 0.5 .* (xi[1] .+ 1) ;
			y01 = 0.5 .* (xi[2] .+ 1) ;
			z01 = 0.5 .* (xi[3] .+ 1)
			# get the local coordinates
			xyz = [x01 .* (v[2][j] - v[1][j]) .+
		 		   y01 .* (v[3][j] - v[1][j]) .+
				   z01 .* (v[4][j] - v[1][j]) .+ v[1][j] for j = 1 : d]
		end
		return xyz, nothing
	end
end


"""
Given a value in reference coordinates, computes the physical map
Params:
-------
z: quadrature points ( rectangular )
v: element vertices
Returns:
--------
Physical coordinates
"""
function local_coords(z::vecD, v::Vector{Vector{dbl}})
	d  = length(v[1])
	mz = 0.5 * (1.0 .- z) ; pz = 0.5 *(1.0 .+ z)
	lv = length(v)
	if lv == 4
		return (v[1]*mz[1]*mz[2] + v[2]*pz[1]*mz[2] + v[4]*mz[1]*pz[2] + v[3]*pz[1]*pz[2])[1:d]
	elseif lv == 3
		return (v[1]*mz[1]*mz[2] + v[2]*pz[1]*mz[2] + v[3]*pz[2])[1:d]
	elseif lv == 8
		return (v[1]*mz[1]*mz[2]*mz[3] + v[2]*pz[1]*mz[2]*mz[3] + v[3]*pz[1]*pz[2]*mz[3] + v[4]*mz[1]*pz[2]*mz[3] +
	            v[5]*mz[1]*mz[2]*pz[3] + v[6]*pz[1]*mz[2]*pz[3] + v[7]*pz[1]*pz[2]*pz[3] + v[8]*mz[1]*pz[2]*pz[3])[1:d]
	else
		return (z[1] * (v[2] - v[1]) + z[2] * (v[3] - v[1]) + z[3] * (v[4] - v[1]) + v[1])[1:d]
	end
end

local_coords_by_element(data::MESH, z::vecD, e::Int) = local_coords(z, element_coords(data,ele))

#= this functions are used to evaluate analytic
functions: gets the corresponding parametric value
=#
function global_to_local(z::vecD, v::Vector{Vector{dbl}})
	d     = length(v[1]) ; type = length(v)
	if type == 8
		i1,i2,i3 = tensor_index(length(z1),length(z2),n3=length(z3))
		h        = [ norm(v[2] - v[1]), norm(v[4] - v[1]) ,norm(v[5] - v[1]) ]  # 1 hori 2 vert
		t        = zeros(3, length(i1))
		t[1,:]   = v[1][1].+ 0.5 .* (z1[i1] .+ 1.0) .*h[1]
		t[2,:]   = v[1][2].+ 0.5 .* (z2[i2] .+ 1.0) .*h[2]
		t[3,:]   = v[1][3].+ 0.5 .* (z3[i3] .+ 1.0) .*h[3]
		return t
	else
		i1,i2 = tensor_index(length(z1),length(z2))
		l2    = 0.5 .* (z2[i2] .+ 1.0)
		h     = [ norm(v[2] - v[1]), norm(v[type] - v[1]) ]  # 1 hori 2 vert
		t     = zeros(2, length(i1))
		if type == 4
			l1 = 0.5 .* (z1[i1] .+ 1.0) ;
		else
			aux = 0.5 .*(1.0 .+ z1[i1]).*(1.0 .- z2[i2]) .- 1.0
			l1  = 0.5 .* (aux .+ 1.0) ;
		end
		t[1,:] = v[1][1].+l1.*h[1]
		t[2,:] = v[1][2].+l2.*h[2]
		return t
	end
end

function mesh_global_to_local(data::MESH, zw::Vector{Vector{dbl}})
	q  = data.dim == 2 ? length(zw[1]) * length(zw[2]) : length(zw[1]) * length(zw[2]) * length(zw[3])
	t  = zeros(data.dim,data.N * q)
	je = 1
	for e = 1:data.N
		aux,dumb = local_region(zw,element_coords(data,e))
		t[1,je:je+q-1] = aux[1]
		t[2,je:je+q-1] = aux[2]
		data.dim == 3 && (t[3,je:je+q-1] = aux[3])
        je += q
	end
	data.dim == 2 ? (return [t[1,:],t[2,:]]) : (return [t[1,:],t[2,:],t[3,:]] )
end

"""
returns global coordinates local element --> standard element.
For triangles, returns rectangular coordinates [eta1,eta2] = [-1,1]^2 (uncollapsed) "
"""
function global_coords(p::Vector,v::Vector{Vector{dbl}})
	n = length(v)
	if n == 8
		h = [norm(v[2]-v[1]) , norm(v[4]-v[1]) , norm(v[5]-v[1]) ]
		return (2.0 ./ h) .* (p .- v[1]) .- 1
	elseif n == 2
		return norm(p .- v[1]) / norm(v[2]- v[1]) # local coords in edge !  #WARNING THIS IS ASSUMING STRAIGHT ELEMENTS )
	else
		dv1 = v[2]- v[1] ;	dv2 = v[n]- v[1]
		np   = cross(dv1,dv2)
		o1   = cross(np,dv1); o2 = cross(np,dv2)
		z    = [2.0 * dot(p - v[1],o2) / dot(dv1,o2).-1.0,
		        2.0 * dot(p - v[1],o1) / dot(dv2,o1).-1.0]
		n == 4 ? (return z) : abs(1.0 -z[2]) < 1.e-14 ?
							  (return [z[1],z[2]]) : (return [(1 + z[1]) / (1.0 - z[2]) * 2.0 - 1.0, z[2]])
	end
end


## DEBUG FUNCTION
