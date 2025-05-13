#=
*
*      MSIAC: unit testing
*
*      Copyright 2022, Licensed under The GNU Lesser General Public License, version 2.1
*      See http://www.opensource.org/licenses/lgpl-2.1.php
*      Written by: Julia Docampo Sanchez
*      Email: xulia.docampo@pm.me
*
=#

using Test, Printf, MSIAC, Distributed
@info " Running Julia with $(Threads.nthreads()) threads "
### Generate mesh & field files ####

function create_meshes()
	nX = 15 ; nY = 20

	fq = [MSIAC.create_2D_mesh(4,nX,nY,[0.0,1.0],[0.0,1.0],false,joinpath(@__DIR__,"Qmesh1.txt")),
		  MSIAC.create_2D_mesh(4,nX,nY,[0.0,2.0],[0.0,2.0],false,joinpath(@__DIR__,"Qmesh2.txt"),pert = [0.1,0.15], structured = true),
		  MSIAC.create_2D_mesh(4,nX,nY,[0.0,2.0],[0.0,2.0],true, joinpath(@__DIR__,"Qmesh3.txt"))
		 ]

 	st = [true,false,true]

	nX = 12 ; nY = 12
	fq = [MSIAC.create_2D_mesh(3,nX,nY,[0.0,2.0],[0.0,2.0],false,joinpath(@__DIR__,"Tmesh1.txt"), pert = [0.0,0.0],structured = st[1]),
	      MSIAC.create_2D_mesh(3,nX,nY,[0.0,2.0],[0.0,2.0],false,joinpath(@__DIR__,"Tmesh2.txt"), pert = [0.2,0.2],structured = st[2]),
	      MSIAC.create_2D_mesh(3,20,20,[0.0,2.0],[0.0,2.0],true, joinpath(@__DIR__,"Tmesh3.txt"),                 structured = st[3])
		 ]

	nX  = 10 ; nY = 10 ; nZ = 10
 	dX  = [0.0,2.0]
 	dY  = [0.0,2.0]
 	dZ  = [0.0,2.0]
	# GMSH FILES
	fm  = MSIAC.gmsh_3Dmesh("hex" , nX, nY, nZ, dX, dY, dZ,joinpath(@__DIR__,"Hmesh1"))

	return
end

function test_quads()
	@info " Testing filters with quad meshes, BCs: wall & periodic"
	ps = [3,4,2]
	for k = 1 : 3
		msh = joinpath(@__DIR__,"Qmesh"*string(k)*".txt")
		if !isfile(msh)
			create_meshes()
		end
		sol = k == 1 ? "poly_xy" : k == 2 ? "poly_x3y3" : "sincos"
		ff   = MSIAC.create_analytic_field([ps[k],ps[k]],["legendre","legendre"],msh,joinpath(@__DIR__,"field.txt"),type=sol)
		data = MSIAC.load_data(msh,ff, structured = true, modalExp = "legendre")
		p    = MSIAC.polynomial_degree(data)
		Q1   = 1 ; Q2 = 2
		MSIAC.update_phys!(data, ["legendre","legendre"],[Q1,Q2])
		p1,p2 = MSIAC.get_evaluation_points(data)
		exact = MSIAC.get_exact_data(data, sol, [p1,p2])
		fA    = MSIAC.get_field(data)

		fL    = MSIAC.filter_data(data,"line")
		fT    = MSIAC.filter_data(data,"tensor")
		eL    = maximum(abs.(fL .- exact))
		eT    = maximum(abs.(fT .- exact))
		eA    = maximum(abs.(fA .- exact))
		if k < 3
			@test eL < 1.e-10 && eT < 1.e-10 ;
		else
			@test eL < eA && eT < eA ;
		end
	end
	return true
end

function test_tris()

	@info " Testing filters with tri meshes, BCs wall & periodic ! "
	st = [true,false,true]
	ps = [3,1,2]
	for k = 1 : 3
		msh = joinpath(@__DIR__,"Tmesh"*string(k)*".txt")
		if !isfile(msh)
			create_meshes()
		end
		sol = k == 1 ? "poly_xy" : k == 2 ? "ctt" : "sine_xpy"
		ff   = MSIAC.create_analytic_field([ps[k],ps[k]],["legendre","legendre"],msh,joinpath(@__DIR__,"field.txt"),type=sol)
		data = MSIAC.load_data(msh,ff, structured = st[k])
		p    = MSIAC.polynomial_degree(data)
		Q1   = 2 ; Q2 = 2


		MSIAC.update_phys!(data, ["lobatto","lobatto"],[Q1,Q2])
		p1,p2 = MSIAC.get_evaluation_points(data)
		exact = MSIAC.get_exact_data(data, sol, [p1,p2])
		fA    = MSIAC.get_field(data)
		eA    = MSIAC.maximum(abs.(fA .- exact))
		fL    = MSIAC.filter_data(data,"line")
		fT    = MSIAC.filter_data(data,"tensor")#, parallel = false, rx = 3, lx = p+1, ry = 3, ly = p+1)
		eL    = maximum(abs.(fL .- exact))
		eT    = maximum(abs.(fT .- exact))
		eA    = maximum(abs.(fA .- exact))
		if k == 3
			@warn " For Triangles we are not checking the error for the non polynomial approximations"
			return true
		end
		@test eL < 1.e-10 ;#&& eT < 1.e-10 ;
	end
end

function test_hexes()
	@info " Testing line filters with hex meshes, BCs wall  ! "
	msh = joinpath(@__DIR__,"Hmesh1.msh")

	if !isfile(msh)
		create_meshes()
	end

	sol  = "poly_xyz"
	q    = 3
	ff   = MSIAC.create_analytic_field([q,q,q],["legendre","legendre","legendre"],msh,joinpath(@__DIR__,"field.txt"),type=sol)
	data = MSIAC.load_data(msh,ff, structured = true)
	p    = MSIAC.polynomial_degree(data)
	Q    = 1
	qt   = "legendre"
	MSIAC.update_phys!(data, [qt,qt,qt],[Q,Q,Q])
	p123  = MSIAC.get_evaluation_points(data)
	exact = MSIAC.get_exact_data(data, sol, [p123[1], p123[2], p123[3]])

	fA    = MSIAC.get_field(data)
	eA    = MSIAC.maximum(abs.(fA .- exact))
	fL    = MSIAC.filter_data(data,"line")
	eL    = maximum(abs.(fL .- exact))
	eA    = maximum(abs.(fA .- exact))
	@test eL < 1.e-10;
	return true
end

function test_tets()

	@info " Testing line filters with hex meshes, BCs wall & periodic ! "
	msh = joinpath(@__DIR__,"Tets.msh")

	if !isfile(msh)
		create_meshes()
	end

	sol  = "poly_xyz"
	q    = 3
	ff   = MSIAC.create_analytic_field([q,q,q],["legendre","legendre","legendre"],msh,joinpath(@__DIR__,"field.txt"),type=sol)
	data = MSIAC.load_data(msh,ff, structured = true)
	p    = MSIAC.polynomial_degree(data)
	Q    = 1
	qt   = "legendre"
	MSIAC.update_phys!(data, [qt,qt,qt],[Q,Q,Q])
	p123  = MSIAC.get_evaluation_points(data)
	exact = MSIAC.get_exact_data(data, sol, [p123[1], p123[2], p123[3]])
	fA    = MSIAC.get_field(data)
	eA    = MSIAC.maximum(abs.(fA .- exact))
	@info " filter convolution for tets not implemented yet. Only projection is done "
	return true
	fL    = MSIAC.filter_data(data,"line")
	eL    = maximum(abs.(fL .- exact))
	eA    = maximum(abs.(fA .- exact))
end

function test_parallel()
	@info " Testing filters with quad meshes, BCs: wall & periodic"
	ps = [3]
	k = 1
	msh = joinpath(@__DIR__,"Qmesh"*string(k)*".txt")
	if !isfile(msh)
		create_meshes()
	end
	sol = "poly_xy"
	ff   = MSIAC.create_analytic_field([ps[k],ps[k]],["legendre","legendre"],msh,joinpath(@__DIR__,"field.txt"),type=sol)
	data = MSIAC.load_data(msh,ff, structured = true, modalExp = "legendre")
	p    = MSIAC.polynomial_degree(data)
	Q1   = 4 ; Q2 = 4
	MSIAC.update_phys!(data, ["legendre","legendre"],[Q1,Q2])
	p1,p2 = MSIAC.get_evaluation_points(data)
	exact = MSIAC.get_exact_data(data, sol, [p1,p2])
	fA    = MSIAC.get_field(data)

	t1    = @elapsed f1 = MSIAC.filter_data(data,"line")
	t2    = @elapsed f2 = MSIAC.filter_data_hpc(data,"line")

	e1    = maximum(abs.(f1 .- exact))
	e2    = maximum(abs.(f2 .- exact))

	@info " Thread vs Distributed: error  TIMES: $t1  $t2 "
	@test e1 < 1.e-10 && e2 < 1.e-10 ;
	return true
end

function test_gmsh()
	locFdr = @__DIR__
    dX  = [0.0,2.0]
    dY  = [0.0,1.0]
    dZ  = [0.0,3.0]
    nX  = 8 ; nY = 6 ; nZ = 4

    # GMSH FILES
    f2Q = MSIAC.gmsh_2Dmesh("quad", nX, nY, dX, dY,joinpath(locFdr,"gm_quad"))
    f2T = MSIAC.gmsh_2Dmesh("tris", nX, nY, dX, dY,joinpath(locFdr,"gm_tris"))
    f3  = MSIAC.gmsh_3Dmesh("hex" , nX, nY, nZ, dX, dY, dZ,joinpath(locFdr,"gm_hex"))

    mq = MSIAC.load_mesh(f2Q)
    mt = MSIAC.load_mesh(f2T)
    m3 = MSIAC.load_mesh(f3)

    MSIAC.plot_mesh(mq,labels=true,block=false)
	MSIAC.plot_boundary(mq,labels=true,block=false)
    MSIAC.plot_mesh(mt,labels=true,block=false)
	MSIAC.plot_boundary(mt,labels=true,block=false)
    MSIAC.plot_mesh(m3,block=false)
	MSIAC.plot_boundary(m3,labels=true,block=false)
	return true
end

function test_mymesh()
    dX  = [0.0,2.0]
    dY  = [0.0,1.0]
    nX  = 10 ; nY = 15 ;
	locFdr = @__DIR__
    # GMSH FILES
    f2Q = MSIAC.create_2D_mesh(4, nX, nY, dX, dY,true,joinpath(locFdr,"my_quad.txt"),pert = [0.1,0.2], structured = true)
    f2T = MSIAC.create_2D_mesh(3, nX, nY, dX, dY,false,joinpath(locFdr,"my_tri.txt"),pert = [0.2,0.3], structured = false)


    mq = MSIAC.load_mesh(f2Q)
    mt = MSIAC.load_mesh(f2T)


    MSIAC.plot_mesh(mq,labels=true)
    MSIAC.plot_mesh(mt,labels=true)

	return true
end


create_meshes()
t = @elapsed @testset " MSIAC API Unit Testing MESHES " begin
	@test test_gmsh()
	@test test_mymesh()
	@test test_quads()
	@test test_tris()
	@test test_hexes()
	if nworkers() == 1
		Threads.nthreads() == 1 ? addprocs(10) : addprocs(Threads.nthreads())
	end
	@info " TOTAL WORKERS $(nworkers())"
	@everywhere using MSIAC
	test_parallel()
end
