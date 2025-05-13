
using OrdinaryDiffEq
using Trixi, MSIAC

###############################################################################
# semidiscretization of the compressible Euler equations
gamma = 1.4
equations = CompressibleEulerEquations2D(gamma)

"""
    initial_condition_kelvin_helmholtz_instability(x, t, equations::CompressibleEulerEquations2D)

A version of the classical Kelvin-Helmholtz instability based on
- Andrés M. Rueda-Ramírez, Gregor J. Gassner (2021)
  A Subcell Finite Volume Positivity-Preserving Limiter for DGSEM Discretizations
  of the Euler Equations
  [arXiv: 2102.06017](https://arxiv.org/abs/2102.06017)
"""
function initial_condition_kelvin_helmholtz_instability(x, t, equations::CompressibleEulerEquations2D)
  # change discontinuity to tanh
  # typical resolution 128^2, 256^2
  # domain size is [-1,+1]^2
  slope = 15
  amplitude = 0.02
  B = tanh(slope * x[2] + 7.5) - tanh(slope * x[2] - 7.5)
  rho = 0.5 + 0.75 * B
  v1 = 0.5 * (B - 1)
  v2 = 0.1 * sin(2 * pi * x[1])
  p = 1.0
  return prim2cons(SVector(rho, v1, v2, p), equations)
end

initial_condition = initial_condition_kelvin_helmholtz_instability

surface_flux = flux_lax_friedrichs
volume_flux  = flux_ranocha
polydeg = 2
basis = LobattoLegendreBasis(polydeg)
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max=0.002,
                                         alpha_min=0.0001,
                                         alpha_smooth=true,
                                         variable=density_pressure)

volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)


mesh_file = joinpath(@__DIR__, "reference_box.mesh")

mesh = UnstructuredMesh2D(mesh_file, periodicity=true)


semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0,  3.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=20,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

stepsize_callback = StepsizeCallback(cfl=1.3)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)


###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary

solL = deepcopy(sol)
solT = deepcopy(sol)

using Plots

function dump_SIAC_files(mFile::String, msh::Any, sFile::String, sol::Matrix{Float64}, q::Vector{Int}, eQT::Vector{String})
	# THE MESH
	aux = joinpath(@__DIR__, mFile)
	open(aux, "w") do io
		println(io,"\$Nodes")
		println(io,msh.n_corners)
		for n = 1:msh.n_corners
			println(io,n," ", msh.corners[1,n]," ",msh.corners[2,n]," ",msh.corners[3,n])
		end
		println(io,"\$EndNodes")
		## ELEMENT MAP
		println(io,"\$Elements")
		for k = 1 : msh.n_elements
			ele = msh.element_node_ids[:,k]
			println(io,k," ",length(ele))
			for i in ele	println(io," ",i) ;
			end
		end
		println(io,"\$EndElements")
		#For now, I am ignoring boundaries
		@info " Written in $name"
	end
	@info " Written in mesh file $aux"
	aux = joinpath(@__DIR__,sFile)
	eType = length(msh.element_node_ids[:,1]) #WARNING assuming all elements one kind for now
	open(aux, "w") do io
			println(io,"\$Fields")
			if qmesh.dim  == 2
				println(io,msh.n_elements," ",q[1]," ",eQT[1]," ",q[2]," ",eQT[2])
			else
				println(io,msh.n_elements," ",q[1]," ",eQT[1]," ",q[2]," ",eQT[2]," ",q[3]," ",eQT[3])
			end

			for j = 1:msh.n_elements
				print(io,j," ",eType," ")
				f = sol[:,j]
				[ print(io,f[k]," ") for k = 1:length(f)]
				print(io,"\n")
			end
			println(io,"\$EndFields")
	end
	@info " Written in field file $aux"
end


msh   = semi.mesh
N     = msh.n_elements
qmesh = MSIAC.load_trixi_mesh(msh)
totFun = Trixi.wrap_array_native(sol.u[2], semi)

q = [polydeg+1,polydeg+1] ; eQT = ["legendre","legendre"]
zw = [MSIAC.GZW(eQT[1], q[1]), MSIAC.GZW(eQT[2], q[2]) ]

sn  = ["rho","u","v","k"]
fu  = zeros(q[1]*q[2],qmesh.N)
fld = Vector{MSIAC.FIELD_2D}(undef,4)

for s = 1 : length(sn)
	for j = 1 : qmesh.N
		fu[:,j] = hcat(totFun[s,:,:,j]...)
		for k = 1:length(fu[:,1])
			@printf("%1.16f\t", fu[k,j])
		end
		println()
	end
	fld[s] = MSIAC.FIELD_2D("legendre", polydeg,qmesh.type == 3, zw = zw, name=sn[s], f = deepcopy(fu), reverse=true)
end

data = MSIAC.DATAFIELD(fld[1:end], qmesh)

rx = 2*polydeg+1 ; ry = polydeg+1;
lx = 2*polydeg+1 ; ly = polydeg+1;

i1,i2 = MSIAC.tensor_index(q[1],q[2], reverse = true)

for s = 1 : 4
	fL = MSIAC.filter_data(data,"line",f =sn[s],rx = rx, lx = lx)
	fT = MSIAC.filter_data(data,"tensor",f =sn[s],rx = rx, lx = lx,ry = ry, ly = ly )
	i  = s
	for j = 1 : qmesh.N
		for k = 1 : length(i1)
			solL.u[2][i] = fL[k,j]
			solT.u[2][i] = fT[k,j]
			i += 4
		end
	end
	@info " FILTERED QUANTITY $(sn[s]) "
end


Plots.plot(sol)
Plots.plot(solL)
Plots.plot(solT)
