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

using Test, Printf, MSIAC
@info " Running Julia with $(Threads.nthreads()) threads "
### Generate mesh & field files ####
using PyCall, PyPlot
pygui(true)

nX = 7 ; nY = 7 ; nZ = 7
dX = [0.0,2.0]
dY = [0.0,2.0]
dZ = [0.0,2.0]
msh = joinpath(@__DIR__,"Tets.msh")
if !isfile(msh)
	fm  = MSIAC.gmsh_3Dmesh("tet" , nX, nY, nZ, dX, dY, dZ,joinpath(@__DIR__,"Tets"), true)
end

sol  = "poly_xyz"
q    = 3
ff   = MSIAC.create_analytic_field([q,q,q],["legendre","legendre","legendre"],msh,joinpath(@__DIR__,"field.txt"),type=sol)
data = MSIAC.load_data(msh,ff, structured = true)


fig = plt.figure(figsize=(10, 10))
ax  = fig.add_subplot(111,projection="3d")

vs  = MSIAC.element_coords(data.mesh,1)
vsf = hcat(vs...)

zw  = repeat([MSIAC.GZW("legendre",6)] , 3)
xyz,_ = MSIAC.local_region(getfield.(zw,:nodes),vs)
pts = hcat(xyz...)
ax.scatter(pts[:,1],pts[:,2],pts[:,3], c = "r")

for j = 1 : length(vs)
	ax.text(vs[j][1],vs[j][2],vs[j][3],string(j), fontsize=12)
end
for j in data.mesh.FMAP
	fl = vcat(j,j[1])
	ax.plot(vsf[1,fl],vsf[2,fl],vsf[3,fl])
end

show(block=true)
