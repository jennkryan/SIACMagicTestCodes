using Libdl
using PyCall, PyPlot
pygui(true)

const symqLIB = "libsymq.so"
libpath = joinpath(@__DIR__,"triangle_symq")
C_symqlib = Libdl.find_library(symqLIB, [libpath])
rule_full_size(degree::Int) = ccall((:rule_full_size, C_symqlib),Cint,(Cint,), degree );

triangle_area(v1::Vector, v2::Vector, v3::Vector) = ccall((:triangle_area, C_symqlib),
														  Cdouble, (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},), v1,v2,v3 );

function triasymq(n::Int, v1::Vector, v2::Vector, v3::Vector)
  numnodes = rule_full_size(n)
  rnodes  = zeros(2*numnodes)
  weights = zeros(numnodes)
  ccall((:triasymq, C_symqlib), Cvoid,
		 (Cint, Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Cint),
		  n,v1,v2,v3,rnodes,weights,numnodes)
  return reshape(rnodes,2,Int.(numnodes)), weights
end

r8vec_sum(a::Vector) = ccall((:r8vec_sum, C_symqlib),Cdouble,(Cint,Ptr{Cdouble}), length(a), a);

degree = 3
vert1 = [1.0,0.0];
vert2 = [4.0,4.0];
vert3 = [0.0,3.0];
header = "user08";

v = hcat(vert1,vert2,vert3)
c3 = [1,2,3,1]

z,w = triasymq(degree,vert1,vert2,vert3)

a1 = triangle_area(vert1,vert2,vert3)
a2 = r8vec_sum(w)
plot(v[1,c3], v[2,c3])
scatter(z[1,:], z[2,:])
