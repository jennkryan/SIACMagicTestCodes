using Libdl
const symqLIB = "libsymq.so"
libpath = joinpath(@__DIR__,"triangle_symq")
C_symqlib = Libdl.find_library(symqLIB, [libpath])
rule_full_size(degree::Int) = ccall((:rule_full_size, C_symqlib),Cint,(Cint,), degree );
triangle_area(v1::Vector, v2::Vector, v3::Vector) = ccall((:triangle_area, C_symqlib),
														  Cdouble, (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},), v1[1:2],v2[1:2],v3[1:2]);

function triasymq(n::Int, v1::Vector, v2::Vector, v3::Vector)
  numnodes = rule_full_size(n)
  rnodes  = zeros(2*numnodes)
  weights = zeros(numnodes)
  ccall((:triasymq, C_symqlib), Cvoid,
		 (Cint, Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Cint),
		  n,v1[1:2],v2[1:2],v3[1:2],rnodes,weights,numnodes)
  return reshape(rnodes,2,numnodes), weights
end

function split_polygon(n::Int)
    piv = Vector{Int}[]
    for j = 1:2:n
        vec = filter( x-> x <= n, j .+ [1,2,3])
        length(vec) < 2 && continue
        push!(piv, vcat(1,vec))
    end
    return piv
end


function support2D_fits!(data::MESH, pC::BP, dir::vecD, T::vecD;
	                   breaks = BP[]::Vector{BP})

	(length(T) == 0) && (return BP[], 0.0)
	ndir      = dir ./ norm(dir)
	tDir      = sign(T[end]) .* ndir
	largeStep = max(10.0, abs(T[end])) * max_length(data, pC.e[1]) * tDir
	cN        = pC.z .+ largeStep
	tPN       = 0.0;
	cPN       = (pC.z,cN)
	#@info " CHECKING POINT $pC DIR $dir KNOTS $T "
	# Check if pp point is in element hull: need to choose the appropriate starting element
	vs  = element_coords(data,pC.e[1])
	ne  = data.ELEMENTS[pC.e[1]]
	pcb = pC.e[1]

	# check if it is vertex
	inPol = in_face(vs, pC.z, data.toler)
	eAux  = false
	if inPol[1] == 1  # cN is vertex
		pC.v = ne[inPol[2]]
		pull_vertex!(data, pC, largeStep)
		eAux = true
	elseif inPol[1] == 0 # cN is in edge
		pC.f = sort([ ne[inPol[2]] , ne[mod1(inPol[2]+1, length(ne)) ]])
		pC = pull_edge!(data, pC, dir = largeStep)
		pcb != pC.e[1] && (eAux = true)
	end
	if pC.e[1] == -1
		# we hit a wall. Return and come back with shifted
		pC.e[1] = pcb
		dist, breaks = wall_distance!(data, breaks, pC, cN, tDir)
		return breaks, abs(dist - T[end]) < data.toler ? 0.0 : dist-T[end]
	end

    eAux && (breaks = add_break!(breaks,pC))
	eP   = zeros(3)

	# now start geometry loop: go from t0 -- tN and collect intersections
	#@info " WE START AT $pC  KNOTS $T "
	for knot in T
		cN = pC.z .+ (knot - tPN) .* dir
		#@info " KNOT SEQ $knot next $cN "
		do_while = !(abs(knot) < data.toler && eAux)

		while (do_while)
			# reset point: is a vertex = 0 is in edge = [-1,-1]
			pC.v = 0; pC.f = [-1,-1];
			cPN  = (pC.z,cN)

			# loop around element: find intersection point
			le     = data.ELEMENTS[pC.e[1]]
			lv     = element_coords(data,pC.e[1])
			inPol  = in_face(lv, cN, data.toler)

			eP = deepcopy(cN)
			if inPol[1] == 1  # cN is vertex
				pC.v = le[inPol[2]]
			elseif inPol[1] == 0 # cN is in edge
				pC.f = sort([ le[inPol[2]] , le[mod1(inPol[2]+1, length(le)) ]])
			elseif inPol == [-1,0] #cN is strictly oustide polygon. Loop around edges and find intersection
				for j in data.FMAP
					eP = line_intersection(cPN, (lv[j[1]], lv[j[2]]), data.toler)
					isempty(eP) && continue
					(norm(eP .- pC.z) < data.toler) && continue # could happen is previous point was vertex. Avoid getting stuck
					# store point type : either save the vertex  or edge ID
					pC.v = norm(lv[j[1]]- eP) < data.toler ? le[j[1]] :
						   norm(lv[j[2]]- eP) < data.toler ? le[j[2]] : 0
					pC.v == 0 && (pC.f = sort(le[j]))
					break
				end
			end

			nt       = norm(pC.z .- eP)
			tACC     = nt * sign(T[end]) / norm(dir) # kernel coordinates
			pC.t    += nt .* tDir
			tPN     += tACC
			pC.z    .= eP
			breaks   = add_break!(breaks,pC) #add break from current element
			abs(tPN - T[end]) < data.toler && break  #sto pif we are at the end of our loop
			do_while = norm(pC.z - cN) > data.toler

			back = pC.e[1]
			if inPol != [-1,1] # Point cN was at element bounds or outside
				# check that point is actually in element (can be in adjacent for unstructured meshes )
				cNP = norm(pC.z - cN) #save distance in case we apply periodic BCs
				if pC.v > 0
					pull_vertex!(data, pC, largeStep) #vertices have multiple elements. Find the next element according to the search direction
				else
					pC.f == [-1,-1] && throw(@error " GEOMETRY ERROR: UNEXPECTED RELATIVE POSITION !!!!! ")
					pC = pull_edge!(data, pC, dir = largeStep, plotFig = false) #  get adjacent even if it is periodic and update coordinates if so
				end
				cN  .= pC.z + cNP * tDir  # if we have periodic bounds, need to update coordinates from per. vertex
			end

			if pC.e[1] == -1  # we have hit a wall. Return shift for kernel
				pC.e[1] = back
				dist, breaks = wall_distance!(data, breaks, pC,cN,tDir)
				tPN += dist
				return breaks, abs(tPN - T[end]) < data.toler ? 0.0 : tPN-T[end]
			end
			# if point wasn't inside, add the break from the new break. Basically, add the adjacent element ID
			inPol != [-1,1] && add_break!(breaks,pC)
		end
	end
	return breaks, abs(tPN - T[end]) < data.toler ? 0.0 : tPN-T[end]
end


"""
This is a high-level function: computes the convolution for the tensor filter & quad elements
"""
function tensor_convolution(msh::MESH,fld::FIELD_2D, kernels::Matrix{KERNEL}, eID::Int,
	                        zeta::vecD, theta::Union{vecD,String}, scaling::Union{vecD,String}, inexact::Bool, track_footprint::Bool)


	# FRAME : USE BDRY 1 AS "SOUTH" AND BDRY 2 AS "WEST"
	frame = get_frame(msh)
	if theta == "auto" || scaling == "auto"
		vs = element_coords(msh, eID)
		he = element_lengths(msh,eID)
		ca = CircularArray([1:length(he);])
		if msh.structured
			hx =[ abs(ndot(frame[:,1],vs[ca[j+1]]-vs[ca[j]])) for j = 1 : length(he)]
			hy =[ abs(ndot(frame[:,2],vs[ca[j+1]]-vs[ca[j]])) for j = 1 : length(he)]
			scaling = [he[findfirst(x-> isapprox(x,1.0),hx)], he[findfirst(x-> isapprox(x,1.0),hy)],0.0]
		else
			scaling = [maximum(he),maximum(he),0.0]
		end
		kDIR = [frame[:,1] .* scaling[1],frame[:,2] .* scaling[2],0.0]
	else
		kDIR = [[cos(theta[1]),sin(theta[1])] .* scaling[1] , [cos(theta[2]),sin(theta[2])] .* scaling[2],0.0]
	end
	centre = BP(local_coords(zeta, element_coords(msh,eID)),[0.0,0.0,0.0], [eID],0)
	# MAIN AXIS

	kT = [eSYM, eSYM]
	kT[1], mA = get_element_list!(kernels[:,1], msh, centre, kDIR[1], skip_knots=inexact)
	isempty(mA) && (return evaluate_2Dpoint(fld,eID,zeta))

	convolution = 0.0


	if track_footprint
		for k in mA
			@info " MY FIRST AX $(k.z) $(k.t[1:2] ./ scaling[1:2])"
		end
		plot_mesh(msh)
		fig = figure(1)
		af = fig.gca()

		ca = vcat([1:msh.type;],1)
		for j = 1 : msh.N
			coords = hcat(element_coords(msh,j)...)
			af.plot(coords[1,ca], coords[2,ca], c = "indigo")
			cen = element_centre(msh,j)
			af.text(cen[1], cen[2], string(j))
		end
		fig2 = figure(2)
		ax = [fig2.add_subplot(121), fig2.add_subplot(122)]
		for (k,v) in msh.NODES
			af.text(v[1][1], v[1][2], string(k), c = "red")
		end
		for j in mA
			af.scatter(j.z[1], j.z[2])
		end
	end

	nQ = fld.basis.degree + maximum(getfield.(kernels,:l))+2
	kx = kernels[kT[1],1]

	for l = 1:length(mA) - 1
		aux = intersect(mA[l].e, mA[l+1].e)
		isempty(aux) && continue

		kT[2], b1 = get_element_list!(kernels[:,2], msh, BP(mA[l].z, mA[l].t, [aux[1]], mA[l].v), kDIR[2], secondAxis = true , skip_knots=inexact)
		isempty(b1) && (return evaluate_2Dpoint(fld,eID,zeta))


		k2, b2 = get_element_list!(kernels[:,2], msh, BP(mA[l+1].z, mA[l+1].t, [aux[1]], mA[l+1].v), kDIR[2],kT = [kT[2]], secondAxis = true, skip_knots=inexact)
		isempty(b2) && (return evaluate_2Dpoint(fld,eID,zeta))
		if k2 != kT[2]
			throw(@error " Unsupported domain type ")
		end

		[kT[d] == ePOS && (kernels[ePOS,d].c = kernel_coefficients(kernels[ePOS,d].T, kernels[ePOS,d].zw)) for d = 1:2]
		lb1 = length(b1) ; lb2 = length(b2) ;

		c1 = 1 ; c2 = 1;
		ky = kernels[kT[2],2]
		supp_box = kernel_width(ky) * scaling[2]
		if track_footprint
			@info  " LEFT LIST "
			for k in b1
				println("$k")
				if track_footprint
					af.scatter(k.z[1],k.z[2], s = 10, zorder = 5 , color="red")
				end
			end
			@info  " RIGHT LIST "
			for k in b2
				println("$k")
				if track_footprint
					af.scatter(k.z[1],k.z[2], s = 10, zorder = 5 , color="red")
				end
			end
		end

		while (c1 < lb1 || c2 < lb2)

			d1 = min(lb1, c1+1) ; 	d2 = min(lb2, c2+1)
			(d1 == c1 && d2 == c2 ) && break
			norm(b1[d1].z .- b1[c1].z) > supp_box && (d1 = c1)
			norm(b2[d2].z .- b2[c2].z) > supp_box && (d2 = c2)
			if d1 == c1 && d2 == c2
				c1 = min(lb1, c1+1) ; 	c2 = min(lb2, c2+1)
				continue
			end
			## for break points that are vertices, collect all surrounding elements


			ec1 = b1[c1].v >0 ? msh.NODES[b1[c1].v][2] : b1[c1].e
			ed1 = b1[d1].v >0 ? msh.NODES[b1[d1].v][2] : b1[d1].e
			ec2 = b2[c2].v >0 ? msh.NODES[b2[c2].v][2] : b2[c2].e
			ed2 = b2[d2].v >0 ? msh.NODES[b2[d2].v][2] : b2[d2].e

			## try finding common elements to see if we can integrate directly
			i1  = intersect(ec1, ed1) ;	i2  = intersect(ec2, ed2)
			i12 = intersect(i1,i2)
			if isempty(i12) && d1 > c1 && d2 > c2
				if !isempty(intersect(i1, ec2))
					d2 = c2 ;i2 = ec2
				elseif !isempty(intersect(i2, ec1))
					d1 = c1; i1 = ec1
				end
				i12 = intersect(i1,i2)
			end
			#draw kernel box in physical & kernel space
			kBOX = [b1[c1].t,b2[c2].t]
			kXYZ = [b1[c1].z,b2[c2].z]

			if d2 == c2 + 1
				push!(kBOX,b2[d2].t)
				push!(kXYZ,b2[d2].z)
			end
			if d1 == c1 + 1
				push!(kBOX,b1[d1].t)
				push!(kXYZ,b1[d1].z)
			end

			if track_footprint
				for k = 1: size(kBOX,1)
					af.scatter(kXYZ[k][1],kXYZ[k][2], color = "r")
					@info " T $(kBOX[k][1:2] ./ scaling[1:2])    POS $(kXYZ[k])   "
				end
			end

			h12 = BP[]
			## if no intersection, draw a box to find all possible extra cuts with the mesh
			if isempty(i12)
				if isempty(intersect(ec1,ec2))
					p1    = BP(b1[c1].z, b1[c1].t, [b1[c1].e[1]], b1[c1].v)
					#track_footprint && 	@info " KNOT LIST FirsT  $kt   DIRECTION $(p1.z) towards $(b2[c2].z ) "
					h12,dumb = support2D_fits!(msh, p1,b2[c2].z .- p1.z,[0.0,1.0])
				end
				if isempty(intersect(ed1,ed2))
					p2    = BP(b2[d2].z, b2[d2].t, [b2[d2].e[1]], b2[d2].v)
					#track_footprint && 	@info " KNOT LIST SEcOND  DIRECTION $(p2.z) towards $(b1[d1].z ) "
					h12,dumb = support2D_fits!(msh, p2, b1[d1].z .- p2.z,[0.0,1.0],breaks = h12)
				end
			end

			## Now add the bounding box as breaks for the integrals
			h12 = add_break!(h12, b1[c1]) ;
			d1 > c1 && (h12  = add_break!(h12, b1[d1])) ;
			h12 = add_break!(h12, b2[c2])
			(d2 > c2) && (h12  = add_break!(h12, b2[d2])) ;

			vl  = findall(getfield.(h12,:v) .> 0)
			[h12[j].e = msh.NODES[h12[j].v][2] for j in vl ]
			totE = unique!(vcat(getfield.(h12,:e)...))

			if track_footprint
				@info " MY ELEMENT LIST FOR THIS BREAKS  $totE **************************\n\n "
				for k in h12
					@info " BREAKS $k"
				end
			end

			for e in totE
				track_footprint && @info " INTEGRATING ELEMENT $e in $totE"
				# look around the element vertices and check if vertex lives in kernel box
				for k in msh.ELEMENTS[e]
					aux = findfirst(x->x == k, getfield.(h12,:v))
					if aux != nothing
						e ∉ h12[aux].e && push!(h12[aux].e,e)
						continue
					end
					zt = dbl[]
					tris = length(kXYZ) == 4 ? [[1,2,3],[1,3,4]] : [[1,2,3]]
					nXYZ = msh.NODES[k][1]
					for t in tris
						if ccwV(kXYZ[t[1]],kXYZ[t[2]],kXYZ[t[3]]) > 0 &&
							in_face(kXYZ[t],nXYZ, msh.toler)[2] != 0
							z12 = global_coords(nXYZ,kXYZ[t])
							zt  = local_coords(z12, kBOX[t])
							break
						end
					end
					!isempty(zt) && (h12 = vcat(h12,BP(nXYZ, zt, [e], k)))
				end
				# Now collect all breaks corresponding to same element
				eles = findall(x-> e in x, getfield.(h12,:e))
				length(eles) < 3 && continue # not enough to integrate. SKIP !

				vs  = getfield.(h12[eles],:z)
				vID = convex_polygon(vs)
				n0  = size(vs,1)
				vO  = vs[vID]
				tO  = [ [h12[j].t[1], h12[j].t[2], 0.0]  for j in eles[vID] ]

				if track_footprint
					@info " MY BREAK LIST "
					for k = 1:n0
						println("$(tO[k]) SCALED $(tO[k][1:2] ./ scaling[1:2] ) -------------- $(vO[k])")
					end
				end

				# now check that spline breaks are not crossed within the same element. Otherwise, we need to split the integral parts for exact quadratures
				pI  = Int[] ; pT = Vector{dbl}[] ; pV = Vector{dbl}[] ;
				pos = [1:n0;]

				for k = 1 : n0
					kp    = mod1(k+1,n0)
					yLIMS = [ tO[k][2] / scaling[2], tO[kp][2] / scaling[2] ]  #normalize for knot check
					yPIV  = sortperm(yLIMS)
					# Check for internal spline knots. Excludee yLIMS, already accounted for !!
					spB   = unique!(filter(x->( x > yLIMS[yPIV[1]]+msh.toler && x < yLIMS[yPIV[2]]-msh.toler), sort(vcat(ky.T...))))
					isempty(spB) && continue
					pos[k+1:end] .+= length(spB)
					for s = 1:length(spB)
						push!(pI, pos[k] + s)
						t2 = (spB[s] * scaling[2] - tO[k][2]) / (tO[kp][2]-tO[k][2])
						push!(pT,tO[k] .+ t2 .*  (tO[kp] - tO[k]))
						push!(pV,vO[k] .+ t2 .*  (vO[kp] - vO[k]))
					end
				end

				if length(pI) > 0
					# merge all integration points (ordered) into a single Vector
					n  = n0 + length(pI)
					aT = zeros(3,n) ; aV = zeros(3,n)
					for k = 1:n0
						aV[:,pos[k]] = vO[k] ; aT[:,pos[k]] = tO[k] ;
					end
					for k = 1:length(pI)
						aV[:,pI[k]] = pV[k]; aT[:,pI[k]]  = pT[k]
					end
					vO = [ aV[:,k] for k = 1 : n]
					tO = [ aT[:,k] for k = 1 : n]
				end

				if track_footprint
					for j = 1 : size(vO,1)
					@info " $j --> $(tO[j])  $(vO[j])"
					end
				end

				ni      = size(vO,1)
				ca      = CircularArray([1:ni;])
				pAngles = [ndot(vO[ca[j-1]]  .- vO[ca[j]],vO[ca[j+1]] .- vO[ca[j]]) for j = 1:ni]
				ir      = Vector{Int}[]
				if ni == 4 && all( x-> abs(x.- pAngles[1]) .< fTOL, pAngles[2:end]) # squared / rectangular. Directly
					ir = [[1:4;]]
				elseif ni == 3
					ir = [[1:3;]]
				else
					tf    = hcat(tO...)
					yLIMS = [minimum(tf[2,:])/scaling[2], maximum(tf[2,:])/scaling[2]]
					spB   = unique!(filter(x->( x >= yLIMS[1]-msh.toler&& x <= yLIMS[2]+msh.toler), sort(vcat(ky.T...))))
					ir1 = split_polygon(ni)
					if  isempty(spB)
						ir = split_polygon(ni) # split the polygon and integrate
					else
						tmin = yLIMS[1] .* scaling[2]
						bl = unique!( vcat(spB,yLIMS[2]).* scaling[2])
						for j in bl
							aux = findall(x -> x >= tmin-msh.toler && x <= j+msh.toler, tf[2,:])
							idx = split_polygon(length(aux))
							[push!(ir, aux[k]) for k in idx]
							tmin = j
						end
					end
				end
				track_footprint && @info " INTEGRAL GROUPS $ir $n0 $ni == $(length(h12)) ? "
				for pv in ir  # triangulate polygon

					nv = length(pv)
					(nv != 3 && nv != 4) && throw(@error " NOT AN INTEGRAL ")
					ve = element_coords(msh,e)

					# this lives in rectangular coordinates (xi for quads, etas for tris )
					sR = [ global_coords(vO[j],ve) for j in pv]
					kR = [ [tO[j][1] / scaling[1], tO[j][2] / scaling[2] ]   for j in pv]

					if nv == 3
						z,w   = triasymq(nQ, kR[1], kR[2], kR[3])
						kXY   = [z[1,:], z[2,:]]
						z,wR  = triasymq(nQ, sR[1], sR[2], sR[3])
						sXY   = [z[1,:], z[2,:]]
						_,jac = triasymq(nQ, vO[pv[1]], vO[pv[2]], vO[pv[3]])
					else
						zw     = GZW("legendre",nQ)
						kXY,j2 = local_region([zw.nodes, zw.nodes],kR)
						sXY,j1 = local_region([zw.nodes, zw.nodes],sR)      # coordinates for evaluating function
						_,jac  = local_region([zw.nodes, zw.nodes],vO[pv])  # jacobian info
						ix,iy  = tensor_index(nQ,nQ)
						jac  .*= zw.weights[ix] .* zw.weights[iy]
					end

					if track_footprint
						piv = vcat([1:length(kR);],1)
						lx  = [kR[j][1] for j in piv]
						ly  = [kR[j][2] for j in piv]

						# kernel coordinates
						ax[1].fill(lx,ly, edgecolor="white", lw = 0.2, alpha= 0.2)
						aux = hcat(kXY...)
						ax[1].scatter(kXY[1], kXY[2], c = "k", s = 0.5)
						tX = unique!(vcat(kx.T...))
						tY = unique!(vcat(ky.T...))
						xlims = [minimum(tX), maximum(tX)]
						ylims = [minimum(tY), maximum(tY)]
						for k in tX
							ax[1].plot([k,k],ylims,  c = "orange", ls = "--", zorder = 3)
						end
						for k in tY
							ax[1].plot(xlims, [k,k],  c = "orange", ls = "--", zorder = 3)
						end
					end


					kval = [evaluate_kernel(kx.T,kx.c,kXY[1]),evaluate_kernel(ky.T, ky.c, kXY[2])]
					if track_footprint
						poly = hcat(vO[pv]...)
						cnv  = vcat([1:length(pv);],1)
						af.fill(poly[1,cnv], poly[2,cnv], edgecolor="white", lw = 0.2, alpha = 0.5)
						# kernel coordinates
						ax[2].plot(kXY[2], kval[2])
						tf = hcat(tO[pv]...)
						ax[1].plot(tf[1,:] /scaling[1], tf[2,:]/scaling[2], c = "indigo", zorder = 4, lw = 0.4)
						ax[2].scatter(kXY[2], kval[2], marker=".")
						kerT = unique!(vcat(ky.T...))
						for k in kerT
							ax[1].plot( [minimum(vcat(kx.T...)), maximum(vcat(kx.T...))], [k,k], c = "red", zorder = 4, lw = 0.4, ls = "--")
						end
					end

					pLeg = eval_basis(fld.basis,sXY[1], sXY[2])
					proj = [sum(fld.modes[:, e] .* pLeg[i,:]) for i = 1:length(jac) ]
					convolution += sum(kval[1] .* kval[2] .* jac .* proj)
				end
				n0 == length(h12) && break # already done loop around any more breaks in h. all have been accounted for
			end
			c1  = d1; c2 = d2
		end
	end
	convolution /= prod(scaling[1:2])
	#@info " CONVOLUTION $(convolution - 1.0) "
	#@info " CONVOLUTION $convolution center $centre"
	return convolution
end

"""
This is a high-level function: computes the convolution for the line filter for quad & tris elements
"""
function line_convolution(msh::MESH,fld::FIELD_2D, kernel::Vector{KERNEL}, eID::Int, zeta::vecD,
	                      theta::Union{dbl,String},scaling::Union{dbl,String}, track_footprint::Bool)


	zw = GZW("legendre", fld.basis.degree+1+maximum(getfield.(kernel,:l)))
	z  = zw.nodes ; w = zw.weights
	#@info " POST PROCESSING $eID ZETA $zeta "
	# MAIN AXIS
	he    = element_lengths(msh, eID)
	nodes = element_coords(msh, eID)
	if theta == "auto" || scaling == "auto"
        if msh.type == 3
    		ca    = [1,2,3,1]
			vj    = [nodes[ca[j+1]] .- nodes[ca[j]] for j = 1:3]
			va    = [acos(clamp( dot(vj[ca[j+1]],vj[ca[j]])/(he[ca[j+1]] * he[ca[j]]), -1.0, 1.0)) for j = 1:3]
            (theta == "auto") && (theta = minimum(vcat(va, pi .- va)))
			(scaling == "auto") && (scaling = maximum(he))
        else msh.type == 4
			if scaling == "auto"
                scaling = theta == 0.0 ? he[1] : theta == 0.5 * pi ? he[2] : sqrt(he[1]^2+he[2]^2)
			end
			theta == "auto" && (theta = atan(he[2]/he[1]))
        end
    end
	kDIR   = [cos(theta), sin(theta), 0.0] .* scaling
	centre = BP(local_coords(zeta, nodes),[0.0,0.0,0.0], [eID],0)
	kT,kB  = get_element_list!(kernel, msh, centre, kDIR)
	if isempty(kB)
		theta = pi - theta
		kDIR = [cos(theta), sin(theta),0.0] .* scaling
		kT,kB = get_element_list!(kernel, msh, centre, kDIR)
		if isempty(kB)
			theta = pi - theta
			scaling *= 0.5
			kDIR = [cos(theta), sin(theta),0.0] .* scaling
			kT,kB = get_element_list!(kernel, msh, centre, kDIR)
			if isempty(kB)
				kDIR = [cos(theta), sin(theta),0.0] .* scaling
				kT,kB = get_element_list!(kernel, msh, centre, kDIR)
				isempty(kB) && return evaluate_2Dpoint(fld,eID,zeta)
			end
		end
	end

	kT == ePOS && (kernel[ePOS].c = kernel_coefficients(kernel[ePOS].T, kernel[ePOS].zw))

	if track_footprint
		fig = figure(3)
		af = fig.gca()
		for j = 1 : msh.N
			coords = hcat(element_coords(msh,j)...)
			ca     = [1:length(coords[1,:]);1]
			af.plot(coords[1,ca], coords[2,ca], c = "indigo")
			cen    = element_centre(msh,j)
			af.text(cen[1], cen[2], string(j))
		end
		for j = 1 : length(kB)
			af.scatter(kB[j].z[1], kB[j].z[2], c = "r", zorder = 4)
			j == length(kB) && continue
			af.plot([kB[j+1].z[1],kB[j].z[1]],[kB[j+1].z[2],kB[j].z[2]])
		end
		af.scatter(centre.z[1],centre.z[2], c = "b", zorder = 5)
	end

	convolution = 0.0
	for l = 1:length(kB)-1
		### KERNEL COORDINATES t: breaks in kernel space ( t = b/H )
		le = intersect(kB[l].e,kB[l+1].e)
		isempty(le) && continue

		t2 = norm(kB[l+1].t) / scaling * sign(dot(kDIR,kB[l+1].t));
		t1 = norm(kB[l].t)   / scaling * sign(dot(kDIR,kB[l].t));
		dk = 0.5 * (t2-t1)#nt   # parametric space
		# Map for gauss quadrature
		kerT  = dk.* z .+ 0.5 * (t2 +t1)

		kval  = evaluate_kernel(kernel[kT].T,kernel[kT].c, kerT)
		### DG COORDINATES z (phys): map to local (x) -> chi -> eta
		lxy   = [0.5.*((kB[l+1].z - kB[l].z) .* j .+ kB[l+1].z .+ kB[l].z) for j in z]
		ve    = element_coords(msh,le[1])
		legXY = hcat([global_coords(lxy[j],ve) for j = 1:length(z)]...)
		pLeg  = eval_basis(fld.basis,legXY[1,:], legXY[2,:])
		proj  = [sum(fld.modes[:, le[1]] .* pLeg[i,:]) for i = 1:length(z)]
		convolution  += sum(proj .* kval .* w) * dk
	end
	return  convolution
end
