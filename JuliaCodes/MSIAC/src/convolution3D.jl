function plot_local_map(data::MESH, pC::BP, cN::vecD)
	fig = figure(3)
	ax  = fig.add_subplot(111,projection="3d")
	eles = [pC.e[1]]
	for j in data.FMAP
		push!(eles, adjacent(data, eles[1], j) )
	end
	for j in eles
		j == -1 && continue
		coords = hcat(element_coords(data,j)...)
		nID = data.ELEMENTS[j]
		for k in length(nID)
			ax.text(coords[1,k], coords[2,k],coords[3,k], string(nID[k]), c = "indigo")
		end

		cen    = element_centre(data,j)
		ax.text(cen[1], cen[2],cen[3], string(j))
		for f in data.FMAP
			ca = vcat(f,f[1])
			ax.plot3D(coords[1,ca], coords[2,ca],coords[3,ca], c = "indigo")
		end
	end
	ax.plot3D([pC.z[1],cN[1]],[pC.z[2],cN[2]],[pC.z[3],cN[3]])
	ax.scatter(pC.z[1],pC.z[2],pC.z[3],color="red")
	ax.scatter(cN[1],cN[2],cN[3],color="blue")
end

function support3D_fits!(data::MESH, pC::BP, dir::vecD, T::vecD; breaks = BP[]::Vector{BP}, plotStuff = false)

	(length(T) == 0) && (return BP[], 0.0)

	ndir      = dir ./ norm(dir)
	tDir      = sign(T[end]) .* ndir
	largeStep = max(10.0, abs(T[end])) * max_length(data, pC.e[1]) * tDir
	cN        = pC.z .+ largeStep
	tPN       = 0.0;
	cPN       = (pC.z,cN)

	# Check if pp point is in element hull: need to choose the appropriate starting element
	#@info "\n\n\n Begin support 3d with point $pC direction $dir KNOTS $T "
	vs  = element_coords(data,pC.e[1])
	ne  = data.ELEMENTS[pC.e[1]]
	pcb = pC.e[1]
	inVol = in_volume(vs, data.FMAP,pC.z,data.toler)
	eAux = false
	if inVol != [-1,1]
		face = inVol[1]
		if inVol[2][1] == 1
			pC.v = ne[face[ inVol[2][2]] ]  #inVol[1] = [a,b,c,d] inVol[2][2] = index j in [1:4]
			pull_vertex!(data, pC, largeStep)
			eAux = true
		elseif inVol[2][1] == 0
			pC.f = sort([ ne[ face[inVol[2][2]]] , ne[ face[mod1(inVol[2][2]+1, length(ne)) ]]])
			pC   = pull_edge!(data, pC, dir = largeStep)
			pcb != pC.e[1] && (eAux = true)
		else
			pC   = pull_face!(data, pC, largeStep, ne[face])
			pcb != pC.e[1] && (eAux = true)
			## PULL FACE
		end
	end
	if pC.e[1] == -1
		# we hit a wall. Return and come back with shifted
		pC.e[1] = pcb
		dist, breaks = wall_distance!(data, breaks, pC, cN, tDir)
		return breaks, abs(dist - T[end]) < data.toler ? 0.0 : dist-T[end]
	end
	eP   = zeros(3)
	eAux && (breaks = add_break!(breaks,pC))

	for knot in T
		cN  = pC.z .+ (knot - tPN) .* dir
		do_while = !(abs(knot) < data.toler && eAux)
		#find all intersections between two spline knots

		while (do_while)

			pC.v = 0; pC.f = [-1,-1]
			cPN = (pC.z,cN)

			# loop around element: find intersection point
			le    = data.ELEMENTS[pC.e[1]]
			lv    = element_coords(data,pC.e[1])

			if plotStuff
				fig = figure(1)
				ax  = fig.add_subplot(111,projection="3d")
				ax.axis("off") ; ax.grid(b=nothing)
				vs = hcat(lv...)
				for j = 1 : length(lv)
					ax.text(lv[j][1],lv[j][2],lv[j][3], string(le[j]))
				end
				for f in data.FMAP
					fl = vcat([1:length(f);],1)
					vs = hcat(lv[f]...)
					ax.plot(vs[1,fl], vs[2,fl], vs[3,fl], c = "r")
					neighbor = adjacent(data, pC.e[1],f)
					neighbor == -1 && continue
					aux = element_coords(data,neighbor)
					vs  = hcat(aux...)
					ae  = data.ELEMENTS[neighbor]
					for j = 1 : length(aux)
						ax.text(aux[j][1],aux[j][2],aux[j][3], string(ae[j]))
					end
					for f in data.FMAP
						fl = vcat([1:length(f);],1)
						vs = hcat(aux[f]...)
						ax.plot(vs[1,fl], vs[2,fl], vs[3,fl], c = "r")
					end
				end
				ax.plot3D([pC.z[1], cN[1]], [pC.z[2], cN[2]], [pC.z[3], cN[3]], lw = 3, c = "k" )
			end

			inVol = in_volume(lv, data.FMAP,cN,data.toler, plotFig = true)

			#@info " POINT SEGMENT $(pC.z) -- $cN  IN $lv in volume $inVol "
			#plot_local_map(data,pC, cN)
			interType = false
			iP = deepcopy(cN)
			if inVol == [-1,0] # pN not in volume
				#@info " POINT PN OUTSIDE VOLUME. LOOK FOR INTERSECTINS  "
				for j in data.FMAP
					res = plane_intersection(lv[j],cPN, data.toler)
					isempty(res) && continue
					norm(res[3] - pC.z) < data.toler && continue
					iP = deepcopy(res[3])
					#@info "INTERSECTION GAVE $res $le   $j  $(res[2]) and $(mod1(res[2]+1, length(j)))  "
					interType = res[1] == 2 ? le[j][res[2]] :
				                res[1] == 1 ? [le[j][res[2]], le[ j[mod1(res[2]+1, length(j))] ]] : le[j]
					break
				end
			elseif inVol == [-1,1]
				interType = true
			elseif inVol[2][1] == 1
				interType = le[ inVol[1][ inVol[2][2] ] ]
			elseif inVol[2][1] == 0
				v1 = inVol[1][ inVol[2][2] ] ; v2 = inVol[1][ mod1(inVol[2][2]+1,length(inVol[1])) ] ;
				interType = [ le[v1], le[v2]]
				edge = [coords(data,interType[1]), coords(data,interType[2])]
			else
				interType = le[inVol[1]]
			end

			#@info " DONE CHECKING FOR POINT $pC --> $cN  intersection piont $iP  type $interType invol $inVol "
			# check that point is actually in neighbor element (could have skipped one element)
			nt    = norm(pC.z .- iP)
			tACC  = nt * sign(T[end]) / norm(dir) # kernel coordinates
			pC.t += nt .* tDir
			tPN  += tACC
			pC.z .= iP
			breaks = add_break!(breaks,pC) #add break from current element
			abs(tPN - T[end]) < data.toler && break
			do_while = norm(pC.z - cN) > data.toler

			cNP = norm(pC.z - cN)
			back = pC.e[1]
			if interType == true
				continue
			elseif interType == false # we have hit a wall. Return shift for kernel
				pC.e[1] = -1
			elseif typeof(interType) == Int
				pC.v = interType
				pull_vertex!(data, pC, largeStep)
			elseif length(interType) == 2
				pC.f .= sort(interType)
				pC = pull_edge!(data, pC, dir=largeStep)
			else
				pC = pull_face!(data, pC, largeStep, interType)
			end
			cN  .= pC.z + cNP * tDir  # if we have periodic bounds, need to update coordinates from per. vertex
			if pC.e[1] == -1  # we have hit a wall. Return shift for kernel
				pC.e[1] = back
				dist, breaks = wall_distance!(data, breaks, pC,cN,tDir)
				tPN += dist
				return breaks, isapprox(tPN,T[end],atol=data.toler) ? 0.0 : tPN-T[end]
			end
			add_break!(breaks,pC)
		end
	end
	return breaks, isapprox(tPN,T[end],atol=data.toler) ? 0.0 : tPN-T[end]
end


"""
This is a high-level function: computes the convolution for the line filter for hexes
"""
function line_3Dconvolution(msh::MESH,fld::FIELD_3D, kernel::Vector{KERNEL}, eID::Int, zeta::vecD,
	                        theta::Union{dbl,String},scaling::Union{dbl,String}, track_footprint::Bool)

	sides = element_lengths(msh, eID)
	nodes = element_coords(msh, eID)

	if theta == "auto" || scaling == "auto"
		hx    = sides[1] ; hy = sides[2] ; hz = sides[9]
		theta == "auto"   && (theta = [atan(hy/hx),acos(hz / sqrt(hx^2+hy^2+hz^2))])
		scaling == "auto" && (scaling = hx / ( cos(theta[1]) * sin(theta[2]) ))
	end
	#@info " KNOTS INIT  $(getfield.(kernel, :T))"
	kDIR   = [cos(theta[1]) * sin(theta[2]), sin(theta[1]) * sin(theta[2]), cos(theta[2])] .* scaling
	centre = BP(local_coords(zeta, nodes),[0.0,0.0,0.0], [eID],0)
	#@info " CENTRE $centre ACTUAL ELEMENT cENTRE $(element_centre(msh, eID)) "
	kT,kB  = get_element_list!(kernel, msh, centre, kDIR)

	if isempty(kB)

		#@info " SUPPORT DIDNT FIT ! Try moving around "
		theta[1] = pi - theta[1]
		kDIR = [cos(theta[1]) * sin(theta[2]), sin(theta[1]) * sin(theta[2]), cos(theta[2])] .* scaling
		kT,kB = get_element_list!(kernel, msh, centre, kDIR)

		if isempty(kB)
			theta[1] = pi - theta[1]
			kDIR = [cos(theta[1]) * sin(theta[2]), sin(theta[1]) * sin(theta[2]), cos(theta[2])] .* scaling
			kT,kB = get_element_list!(kernel, msh, centre, kDIR)
			if isempty(kB)
				theta[2] = 0.5 * pi - theta[2]
				kDIR = [cos(theta[1]) * sin(theta[2]), sin(theta[1]) * sin(theta[2]), cos(theta[2])] .* scaling
				kT,kB = get_element_list!(kernel, msh, centre, kDIR)
				isempty(kB) && return evaluate_3Dpoint(fld,eID,zeta)
			end
		end
	end

	kT == ePOS && (kernel[ePOS].c = kernel_coefficients(kernel[ePOS].T, kernel[ePOS].zw))
	#@info " KNOTS END  $(getfield.(kernel, :T))"

	convolution = 0.0
	zw = GZW("Legendre", fld.basis.degree + maximum(getfield.(kernel,:l)) + 2)
	z = zw.nodes; w = zw.weights
	eles = Int[]
	#=@info " -------------- FOOTPRINT -------------- "
	for l in kB
		@info " BREAK $l "
	end
	=#
	if track_footprint
		pygui(true)
		fig4 = figure(4)
	end
	for l = 1:length(kB)-1
		### KERNEL COORDINATES t: breaks in kernel space ( t = b/H )
		le = intersect(kB[l].e,kB[l+1].e)
		isempty(le) && continue
		eles = vcat(eles,le)
		t2   = norm(kB[l+1].t) / scaling * sign(dot(kDIR,kB[l+1].t));
		t1   = norm(kB[l].t)   / scaling * sign(dot(kDIR,kB[l].t));
		dk   = 0.5 * (t2 - t1)#nt   # parametric space
		# Map for gauss quadrature
		kerT  = dk.* z .+ 0.5 * (t2 + t1)
		#@info " KERT $kerT coeffs $(kernel[kT].c)  KNOTS $(kernel[kT].T)"
		kval  = evaluate_kernel(kernel[kT].T,kernel[kT].c, kerT)
		if track_footprint
			axk = fig4.gca()
			axk.plot(kerT, kval)
		end
		### DG COORDINATES z (phys): map to local (x) -> chi -> eta
		lxy    = [0.5.*((kB[l+1].z - kB[l].z) .* j .+ kB[l+1].z .+ kB[l].z) for j in z]
		ve     = element_coords(msh,le[1])
		legXYZ = hcat([global_coords(lxy[j],ve) for j = 1:length(z)]...)
		pLeg   = eval_basis(fld.basis,legXYZ[1,:], legXYZ[2,:], legXYZ[3,:])

		#@info " LEGENDRE BASIS $pLeg  POINTS $legXYZ  modes $(fld.modes[:, le[1]])"
		proj  = [sum(fld.modes[:, le[1]] .* pLeg[i,:]) for i = 1:length(z)]
		convolution  += sum(proj .* kval .* w) * dk
	end
	#@info " CONVOLUTION $convolution"
	if track_footprint
		pygui(true)
		fig = figure(3)
		ax  = fig.add_subplot(111,projection="3d")
		for j in eles
			coords = hcat(element_coords(msh,j)...)
			cen    = element_centre(msh,j)
			ax.text(cen[1], cen[2],cen[3], string(j))
			for f in msh.FMAP
				ca = vcat(f,f[1])
				ax.plot3D(coords[1,ca], coords[2,ca],coords[3,ca], c = "indigo")
			end
		end
		for j = 1 : length(kB)
			ax.scatter(kB[j].z[1], kB[j].z[2], kB[j].z[3], c = "r", zorder = 4)
			j == length(kB) && continue
			ax.plot3D([kB[j+1].z[1],kB[j].z[1]],[kB[j+1].z[2],kB[j].z[2]],[kB[j+1].z[3],kB[j].z[3]])
		end
		ax.scatter(centre.z[1],centre.z[2],centre.z[3], c = "b", zorder = 5)

		show(block=true)
	end
	#@info " **************************************** CONVOLUTION $convolution"
	return  convolution
end
