struct POINT
	 x::dbl
	 y::dbl
	 z::dbl
	 POINT(x::Vector) = length(x) == 2 ? new(x[1],x[2],0.0) : new(x[1],x[2],x[3])
	 POINT(x1::dbl,x2::dbl,x3=0.0::dbl) = new(x1,x2,x3)
 end

cross2D(v1::Vector, v2::Vector) = v1[1] * v2[2] - v1[2] * v2[1]
ndot(v1::Vector,v2::Vector) = dot(v1./norm(v1),v2 ./ norm(v2))
const SEG = Union{Tuple{POINT,POINT},Tuple{Vector,Vector}}
Base.:+(p::POINT,q::POINT) = (p.x+q.x,p.y+q.y,p.z+q.z)
Base.:-(p::POINT,q::POINT) = (p.x-q.x,p.y-q.y,p.z-q.z)
Base.:+(p::POINT,q::Vector)  = (p.x+q[1],p.y+q[2],p.z+q[3])
Base.:-(p::POINT,q::Vector)  = (p.x+q[1],p.y+q[1],p.z-q[3])

is_vertex(p::Vector, v::Vector{Vector{dbl}}, tol::dbl) = findfirst(x-> norm(x - p) < tol,v)

in_segment(p::Vector, s::Union{SEG,Vector{Vector{dbl}}},tol::dbl) = abs(norm(p.-s[1]) + norm(s[2].-p) - norm(s[2]-s[1])) < tol


"""
Returns the relative poisition of a point and an element

Parameters:
===========
v  vertex coordinates (ordered!)
p  point to check

Optional:
=========
tol  tolerance deciding if the test passes or not

Return:
=======
a vector [i,j] indicating the position
i = -1 (interior (edge) / exterior) , 0 edge, 1 vertex
j > 0 vertex / edge ID , 1 inside, 0 outside

Examples: [1, 2] point is vertex with id 2
          [-1,1] point is strictly inside face
		  [-1,0] point is strictly outside face
		  [0, 3] point belongs to edge 3 and is neither of its ends
"""
function in_face(v::Vector{Vector{Float64}}, p::Vector,tol::dbl; plotFig=false::Bool)
	# check if p is vertex
	aux = is_vertex(p,v,tol)
	aux != nothing && (return [1,aux])

	nv  = length(v) ;
	res = [-1,-1]
	cp  = 1
	if plotFig
		fig = figure(1)
		fig.add_subplot(projection="3d")
		vs = hcat(v...)
		plot(vs[1,fl], vs[2,fl], vs[3,fl])
		scatter(p[1],p[2],p[3], color = "red")
	end

	normal = cross(v[1]-p, v[2]-p)
	anglesum=0.0
	for j = 1:nv
		in_segment(p, (v[j], v[mod1(j+1,nv)]), tol) && return [0,j] # p is in edge. Save it
		u1 = v[j] - p ; u2 = v[mod1(j+1, length(v))] - p
		nd = ndot(u1,u2)
		cosdot = nd > 1.0 ? 1.0 : nd < -1.0 ? -1.0 : nd
      	anglesum += acos(cosdot)
    end
	return [-1,isapprox(anglesum, 2*pi, atol=1.e-8) ? 1 : 0]
end


in_plane(v::Vector{Vector{Float64}}, pq::SEG,tol::dbl) = isapprox(abs(dot(cross(v[2]-v[1], v[3]-v[1]), pq[2]-pq[1])),0.0, atol = tol)



"""
Returns the relative poisition of a point and an element

Parameters:
===========
v  vertex coordinates (ordered!)
p  point to check

Optional:
=========
tol  tolerance deciding if the test passes or not

Return:
=======
a vector [i,j] indicating the position
i = [a,b,c,d] face to whic it belongs, j -> intersection type in face ( j = [1,1] -> intersection is vertex 1 (see in_face for types))
i = -1 --> j = 1 (inside) 0 (outside)
"""
function in_volume(v::Vector{Vector{Float64}}, flist::Vector{Vector{Int}}, p::Vector, tol::dbl; plotFig=false::Bool)
	# check if p is vertex
	ce = element_centre(v)
	if plotFig == 10
		fig = figure(1)
		ax  = fig.add_subplot(111,projection="3d")
		for f in flist
			fl = vcat([1:length(f);],1)
			vs = hcat(v[f]...)
			ax.plot(vs[1,fl], vs[2,fl], vs[3,fl])
			ax.scatter(p[1],p[2],p[3], color = "red")
		end
		readline()
	end

	inters = 0
	isVert = is_vertex(p,v,tol)
	if isVert != nothing
		idx = findfirst( x-> issubset(isVert, x), flist )
		return [flist[idx], [1,findfirst(x -> x == isVert, flist[idx]) ] ]
	end
	for f in flist
		face = in_face(v[f],p,tol)
		face != [-1,0] && ( [f,face] )
		isempty(plane_intersection(v[f],(ce,p), tol )) && continue
		inters = 1
	end
	return [-1,inters == 0]
end

function line_intersection(ls1::SEG, ls2::SEG, tol::dbl; witness::Bool=false)
	witness = true
	r = ls1[2] - ls1[1]
	s = ls2[2] - ls2[1]

	#a,b = line_intersection2(ls1,ls2,witness = true, tol = tol)
	# first or second line segments is a point

	disjoint = true ; point = dbl[]

    if norm(r) < tol && in_segment(ls1[2],ls2, tol)
		return ls1[2] #witness ? (return false, ls1[2]) : (return false)
	end
	if norm(s) < tol && in_segment(ls2[2],ls1, tol)
		return ls2[2] #witness ? (return false, ls2[2]) : (return false)
    end
	if in_segment(ls1[2], ls2, tol)
		return ls1[2] #witness ? (return false, ls1[2]) : (return false)
	elseif in_segment(ls1[1], ls2, tol)
		return ls1[1] #witness ? (return false, ls1[1]) : (return false)
	end

	c1 = cross(r,s)
	c2 = cross(ls2[1]-ls1[1],s)

	cn1 = c1 ./ norm(c1)
	cn2 = c2 ./ norm(c2)

	par_test = dot(cn1,cn2)
	t = sign(par_test) * norm(c2) / norm(c1)
	if isapprox(abs(par_test),1.0,atol=1.e-8) && in_segment(ls1[1] + t.* r, ls1, tol) && in_segment(ls1[1] + t.* r, ls2, tol)
		disjoint = false ; point =ls1[1] + t .* r
		return ls1[1] + t .* r
	else return dbl[]
	end

	#=if a != disjoint || norm(point-b) > fTOL
		plot([ ls1[1][1],ls1[2][1] ],[ ls1[1][2],ls1[2][2]] )
		plot([ ls2[1][1],ls2[2][1] ],[ ls2[1][2],ls2[2][2]] )
		if !disjoint
			scatter(point[1],point[2])
		end
		throw(@error " I GOT DIFFERENT VALUES LINE $ls1  $ls2  DIS $disjoint $par_test $t  --> $point OLD $a  $b")
	end
	witness ? (return disjoint, point) : (return disjoint)
	=#
end

function edges_disjoint2(ls1::SEG, ls2::SEG; witness::Bool=false, tol=fTOL::dbl)

		r = ls1[2] - ls1[1]
		s = ls2[2] - ls2[1]
	    # first line segment is a point
		if norm(r) < tol
	        empty_intersection = !in_segment(ls1[2],ls2)
	        if witness
	            return (empty_intersection, empty_intersection ? dbl[] : ls1[2])
	        else
	            return empty_intersection
	        end
	    end
		if norm(s) < tol
	        # second line segment is a point
	        empty_intersection = !in_segment(ls2[2],ls1)
	        if witness
	            return (empty_intersection, empty_intersection ? dbl[] : ls2[2])
	        else
	            return empty_intersection
	        end
	    end

	    p1p2 = ls2[1] - ls1[1]
	    u_numerator   = cross2D(p1p2, r)
	    u_denominator = cross2D(r, s)

	    if abs(u_denominator) < tol
			empty_intersection = true
	        # line segments are parallel
	        if abs(u_numerator) < tol
	            # line segments are collinear
				if in_segment(ls1[2],ls2)
	                empty_intersection = false
	                witness && (v = ls1[2]) #: [ls1[1], ls1[2]]

				elseif in_segment(ls2[1],ls1)
	                empty_intersection = false
					witness && (v = ls2[1]) #: [ls1[1], ls1[2]]

				elseif in_segment(ls2[2],ls1)
	                empty_intersection = false
					witness && (v = ls2[2]) #: [ls1[1], ls1[2]]

	            elseif in_segment(ls1[1],ls2)
	                empty_intersection = false
					witness && (v = ls1[1]) #: [ls1[1], ls1[2]]
	            end
			end
	    else
	        # line segments are not parallel
	        u = u_numerator / u_denominator
	        if u < -tol || u > 1.0 + tol
	            empty_intersection = true
	        else
	            t = cross2D(p1p2, s) / u_denominator
	            empty_intersection = t < -tol || t > 1.0 + zeroTOL
	            if witness
	                v = ls1[1] + t * r
	            end
	        end
	    end
	    if witness
	        return (empty_intersection, empty_intersection ? dbl[] : v)
	    else
	        return empty_intersection
	    end
	end

# ** Graham scan algorithm : Sorting points ccw  **
ccw(a::POINT, b::POINT, c::POINT) = ((b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x))
ccwV(a::Vector, b::Vector, c::Vector) = ((b[1] - a[1])*(c[2] - a[2]) - (b[2] - a[2])*(c[1] - a[1]))

function graham_scan!(points::Vector{POINT})
	N = length(points)
	# Place the lowest point at the start of the array
    p = sortperm(points, by = item -> item.y)
	# Sort all other points according to angle with that point
    q = sortperm(points[p[2:end]], by = item -> atan(item.y - points[p[1]].y,item.x - points[p[1]].x))

	p[2:end] = p[q.+1]
	#@info " CONVEX POLY INITAL PIV $p "
	# M will be the point on the hull
    M = 2
	for i = 1:N
		do_while = false
        while (ccw(points[p[M-1]], points[p[M]], points[p[i]]) < 0.0)
			do_while = true
            if (M > 2)
				#@info " REDUCE BY 1 $M "
                M -= 1
            # All points are collinear
            elseif (i == N)
				#@info " ALL BPS COLLINEAR "
                break
            else
				#@info " NEXT "
                i += 1
            end
        end
        # ccw point found, updating hull and swapping points
        M += 1
		if do_while
			#@info " M NOW $M i $i swap --- PIV $p"
			p[i],p[M] = p[M], p[i]
		end
		M == N && break
    end
	return p
end

convex_polygon(v::Vector{Vector{dbl}}) = graham_scan!([POINT(j) for j in v])
# *************************************************
function rayintersectseg(p::POINT, edge::SEG)
	a, b = edge
	if a.y > b.y
		a, b = b, a
	end
	#@info " POINT $p"
	px = p.x ∈ (a.x,b.x) ? p.x + eps(p.x) : p.x# == a.x ? p.x + eps(p.x) : p.x == b.x ? p.x - eps(p.x) : p.x
	py = p.y ∈ (a.y,b.y) ? p.y + eps(p.y) : p.y # == a.y ? p.y + eps(p.y) : p.y == b.y ? p.y - eps(p.y) : p.y
	p = POINT([px,py])
	if (p.y > b.y || p.y < a.y) || (p.x > max(a.x, b.x))
		return false
	end

	if p.x < min(a.x, b.x)
		true
	else
		mred = (b.y - a.y) / (b.x - a.x)
		mblu = (p.y - a.y) / (p.x - a.x)
		return mblu ≥ mred
	end
end

function plane_intersection(verts::Vector{Vector{Float64}}, seg::SEG, tol::dbl)

	pos = in_face(verts, seg[2],tol)
	pos[1] >= 0 ? (return (pos[1]+1,pos[2], seg[2]) ) : pos == [-1,1] && (return (0,0,seg[2]) )

	pos = in_face(verts, seg[1], tol)
	pos[1] >= 0 ? (return (pos[1]+1,pos[2], seg[1]) ) : pos == [-1,1] && (return (0,0,seg[1]) )

	ray    = seg[2] - seg[1]
	dir    = cross(verts[2] .- verts[1], verts[3] .- verts[1])
	normal = dir ./ norm(dir)
	dotNR  = dot(normal,ray)
	abs(dotNR) < fTOL && (return () )

	pXYZ  = (verts[1] .+ verts[2].+verts[3]) ./ 3
    w     = seg[1] - pXYZ
	point = w .- ray .* dot(normal, w) / dotNR .+ pXYZ
	!in_segment(point, (seg[1],seg[2]),tol) && (return ()) # intersection point beyond segment length
	pos = in_face(verts,point, tol)
	pos == [-1,0] && (return ()) #throw(@error " intersection point is not in plane !!! ")
	pos[1] > 0 ? (return (2,pos[2],point)) : pos[1] == 0 ? (return (1,pos[2],point) ) :
				 (return (0,0,point))
end
