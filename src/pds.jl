
#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	pds.jl
	Pipe dream module of SchubertPolynomials.jl
	Hugh Dennin, 30 July 2024.
	
	This file contains the definition and methods for the struct PD, a type for (double) pipe dreams.
	A PD object `p` has two fields, `p.x` and `p.y`, each of which has type Vector{Tuple{Int8,Int8}}.
	These vectors record the locations of the crosses of x/row weight and y/column weight in p, respectively.
	
	IMPORTANT NOTE: `p.x` and `p.y` always are sorted in reading word order and do not contain duplicate entries.
	For this reason, only modify p.x and p.y with the included mutators add_x!, delete_x!, add_y!, and delete_y!, or else!
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



	#- - - - - - - - - - - - - - - - - - - -#
	#	Type, Constructors, & Iterators		#
	#- - - - - - - - - - - - - - - - - - - -#

export PD

struct PD
    x::Vector{Tuple{Int8,Int8}}
	y::Vector{Tuple{Int8,Int8}}
	
	function PD(xvect::Vector{<:Tuple{Vararg{Integer}}},
				yvect::Vector{<:Tuple{Vararg{Integer}}})
		xsorted = sort(unique(xvect),by=p->(p[1],-p[2]))
		ysorted = sort(unique(yvect),by=p->(p[1],-p[2]))
		new(xsorted,ysorted)
	end
end

# constructor for the empty pipe dream
PD() = PD(Tuple{Int8,Int8}[],Tuple{Int8,Int8}[])

# constructor with only the x-crosses given as a vector
PD(xvect::Vector{<:Tuple{Vararg{Integer}}}) = PD(xvect,Tuple{Int8,Int8}[])



		# * * Base Methods * * #

# redefine equality for pipe dreams
Base.:(==)(p1::PD, p2::PD) = (p1.x == p2.x) && (p1.y == p2.y)

# redefine Base.show for pipe dream display
function Base.show(io::IO, pd::PD)
    #println(io)
	
	b = bounds(pd)
	i0 = min(0,b[1])
	j0 = min(0,b[3])
	d = max(b[5],0)
	
	NEW_LINE = true
	
	for i in i0:d-j0
		for j in j0:d-i
			# xy-crosses (printed with a purple "Z")
			if (i,j) in pd.x && (i,j) in pd.y
				printstyled(io, "Z"; color=:magenta)
			# x-crosses (printed with a blue "X")
			elseif (i,j) in pd.x
				printstyled(io, "X"; color=:blue)
			# y-crosses (printed with a red "Y")
			elseif (i,j) in pd.y
				printstyled(io, "Y"; color=:red)
			# bump
			else
				# bump at the origin
				if i==0 && j==0
					if i0==0 && j0==0
						print(io, "┌")
					elseif i0==0
						print(io, "┬")
					elseif j0==0
						print(io, "├")
					else
						print(io, "┼")
					end
				# bumps on the x-axis
				elseif i==0
					# if the top row is the x-axis (i.e. i0==0)
					# and there is no cross in position (0,d),
					# then the rightmost "─" in position (0,d) should not be printed
					if i0==0 && j==d
						continue
					else
						print(io, "─")
					end
				# bumps on the y-axis
				elseif j==0
					# if the leftmost column is the y-axis (i.e. j0==0)
					# and there is no cross in position (d,0),
					# then the bottom "│" in position (d,0) should not be printed
					if j0==0 && i==d
						NEW_LINE = false
						continue
					else
						print(io, "│")
					end
				# bumps at positions off the axes
				else
					print(io, "·")
				end
			end
        end
		if NEW_LINE println(io) end
    end
end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--

	julia> pd1 = PD([(1,1), (1,2), (2,2), (4,1)])
	┌────
	│XX··
	│·X·
	│··
	│X
	
	julia> pd2 = PD([(1,1), (1,2), (2,2), (4,1)], [(1,3), (2,2), (2,3)])
	┌────
	│XXY·
	│·ZY
	│··
	│X
	
	julia> pd3 = PD([(-2,4), (-1,0), (-1,-2), (0,2), (2,-2), (2,-1), (2,0), (2,2), (3,1)])
	··│···X··
	X·X·····
	──┼─X──
	··│···
	XXX·X
	··│X
	··│
	··
	·
	
	julia> pd4 = PD([(2,2), (1,2), (1,1), (1,1), (4,1)])
	┌────
	│XX··
	│·X·
	│··
	│X


	julia> pd1 == pd4
	true
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



		# * * Iterators * * #

export PDs, K_PDs

# returns the information (k,is_standard) about the ladder in pd with southwest corner (i,j)
# if (i,j+1) is a bump, k is the smallest value so that (i-k,j) is a bump, otherwise k==0
# is_standard==true unless either k=0 OR (i-l,j+1) is a bump for some 0<l<k, in which case is_standard==false
# (i,j) should be a cross in pd
function ladder_type(pd::PD, (i,j)::Tuple{Vararg{Integer}})
	c = crosses(pd)
	if (i,j+1) in c return (0,false) end
	
	k=1
	is_standard=true
	while (i-k,j) in c
		if !((i-k,j+1) in c) is_standard=false end
		k+=1
	end
	return (k,is_standard)
end

# TO-DO: chute moves?


# iterator for all reduced pipe dreams for the permutation w
struct PDs
	w::Permutation
end
PDs(w::Vector{<:Integer}) = PDs(Permutation(w))

function Base.iterate(iter::PDs, state=([bottom_pipe_dream(iter.w)],[]))
	stack, found = state
	
	# end of iteration
	if isempty(stack) return nothing end
	
	pd = pop!(stack)
	push!(found,pd)
	
	for (i,j) in pd.x
		(k,is_standard) = ladder_type(pd,(i,j))
		if k > 0 && is_standard && i-k > 0
			new_pd = deepcopy(pd)
			delete_x!(new_pd,(i,j))
			add_x!(new_pd,(i-k,j+1))
			if !(new_pd in [stack;found]) push!(stack,new_pd) end
		end
	end
	
	return (pd,(stack,found))
end

Base.length(iter::PDs) = Base.SizeUnknown()
Base.eltype(::Type{PDs}) = PD

# iterator for all (K-theoretic) pipe dreams for the permutation w
struct K_PDs
	w::Permutation
end
K_PDs(w::Vector{<:Integer}) = K_PDs(Permutation(w))

function Base.iterate(iter::K_PDs, state=([bottom_pipe_dream(iter.w)],[]))
	stack, found = state
	
	# end of iteration
	if isempty(stack) return nothing end
	
	pd = pop!(stack)
	push!(found,pd)
	
	for (i,j) in pd.x
		(k,is_standard) = ladder_type(pd,(i,j))
		if k > 0 && is_standard && i-k > 0
			new_pd1 = deepcopy(pd)
			delete_x!(new_pd1,(i,j))
			add_x!(new_pd1,(i-k,j+1))
			if !(new_pd1 in [stack;found]) push!(stack,new_pd1) end
			
			new_pd2 = deepcopy(pd)
			add_x!(new_pd2,(i-k,j+1))
			if !(new_pd2 in [stack;found]) push!(stack,new_pd2) end
		end
	end
	
	return (pd,(stack,found))
end

Base.length(iter::K_PDs) = Base.SizeUnknown()
Base.eltype(::Type{K_PDs}) = PD

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--

	w = [1,4,3,2];

	julia> for pd in PDs(w) print(pd) end
	┌───
	│···
	│XX
	│X
	┌───
	│·X·
	│XX
	│·
	┌───
	│··X
	│X·
	│X
	┌───
	│·XX
	│··
	│X
	┌───
	│·XX
	│·X
	│·


	julia> for pd in K_PDs(w) show(pd) end
	┌───
	│···
	│XX
	│X
	┌───
	│·X·
	│XX
	│X
	┌───
	│·X·
	│XX
	│·
	┌───
	│··X
	│XX
	│X
	┌───
	│·XX
	│XX
	│X
	┌───
	│·XX
	│XX
	│·
	┌───
	│··X
	│X·
	│X
	┌───
	│·XX
	│X·
	│X
	┌───
	│·XX
	│··
	│X
	┌───
	│·XX
	│·X
	│X
	┌───
	│·XX
	│·X
	│·

+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



		# * * Mutators * * #

export add_x!, add_y!, delete_x!, delete_y!

# TO-DO: make p.x & p.y private fields

# adds an x-cross to pd at position (i,j)
function add_x!(pd::PD, (i,j)::Tuple{Vararg{Integer}})
	splice!(pd.x, searchsorted(pd.x, (i,j), by=p->(p[1],-p[2])), [(i,j)])
	return pd
end

# adds a y-cross to pd at position (i,j)
function add_y!(pd::PD, (i,j)::Tuple{Vararg{Integer}})
	splice!(pd.y, searchsorted(pd.y, (i,j), by=p->(p[1],-p[2])), [(i,j)])
	return pd
end

# removes the x-cross in pd at position (i,j) (if any)
function delete_x!(pd::PD, (i,j)::Tuple{Vararg{Integer}})
	deleteat!(pd.x, searchsorted(pd.x, (i,j), by=p->(p[1],-p[2])))
	return pd
end

# removes the y-cross in pd at position (i,j) (if any)
function delete_y!(pd::PD, (i,j)::Tuple{Vararg{Integer}})
	deleteat!(pd.y, searchsorted(pd.y, (i,j), by=p->(p[1],-p[2])))
	return pd
end

# vector versions of the above
function add_x!(pd::PD, crosses::Vector{<:Tuple{Vararg{Integer}}})
	for (i,j) in crosses add_x!(pd,(i,j)) end
	return pd
end
function add_y!(pd::PD, crosses::Vector{<:Tuple{Vararg{Integer}}})
	for (i,j) in crosses add_y!(pd,(i,j)) end
	return pd
end
function delete_x!(pd::PD, crosses::Vector{<:Tuple{Vararg{Integer}}})
	for (i,j) in crosses delete_x!(pd,(i,j)) end
	return pd
end
function delete_y!(pd::PD, crosses::Vector{<:Tuple{Vararg{Integer}}})
	for (i,j) in crosses delete_y!(pd,(i,j)) end
	return pd
end


#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> pd = PD()
	┌

	julia> add_x!(pd,(2,1))
	┌──
	│··
	│X

	julia> add_x!(pd,[(1,2),(1,3),(2,2)])
	┌───
	│·XX
	│XX
	│·

	julia> add_y!(pd,[(1,1),(1,2),(3,2)])
	┌────
	│YZX·
	│XX·
	│·Y
	│·
	
	julia> delete_x!(pd,[(1,2),(1,3)])
	┌────
	│YY··
	│XX·
	│·Y
	│·

	julia> delete_y!(pd,(1,2))
	┌────
	│Y···
	│XX·
	│·Y
	│·
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#




	#- - - - - - - - - - - - - -#
	#	Pipe Dream Methods		#
	#- - - - - - - - - - - - - -#

		# * * General Methods * * #

export crosses, rotate180, translate, swap

# returns the positions of all (both x- and y-) crosses in pd without repetition
function crosses(pd::PD)
	c = deepcopy(pd.x)
	for (i,j) in pd.y
		splice!(c, searchsorted(c, (i,j), by=p->(p[1],-p[2])), [(i,j)])
	end
	return c
end

# returns the min/max row (= i), min/max column (= j), and max diagonal (= i+j) among all crosses (i,j) appearing in pd
function bounds(pd::PD)
	c = crosses(pd)
	
	if isempty(c) return [1,1,1,1,1] end
	
	min_row = minimum([i for (i,j) in c])
	max_row = maximum([i for (i,j) in c])
	
	min_col = minimum([j for (i,j) in c])
	max_col = maximum([j for (i,j) in c])
	
	max_diag = maximum([i+j for (i,j) in c])
	
	return [min_row,max_row,min_col,max_col,max_diag]
end

# returns the transpose of the pipe dream pd
function Base.transpose(pd::PD)
	xvec = [(j,i) for (i,j) in pd.x]
	yvec = [(j,i) for (i,j) in pd.y]
	return PD(xvec,yvec)
end

# returns the pipe dream pd rotated through 180 degrees
function rotate180(pd::PD)
	xvec = [(-i,-j) for (i,j) in pd.x]
	yvec = [(-i,-j) for (i,j) in pd.y]
	return PD(xvec,yvec)
end

# returns the pipe dream pd translated by the offset (i0,j0)
function translate(pd::PD,(i0,j0))
	xvec = [(i+i0,j+j0) for (i,j) in pd.x]
	yvec = [(i+i0,j+j0) for (i,j) in pd.y]
	return PD(xvec,yvec)
end

# returns the pipe dream pd with its x- and y-crosses swapped
function swap(pd::PD)
	return PD(pd.y,pd.x)
end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--

	julia> pd = PD([(1,1), (1,2), (2,2), (4,1)], [(1,3), (2,2), (2,3)])
	┌────
	│XXY·
	│·ZY
	│··
	│X
	
	julia> println(crosses(pd))
	Tuple{Int8, Int8}[(1, 3), (1, 2), (1, 1), (2, 3), (2, 2), (4, 1)]

	julia> println(bounds(pd))
	Int8[1, 4, 1, 3, 5]
	
	julia> transpose(pd)
	┌────
	│X··X
	│XZ·
	│YY
	│·

	julia> rotate180(pd)
	··X│····
	···│···
	YZ·│··
	YXX│·
	───┼
	···
	··
	·
	
	julia> translate(pd,(-5,1))
	│·XXY·
	│··ZY
	│···
	│·X
	├─
	
	julia> swap(pd)
	┌────
	│YYX·
	│·ZX
	│··
	│Y
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



		# * * Permutations * * #

export bottom_pipe_dream, reading_word, permutation

# returns the bottom dream for the permutation w
function bottom_pipe_dream(w::Permutation)
	c = code(w)
	n = length(c)
	x = Tuple{Int,Int}[]
	
    for i in 1:n
		for j in 1:c[i]
			push!(x,(i,j))
		end
	end
	
    return(PD(x))
end
function bottom_pipe_dream(w::Vector{<:Integer}) return bottom_pipe_dream(Permutation(w)) end

# TO-DO: write function top_pipe_dream(w::Permutation)

# returns the reading word of pd
function reading_word(pd::PD) return [(i,i+j-1) for (i,j) in crosses(pd)] end

# returns the permutation of pd
function permutation(pd::PD) return demazure_product([i+j-1 for (i,j) in crosses(pd)]) end

function is_reduced(pd::PD) return len(permutation(pd))==length(pd.x)+length(pd.y) end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> bottom_pipe_dream([1,3,5,2,4])
	┌────
	│····
	│X··
	│XX
	│·
	
	julia> pd1 = PD([(1,2),(2,3),(3,1)])
	┌────
	│·X··
	│··X
	│X·
	│·

	julia> println(reading_word(pd1))
	Tuple{Int8, Int64}[(1, 2), (2, 4), (3, 3)]

	julia> permutation(pd1)
	Permutation(1, 3, 5, 2, 4)

	julia> is_reduced(pd1)
	true

	julia> pd2 = PD([(1,2),(1,4),(2,1),(3,1),(3,2)])
	┌────
	│·X·X
	│X··
	│XX
	│·

	julia> permutation(pd2)
	Permutation(1, 3, 5, 2, 4)

	julia> is_reduced(pd2)
	false
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



		# * * Counting Pipe Dreams * * #

export	pipe_dream_multiplicities,		display_pipe_dream_multiplicities,
		K_pipe_dream_multiplicities,	display_K_pipe_dream_multiplicities

# given a vector "pds" of pipe dreams, returns a vector "mults" of tuples (w,m) as well as another vector of pipe dreams "remainder"
# -	a tuple (w,m) appearing in "mults" consists of a permutation w and a multiplicity m that is the greatest value so that m distinct copies of PDs(w) appear in pds
# - remainder contains the remaining pipe dreams in pds after m copies of PDs(w) are subtracted from pds for every (w,m) in mults
function pipe_dream_multiplicities(pds::Vector{PD})
	stack = deepcopy(pds)
	mults = Tuple{Permutation,Int8}[]
	remainder = PD[]
	
	while !isempty(stack)
		w = permutation(stack[1])
		m = length(stack)
		
		for p in PDs(w)
			indices = findall(q -> q == p, stack)
			m = min(m,length(indices))
		end
		
		push!(mults,(w,m))
		
		for p in PDs(w)
			indices = findall(q -> q == p, stack)
			deleteat!(stack,indices[1:m])
		end
		
		indices = findall(q -> permutation(q) == w, stack)
		append!(remainder,stack[indices])
		deleteat!(stack,indices)
	end
	
	return (mults,remainder)
end

function display_pipe_dream_multiplicities(pds::Vector{PD})
	(mults,remainder) = pipe_dream_multiplicities(pds)
	println("Multiplicies:")
	for (w,m) in mults
		print(m)
		print(" x ")
		println(w)
	end
	println()
	println("Remaining Pipe Dreams:")
	for p in remainder
		show(p)
	end
end

function K_pipe_dream_multiplicities(pds::Vector{PD})
	stack = deepcopy(pds)
	mults = Tuple{Permutation,Int8}[]
	remainder = PD[]
	
	while !isempty(stack)
		w = permsutation(stack[1])
		m = length(stack)
		
		for p in K_PDs(w)
			indices = findall(q -> q == p, stack)
			m = min(m,length(indices))
		end
		
		push!(mults,(w,m))
		
		for p in K_PDs(w)
			indices = findall(q -> q == p, stack)
			deleteat!(stack,indices[1:m])
		end
		
		indices = findall(q -> permutation(q) == w, stack)
		append!(remainder,stack[indices])
		deleteat!(stack,indices)
	end
	
	return (mults,remainder)
end

function display_K_pipe_dream_multiplicities(pds::Vector{PD})
	(mults,remainder) = K_pipe_dream_multiplicities(pds)
	println("Multiplicies:")
	for (w,m) in mults
		print(m)
		print(" x ")
		println(w)
	end
	println()
	println("Remaining Pipe Dreams:")
	for p in remainder
		show(p)
	end
end




	#- - - - - - - - - - - - - - - -#
	#	Reduced Pipe Dream Methods	#
	#- - - - - - - - - - - - - - - -#

# These methods are currently only implemented for reduced pipe dreams with crosses in the positive quadrant
# TO-DO: Add checks to test if inputs are valid

		# * * General Methods  * * #

export south_pipe, west_pipe, north_pipe, east_pipe

# returns the row index of the pipe in pd passing through the south edge of position (i,j)
function south_pipe(pd::PD, (i,j)::Tuple{Vararg{Integer}})
	c = crosses(pd)
	(i0,j0) = (i+1,j)
	dir = (1,0)
	
	while j0>0
		if !((i0,j0) in c)
			if dir==(1,0)		dir=(0,-1)
			elseif dir==(0,-1)	dir=(1,0)
			end
		end
		(i0,j0) = (i0+dir[1],j0+dir[2])
	end
	
	return i0
end

# returns the row index of the pipe in pd passing through the west edge of position (i,j)
function west_pipe(pd::PD, (i,j)::Tuple{Vararg{Integer}})
	c = crosses(pd)
	(i0,j0) = (i,j-1)
	dir = (0,-1)
	
	while j0>0
		if !((i0,j0) in c)
			if dir==(1,0)		dir = (0,-1)
			elseif dir==(0,-1)	dir = (1,0)
			end
		end
		(i0,j0) = (i0+dir[1],j0+dir[2])
	end
	
	return i0
end

# returns the column index of the pipe in pd passing through the north edge of position (i,j)
function north_pipe(pd::PD, (i,j)::Tuple{Vararg{Integer}})
	c = crosses(pd)
	(i0,j0) = (i-1,j)
	dir = (-1,0)
	
	while i0>0
		if !((i0,j0) in c)
			if dir==(-1,0)		dir=(0,1)
			elseif dir==(0,1)	dir=(-1,0)
			end
		end
		(i0,j0) = (i0+dir[1],j0+dir[2])
	end
	
	return j0
end

# returns the column index of the pipe in pd passing through the east edge of position (i,j)
function east_pipe(pd::PD, (i,j)::Tuple{Vararg{Integer}})
	c = crosses(pd)
	(i0,j0) = (i,j+1)
	dir = (0,1)
	
	while i0>0
		if !((i0,j0) in c)
			if dir==(-1,0)		dir = (0,1)
			elseif dir==(0,1)	dir = (-1,0)
			end
		end
		(i0,j0) = (i0+dir[1],j0+dir[2])
	end
	
	return j0
end




	#- - - - - - - - - - - - - -#
	#	Insertion Algorithms	#
	#- - - - - - - - - - - - - -#

		# * * Edelman Greene * * #
		
export edelman_greene

# applies the Edelman Greene map to the reading word of the pipe dream
function edelman_greene(pd::PD) return edelman_greene(reading_word(pd)) end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> pd = PD([(1,4),(1,2),(2,4),(2,2),(2,1),(3,2),(3,1)])
	┌─────
	│·X·X·
	│XX·X
	│XX·
	│··
	│·


	julia> (P,Q) = edelman_greene(pd)
	(
	2 3 4
	3 4 5
	4
	,
	1 1 2
	2 2 3
	3
	)

+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



