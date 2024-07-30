
#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	permutations.jl
	Permutation module of SchubertPolynomials.jl
	Hugh Dennin, 30 July 2024.
	
	This file contains the definition and methods for the struct Permutation, a type for permutations of \mathbb{Z} with finite support.
	A Permutation object `w` consists of a single field `w.rels` which has type Vector{Tuple{Int8,Int8}}.
	w maps i to j whenever (i,j) appears in w.rels, whereas every other input is fixed by w.
	
	IMPORTANT NOTE: `w.rels` is always sorted by the first coordinate and contains no duplicates.
	Do not attempt to make changes to w.rels manually unless you know what you're doing.
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



	#- - - - - - - - - - - - - - - - - - - -#
	#	Type, Constructors, & Iterators		#
	#- - - - - - - - - - - - - - - - - - - -#

export Permutation, Permutations

struct Permutation
    rels::Vector{Tuple{Int8,Int8}}
	
	function Permutation(rels::Vector{<:Tuple{Vararg{Integer}}})
		unique_rels = unique(rels)
		
		# check for validity: each number should appear at most once in each position
		if length(unique(first.(unique_rels))) != length(unique_rels)
			throw(ArgumentError("Invalid permutation (input appears multiple times)"))
		end
		if length(unique(last.(unique_rels))) != length(unique_rels)
			throw(ArgumentError("Invalid permutation (value appears multiple times)"))
		end
		
		# remove fixed points and sort by input
		cleaned_rels = sort(filter(r -> r[1] != r[2], unique_rels), by = first)
		new(cleaned_rels)
	end
end

# constructor for identity permutation
Permutation() = Permutation(Tuple{Int8,Int8}[])

# constructor from a vector of integers
function Permutation(w::Vector{<:Integer})
	if isempty(w) Permutation() end
	i0 = minimum(w)
    rels = [(i0+i-1,w[i]) for i in 1:length(w)]
    Permutation(rels)
end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--

	julia> Permutation([(1,1),(2,4),(3,2),(4,3)])
	Permutation(1, 4, 2, 3)
	
	julia> Permutation([(1,1),(1,4),(3,2),(4,3)])
	ERROR: ArgumentError: Invalid permutation (input appears multiple times)
	
	Permutation([1,3,5,2,4])
	Permutation(1, 3, 5, 2, 4)
	
	julia> Permutation([1,3,3,2,4])
	ERROR: ArgumentError: Invalid permutation (value appears multiple times)
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#

# iterator for permutations of a fixed size n that indexes by the code of the permutation.
struct Permutations
	n::Int8
end

function Base.iterate(iter::Permutations, state=zeros(Int8,iter.n))
	code = state
    w = code_to_perm(code)
	
	# End of iteration
    if code == vcat(zeros(Int8,iter.n-1),Int8[1])
        return nothing
    end
	
    # Find the next code.
	k = 1
	while code[k] == iter.n-k && k != iter.n
		code[k] = 0
		k+=1
	end
	code[k] += 1
	
	return (w,code)
end

Base.length(iter::Permutations) = factorial(iter.n)
Base.eltype(::Type{Permutations}) = Permutation

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> collect(Permutations(3))
	6-element Vector{Permutation}:
	 Permutation(1)
	 Permutation(2, 1)
	 Permutation(3, 1, 2)
	 Permutation(1, 3, 2)
	 Permutation(2, 3, 1)
	 Permutation(3, 2, 1)
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



	# * * Base Methods * * #

# redefine Base.getindex
function Base.getindex(w::Permutation, index::Integer)
    for (i,j) in w.rels
        if i == index
            return j
        end
    end
    return index  # If index is not found above, it's a fixed point
end

function Base.getindex(w::Permutation, indices::Vector{<:Integer})
    return [w[index] for index in indices]
end

# redefine Base.show
function Base.show(io::IO, w::Permutation)
	if isempty(w.rels)
		print(io, "Permutation(1)")
	else
		(p,q) = bounds(w)
		print(io, "Permutation(", join([w[i] for i in min(1,p):q], ", "), ")")
	end
end

# redefine equality
Base.:(==)(w1::Permutation, w2::Permutation) = w1.rels == w2.rels

# composition of permutations
# TO-DO: make composition faster
function *(u::Permutation, v::Permutation)
	(pu,qu) = bounds(u)
	(pv,qv) = bounds(v)
	p, q = min(pu,pv), max(qu,qv)
    return Permutation([(i,u[v[i]]) for i in p:q if u[v[i]] != i])
end

# redefine Base.convert
Base.convert(::Type{Permutation}, w::Vector{<:Integer}) = Permutation(w)

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> w = Permutation([(1,1),(2,4),(3,2),(4,3)])
	Permutation(1, 4, 2, 3)

	julia> u = Permutation([1,4,2,3])
	Permutation(1, 4, 2, 3)
	
	julia> w == u
	true
	
	julia> w*w
	Permutation(1, 3, 4, 2)
	
	julia> w*w*w
	Permutation(1)
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#
	



	#- - - - - - - - - - - - - -#
	#	Permutation Methods		#
	#- - - - - - - - - - - - - -#
	
		# * * Special Permutations * * #

export cycle, t_ij, s_i, longest_perm

# Returns the cycle (c[1] c[2] ... c[k])
function cycle(c::Vector{<:Integer})
	if isempty(c) Permutation() end
	k = length(c)
	rels = [(c[k],c[1])]
	for i in 1:k-1 push!(rels,(c[i],c[i+1])) end
	return Permutation(rels)
end

# Returns the transposition t_{i,j} = (i j)
function t_ij(i::Integer,j::Integer)
	return cycle([i,j])
end

# Returns the simple transposition s_i = t_{i,i+1} = (i i+1)
function s_i(i::Integer)
	return t_ij(i,i+1)
end

function longest_perm(n::Integer)
	return Permutation([(i,n-i+1) for i in 1:n])
end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> cycle([3,4,2])
	Permutation(1, 3, 4, 2)

	julia> t_ij(2,4)
	Permutation(1, 4, 3, 2)
	
	julia> s_i(2)
	Permutation(1, 3, 2)
	
	julia> t_ij(2,4) == s_i(3)*s_i(2)*s_i(3)
	true
	
	julia> longest_perm(5)
	Permutation(5, 4, 3, 2, 1)
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



		# * * General Methods * * #

export support, bounds, inversions, code, len, code_to_perm, descents

# returns the support (= non-fixed set) of w.
function support(w::Permutation)
	return [i for (i,j) in w.rels]
end
function support(w::Vector{<:Integer}) return support(Permutation(w)) end

# returns the least and greatest values in the support of w.
function bounds(w::Permutation)
	return isempty(w.rels) ? (1,1) : (first(w.rels)[1],last(w.rels)[1])
end
function bounds(w::Vector{<:Integer}) return bounds(Permutation(w)) end

# returns the inversion set of w
function inversions(w::Permutation)
	(p,q) = bounds(w)
	Inv = Tuple{Int8,Int8}[]
	for i in p:q-1 for j in i+1:q
		if w[i] > w[j] push!(Inv,(i,j)) end
	end end
	return Inv
end
function inversions(w::Vector{<:Integer}) return inversions(Permutation(w)) end

# returns the code of w
# (only works if w has positive support)
function code(w::Permutation)
	(p,q) = bounds(w)
	code = zeros(Int8,q-1)
	for (i,j) in inversions(w) code[i] += 1 end
	return code
end
function code(w::Vector{<:Integer}) return code(Permutation(w)) end

# returns the Coxeter length of w
function len(w::Permutation)
	return length(inversions(w))
end
# implemented in "permtools.jl"
#function len(w::Vector{<:Integer}) return len(Permutation(w)) end

# returns the permutation with the prescribed code.
function code_to_perm(code::Vector{<:Integer})
	w = Integer[]
	
	k = length(code)
	n = maximum([i+code[i] for i in 1:k])
	avail = collect(1:n)
	
	for i in 1:length(code)
		push!(w,avail[1+code[i]])
		deleteat!(avail,1+code[i])
	end
	
	return Permutation(vcat(w,avail))
end

# Returns the descent set of w
function descents(w::Permutation)
	(p,q) = bounds(w)
	Des = Int8[]
	for i in p:q-1
		if w[i] > w[i+1] push!(Des,i) end
	end
	return Des
end
# implemented in "permtools.jl"
#function descents(w::Vector{<:Integer}) return descents(Permutation(w)) end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--

	julia> w = Permutation([1,3,5,2,7,6,4]);

	julia> println(support(w))
	Int8[2,3,4,5,7]

	julia> bounds(w)
	(2, 7)

	julia> println(inversions(w))
	Tuple{Int8, Int8}[(2, 4), (3, 4), (3, 7), (5, 6), (5, 7), (6, 7)]

	julia> println(code(w))
	Int8[0,1,2,0,2,1]

	julia> len(w)
	6

	julia> code_to_perm([0,1,2,0,2,1])
	Permutation(1, 3, 5, 2, 7, 6, 4)

	julia> println(descents(w))
	Int8[3,5,6]

+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



		# * * Reduced Words * * #
		
export coxeter_product, demazure_product, is_reduced

function coxeter_product(word::Vector{<:Integer})
	w = Permutation()
	for i in word w = w*s_i(i) end
	return w
end

function demazure_product(word::Vector{<:Integer})
	w = Permutation()
	for i in word if w[i] < w[i+1] w = w*s_i(i) end end
	return w
end

function is_reduced(word::Vector{<:Integer})
	w = coxeter_product(word)
	return len(w) == length(word)
end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> coxeter_product([3,2,3,1,3])
	Permutation(4, 1, 2, 3)
	
	julia> demazure_product([3,2,3,1,3])
	Permutation(4, 1, 3, 2)
	
	julia> is_reduced([3,2,3,1,3])
	false
	
	julia> is_reduced([3,2,3,1])
	true
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



		# * * Permutation Patterns * * #

export max_pattern, max_inversion, max_132, max_2143

# finds the lex-final occurance of u as a pattern in w and returns its indices (or [] if w avoids u)
# w must be a permutation with positive support, and is only considered to be a permutation of the positive intergers
# (i.e. w = 21 does not contain a 132-pattern even though it is technically w = ...021...)
# u must be a permutation of the numbers 1,2,...,n with u[n] != n

function max_pattern(w::Permutation,u::Vector{<:Integer})
	# Maybe implement a check to see if u is a good pattern?
	n = bounds(w)[2]
	k = length(u)
	kvec = collect(1:k)
	
	# Deal with cases where n <= k.
	if n<k return Int8[]
	elseif n==k
		w[kvec]==u && return kvec
		return Int8[]
	end
	
	# It would be nice to use a subset iterator here...
	indices = Int8[n-i+1 for i in k:-1:1]
	while indices != kvec
		if invperm(sortperm(w[indices]))==u return indices end
		
		i = 1
		while i != k && indices[k-i+1] == indices[k-i]+1
			indices[k-i+1] = n-i+1
			i += 1
		end
		indices[k-i+1] -= 1
	end
	
	return Int8[]
end
function max_pattern(w::Vector{<:Integer},u::Vector{<:Integer}) return max_pattern(Permutation(w),u) end

# Specific instances of max_pattern(w,u) for the patterns u = 21, 132, and 2143.
function max_inversion(w::Permutation) return max_pattern(w,[2,1]) end
function max_inversion(w::Vector{<:Integer}) return max_pattern(w,[2,1]) end
function max_132(w::Permutation) return max_pattern(w,[1,3,2]) end
function max_132(w::Vector{<:Integer}) return max_pattern(w,[1,3,2]) end
function max_2143(w::Permutation) return max_pattern(w,[2,1,4,3]) end
function max_2143(w::Vector{<:Integer}) return max_pattern(w,[2,1,4,3]) end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> w = Permutation([2,5,6,10,1,8,3,4,9,7]);
	
	julia> println(max_inversion(w))
	Int8[8, 10]
	
	julia> println(max_132(w))
	Int8[7, 8, 10]

	julia> println(max_2143(w))
	Int8[2, 7, 8, 10]
	
	julia> println(max_pattern(w,[1,4,3,2]))
	Int8[2, 3, 8, 10]
	
	# w avoids the pattern [1,5,4,3,2]
	julia> println(max_pattern(w,[1,5,4,3,2]))
	Int8[]
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#


		# * * Transition * * #

export dominant_transition, max_transition

# Do step in Fan-Guo-Sun transition tree.
function dominant_transition(w::Permutation)
	pos = max_132(w)
	if isempty(pos) return [] end
	
	l = len(w)
	r,s = pos[2],pos[3]
	v = w*t_ij(r,s)
	
	transition_perms = Permutation[]
	for i in 1:r-1
		u = v*t_ij(i,r)
		if len(u)==l push!(transition_perms,u) end
	end
	
	return [(r,s), v, transition_perms]
end
# implemented in "permtools.jl"
#function dominant_transition(w::Vector{<:Integer}) dominant_transition(Permutation(w)) end

# Do step in Lascoux-Schutzenberger transition tree.
function max_transition(w::Permutation)
	pos = max_2143(w)
	if isempty(pos) return [] end
	
	l = len(w)
	r,s = pos[3],pos[4]
	v = w*t_ij(r,s)
	
	transition_perms = Permutation[]
	for i in 1:r-1
		u = v*t_ij(i,r)
		if len(u)==l push!(transition_perms,u) end
	end
	
	return [(r,s), v, transition_perms]
end
# implemented in "permtools.jl"
#function max_transition(w::Vector{<:Integer}) max_transition(Permutation(w)) end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> w = Permutation([2,1,5,4,3,8,6,7]);

	julia> dominant_transition(w)
	3-element Vector{Any}:
	 (6, 8)
	 Permutation(2, 1, 5, 4, 3, 7, 6)
	 Permutation[Permutation(2, 1, 7, 4, 3, 5, 6), Permutation(2, 1, 5, 7, 3, 4, 6), Permutation(2, 1, 5, 4, 7, 3, 6)]

	julia> max_transition(w)
	3-element Vector{Any}:
	 (6, 8)
	 Permutation(2, 1, 5, 4, 3, 7, 6)
	 Permutation[Permutation(2, 1, 7, 4, 3, 5, 6), Permutation(2, 1, 5, 7, 3, 4, 6), Permutation(2, 1, 5, 4, 7, 3, 6)]
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#
