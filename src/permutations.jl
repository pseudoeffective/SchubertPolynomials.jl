
#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
    
	permutations.jl
    Permutation module of SchubertPolynomials.jl
    Hugh Dennin, 30 July 2024.
    
    This file contains the definition and methods for the struct Permutation, a type for permutations of \mathbb{Z} with finite support.
    A Permutation object `w` consists of a single field `w.data` of type Vector{Int8}.
	`w.data[2:length(w.data)] = [w[m], w[m+1], ..., w[M]]` is a permutation of some interval `[start,finish]` in \mathbb{Z}
    `w.data[1] = start` records the initial position of this interval.
    
    WARNING: Do not attempt to make changes to w.data manually unless you know what you're doing!

+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



    #- - - - - - - - - - - - - - - - - - - -#
    #    Type, Constructors, & Iterators    #
    #- - - - - - - - - - - - - - - - - - - -#

export Permutation, Permutations

struct Permutation
	data::Vector{Int8}
	
	# constructor taking a vector of integers [w[m], w[m+1], ..., w[M]] that gives a permutation of an interval [m,M]
	# inputs outside [m,M] are fixed
	function Permutation(w::Vector{<:Integer})
		return new(pushfirst!(w,minimum(w)))
	end
	
	function Permutation(w::Vector{<:Integer},m::Integer)
		return new(pushfirst!(w,m))
	end
end

# constructor for identity permutation
Permutation() = Permutation([1],1)

# constructor taking any number of integers w[m], w[m+1], ..., w[M] that give a permutation of an interval [m,M]
Permutation(I::Integer...) = Permutation([i for i in I])

# constructor taking a vector of integer tuples [(i1,w[i1]), (i2,w[i2]), ..., (ir,w[ir])], where inputs not given are fixed
function Permutation(rels::Vector{<:Tuple{Vararg{Integer}}})
	inputs, outputs = first.(rels), last.(rels)
	m, M = minimum(inputs), maximum(inputs)
	w = [i for i in m:M]
	for i in 1:length(inputs) w[inputs[i]-m+1] = outputs[i] end
	return Permutation(w,m)
end

# constructor taking any number of integer tuples (i1,w[i1]), (i2,w[i2]), ..., (ir,w[ir]), where inputs not given are fixed
Permutation(rels::Tuple{Integer,Integer}...) = Permutation([rel for rel in rels])

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> Permutation([1,3,5,2,4])
	Permutation(1, 3, 5, 2, 4)
	
	julia> Permutation(1,3,5,2,4)
	Permutation(1, 3, 5, 2, 4)
	
	julia> Permutation(-2,3,2,-1,0,1,4,5)
	Permutation(3, 2, -1, 0, 1)
	
	julia> Permutation([(2,4),(4,2)])
	Permutation(1, 4, 3, 2)
	
	julia> Permutation((1,1),(2,4),(3,2),(4,3))
	Permutation(1, 4, 2, 3)
	
	julia> Permutation()
	Permutation(1)
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



	# * * Base Methods * * #

# redefine Base.getindex
function Base.getindex(w::Permutation, i::Integer)
	real_i = i-start(w)+1
	if 1 <= real_i && real_i <= length(w.data)-1
		return w.data[1+real_i]
	end
	# otherwise i is a fixed point of w
	return Int8(i)
end
function Base.getindex(w::Permutation, I::Vector{<:Integer}) return Int8[w[i] for i in I] end
function Base.getindex(w::Permutation, I::UnitRange{<:Integer}) return Int8[w[i] for i in I] end

# redefine Base.show
function Base.show(io::IO, w::Permutation)
	m, M = min!(w), max!(w)
	entries = vcat(collect(1:m-1), [w[i] for i in m:M], collect(M+1:0))
	print(io,
		"Permutation(",
		join(entries, ", "),
		")"
	)
end

# redefine equality
function Base.:(==)(u::Permutation, v::Permutation)
	min!(u), max!(u), min!(v), max!(v)
	return u.data == v.data
end

# composition of permutations
function *(u::Permutation, v::Permutation)
	m = min(min!(u),min!(v))
	M = max(max!(u),max!(v))
    return Permutation([u[v[i]] for i in m:M], m)
end

# redefine Base.convert
Base.convert(::Type{Permutation}, w::Vector{<:Integer}) = Permutation(w)

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> w = Permutation(1,4,2,3)
	Permutation(1, 4, 2, 3)

	julia> u = Permutation((1,1),(2,4),(3,2),(4,3)])
	Permutation(1, 4, 2, 3)
	
	julia> v = Permutation(4,2,3,5,6,7,8,9)
	Permutation(1, 4, 2, 3)
	
	julia> w == u && w == v
	true
	
	julia> w*w
	Permutation(1, 3, 4, 2)
	
	julia> w*w*w
	Permutation(1)
	
	julia> Permutation(4,6,5,2,3,1) * Permutation(6,4,5,2,8,1,3,7)
	Permutation(1, 2, 3, 6, 8, 4, 5, 7)
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



		# * * Accessors / Mutators * * #
		
export min!, max!

# These functions are not well-defined on permutation equivalence. Use with caution!
function start(w::Permutation) return w.data[1] end
function entries(w::Permutation) return w.data[2:length(w.data)] end
		
# returns the min value in the support of w, or 1 if w==id
function min!(w::Permutation)
	n = length(w.data)-1
	i0 = start(w)
	for i in i0:i0+n-1
		if i != w[i]
			if i != i0
				w.data[1] = i
				deleteat!(w.data,2:1+i-i0) # ranges from 1 to 1+(i-i0-1) (with +1 data offset)
			end
			return i
		end
	end
	# otherwise w==id
	if w.data != [1,1]
		deleteat!(w.data,1:length(w.data))
		append!(w.data,[1,1])
	end
	return 1
end

# returns the max value in the support of w, or 1 if w==id
function max!(w::Permutation)
	n = length(w.data)-1
	i0 = n+start(w)-1
	for i in i0:-1:i0-n+1
		if i != w[i]
			if i != i0
				deleteat!(w.data,2+n+i-i0:1+n) # ranges from n-(i0-i-1) to n (with +1 data offset)
			end
			return i
		end
	end
	# otherwise w = identity
	if w.data != [1,1]
		deleteat!(w.data,1:length(w.data))
		append!(w.data,[1,1])
	end
	return 1
end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	# This is a technical example, normally w.data should not be accessed by the user.
	
	julia> w = Permutation(1,2,5,8,6,4,3,7,9,10,11);
	
	julia> println(w.data)
	# start index = 1, perm data = [1,2,5,8,6,4,3,7,9,10,11]
	Int8[1, 1, 2, 5, 8, 6, 4, 3, 7, 9, 10, 11]
	
	julia> min!(w)
	# smallest non-fixed point
	3
	
	julia> println(w.data)
	# start index = 3, perm data = [5,8,6,4,3,7,9,10,11] (initial [1,2] cut off)
	Int8[3, 5, 8, 6, 4, 3, 7, 9, 10, 11]
	
	julia> max!(w)
	# largest non-fixed point
	8

	julia> println(w.data)
	# start index = 3, perm data = [5,8,6,4,3,7] (final [9,10,11] cut off)
	Int8[3, 5, 8, 6, 4, 3, 7]
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



		# * * Iterators * * #
		
# iterator for permutations of a fixed size n that indexes by the code of the permutation.
struct Permutations
	n::Int8
end

function Base.iterate(iter::Permutations, state=zeros(Int8,max(iter.n,0)))
	# n<=0 cases
	if isempty(state)
		iter.n < 0 && return nothing
		return (Permutation(),[1])
	end
	
	# end of iteration
    last(state) > 0 && return nothing
	
    w = code_to_perm(state)
		
    # find the next state
	new_state = copy(state)
	k = 1
	while new_state[k] == iter.n-k && k != iter.n
		new_state[k] = 0
		k+=1
	end
	new_state[k] += 1
	
	return (w,new_state)
end

Base.length(iter::Permutations) = iter.n<0 ? 0 : factorial(iter.n)
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
	
	julia> length(Permutations(3))
	6
	
	julia> collect(Permutations(0))
	1-element Vector{Permutation}:
	 Permutation(1)
	 
	julia> collect(Permutations(-1))
	Vector{Permutation}[]
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



    #- - - - - - - - - - - - - -#
    #    Permutation Methods    #
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
function cycle(c::Integer...) return cycle([i for i in c]) end

# Returns the transposition t_{i,j} = (i j)
function t_ij(i::Integer,j::Integer)
	if i == j return Permutation() end
	m, M = min(i,j), max(i,j)
	return Permutation(vcat([M], collect(m+1:M-1), [m]), m)
end

# Returns the simple transposition s_i = t_{i,i+1} = (i i+1)
function s_i(i::Integer)
	return Permutation([i+1,i],i)
end

function longest_perm(n::Integer)
	return Permutation(collect(n:-1:1),1)
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



		# * * Group Operations * * #

export inverse, conjugate, cycle_decomposition, order, support


# returns the inverse permutation of w
function inverse(w::Permutation)
	m, M = min!(w), max!(w)
	return Permutation([(w[i],i) for i in m:M])
end
function inverse(w::Vector{<:Integer}) return inverse(Permutation(w)) end


# returns w*u*w^{-1}
function conjugate(w::Permutation, u::Permutation)
	return w*u*inverse(w)
end


# returns the cycle decomposition of w
function cycle_decomposition(w::Permutation)
	m, M = min!(w), max!(w)
	cycles = Vector{Int8}[]
	
	available = collect(m:M)
	while !isempty(available)
		i = popfirst!(available)
		if w[i] == i continue end
		c = [i]
		while w[i] != c[1]
			i = w[i]
			push!(c,i)
			deleteat!(available,searchsorted(available,i))
		end
		push!(cycles,c)
	end
	
	return cycles
end
function cycle_decomposition(w::Vector{<:Integer}) return cycle_decomposition(Permutation(w)) end


# returns the support (= non-fixed set) of w
function support(w::Permutation)
	return Int8[i for i in min!(w):max!(w) if i != w[i]]
end
function support(w::Vector{<:Integer}) return support(Permutation(w)) end



		# * * Permutation Statistics * * #
		
export inversions, code, len, is_even, code_to_perm, descents, maj


# returns the group order of w
function order(w::Permutation)
	cycles = cycle_decomposition(w)
	isempty(cycles) && return 0
	return lcm([length(c) for c in cycles])
end
function order(w::Vector{<:Integer}) return order(Permutation(w)) end


# returns the inversion set of w
function inversions(w::Permutation)
	Inv = Tuple{Int8,Int8}[]
	m, M = min!(w), max!(w)
	for i in m:M-1 for j in i+1:M
		if w[i] > w[j] push!(Inv,(i,j)) end
	end end
	return Inv
end
function inversions(w::Vector{<:Integer}) return inversions(Permutation(w)) end


# returns the code of w
# NOTE: only use when w has positive support
function code(w::Permutation)
	code = zeros(Int8,max!(w)-1)
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


function is_even(w::Permutation) return mod(len(w),2) == 0 end


# returns the permutation with the prescribed code
function code_to_perm(code::Vector{<:Integer})
	isempty(code) && return Permutation()
	
	w = Int8[]
	
	k = length(code)
	n = maximum([i+code[i] for i in 1:k])
	available = collect(1:n)
	
	for i in 1:length(code)
		push!(w,available[1+code[i]])
		deleteat!(available,1+code[i])
	end
	
	return Permutation(vcat(w,available),1)
end
function code_to_perm(code::Integer...) return code_to_perm([i for i in code]) end


# returns the descent set of w
function descents(w::Permutation)
	Des = Int8[]
	for i in min!(w):max!(w)
		if w[i] > w[i+1] push!(Des,i) end
	end
	return Des
end
# implemented in "permtools.jl"
#function descents(w::Vector{<:Integer}) return descents(Permutation(w)) end

# returns the major index of w
function maj(w::Permutation) return sum(descents(w)) end
function maj(w::Vector{<:Integer}) return maj(Permutation(w)) end



#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--

	julia> w = Permutation(1,3,5,2,8,7,6,4);

	julia> println(support(w))
	Int8[2, 3, 4, 5, 6, 7, 8]
	
	julia> inverse(w)
	Permutation(1, 4, 2, 8, 3, 7, 6, 5)

	julia> w*inverse(w) == Permutation() && inverse(w)*w == Permutation()
	true
	
	julia> println(cycle_decomposition(w))
	Vector{Int8}[[2, 3, 5, 8, 4], [6, 7]]

	julia> order(w)
	10
	
	julia> w*w*w*w*w*w*w*w*w*w == Permutation()
	true

	julia> println(inversions(w))
	Tuple{Int8, Int8}[(2, 4), (3, 4), (3, 8), (5, 6), (5, 7), (5, 8), (6, 7), (6, 8), (7, 8)]

	julia> println(code(w))
	Int8[0, 1, 2, 0, 3, 2, 1]

	julia> len(w)
	9

	julia> code_to_perm([0,1,2,0,3,2,1])
	Permutation(1, 3, 5, 2, 8, 7, 6, 4)

	julia> println(descents(w))
	Int8[3, 5, 6, 7]
	
	julia> maj(w)
	21

+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#




		# * * Reduced Words * * #
		
export coxeter_product, demazure_product, is_reduced

function coxeter_product(word::Vector{<:Integer})
	w = Permutation()
	for i in word w = w*s_i(i) end
	return w
end
function coxeter_product(word::Integer...) return coxeter_product([i for i in word]) end

function demazure_product(word::Vector{<:Integer})
	w = Permutation()
	for i in word if w[i] < w[i+1] w = w*s_i(i) end end
	return w
end
function demazure_product(word::Integer...) return demazure_product([i for i in word]) end

function is_reduced(word::Vector{<:Integer})
	w = coxeter_product(word)
	return len(w) == length(word)
end
function is_reduced(word::Integer...) return is_reduced([i for i in word]) end

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
		
# IMPORTANT NOTE:
# Permutations w in this section must have positive support, and are only considered to be permutations of postitive support.
# For instance, w = [2,1] is regarded as not containing a 132-pattern, even though it is technically w = [...,0,2,1,...]
# Likewise, patterns u in this section must be permutations of the numbers 1,2,...,k with u[k] != k so that they are well-defined as S_\infty statistics.
# Continuing the example above, w = [2,1,3,4,5,...] contains infinitely many 213-patterns, namely at the positions [1,2,3], [1,2,4], [1,2,5], ...

export patterns, max_pattern, max_inversion, max_132, max_2143, is_dominant, is_vexillary


# finds all occurrances of u as a pattern in w and returns the vector of indices (or [] if w avoids u)
function patterns(w::Permutation,u::Vector{<:Integer})
	# implement a check to see if u is a good pattern? or just require u to be a permutation?
	n = max!(w)
	k = length(u)
	
	# n<=k cases
	if n<k return Int8[]
	elseif n==k
		w[1:k]==u && return collect(1:k)
		return Int8[]
	end
	
	pats = Vector{Int8}[]
	for subset in kSubsets(n,k) if invperm(sortperm(w[subset]))==u push!(pats,subset) end end
	return reverse(pats)
end
function patterns(w::Vector{<:Integer},u::Vector{<:Integer}) return patterns(Permutation(w),u) end

# finds the lex-final occurance of u as a pattern in w and returns its indices (or [] if w avoids u)
function max_pattern(w::Permutation,u::Vector{<:Integer})
	# implement a check to see if u is a good pattern? or just require u to be a permutation?
	n = max!(w)
	k = length(u)
	
	# n<=k cases
	if n<k return Int8[]
	elseif n==k
		w[1:k]==u && return collect(1:k)
		return Int8[]
	end
	
	for subset in kSubsets(n,k) if invperm(sortperm(w[subset]))==u return subset end end
	return Int8[]
end
function max_pattern(w::Vector{<:Integer},u::Vector{<:Integer}) return max_pattern(Permutation(w),u) end

# Specific instances of max_pattern(w,u) for the patterns u = [2,1], [1,3,2], and [2,1,4,3].
function max_inversion(w::Permutation) return max_pattern(w,[2,1]) end
function max_inversion(w::Vector{<:Integer}) return max_pattern(w,[2,1]) end
function max_132(w::Permutation) return max_pattern(w,[1,3,2]) end
function max_132(w::Vector{<:Integer}) return max_pattern(w,[1,3,2]) end
function max_2143(w::Permutation) return max_pattern(w,[2,1,4,3]) end
function max_2143(w::Vector{<:Integer}) return max_pattern(w,[2,1,4,3]) end

function is_dominant(w::Permutation) return max_132(w) == Any[] end
function is_dominant(w::Vector{<:Integer}) return is_dominant(Permutation(w)) end
function is_vexillary(w::Permutation) return max_2143(w) == Any[] end
function is_vexillary(w::Vector{<:Integer}) return is_vexillary(Permutation(w)) end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--
	
	julia> w = Permutation(2,5,9,6,1,7,3,4,6,8,10);
	
	julia> println(max_inversion(w))
	Int8[6, 9]
	
	julia> println(max_132(w))
	Int8[5, 6, 9]

	julia> println(max_2143(w))
	Int8[4, 5, 6, 9]
	
	julia> println(max_pattern(w,[1,4,3,2]))
	Int8[2, 3, 6, 9]
	
	# w avoids the pattern [4,3,2,1]
	julia> println(max_pattern(w,[4,3,2,1]))
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
	
	julia> w = Permutation(8,6,3,2,1,4,5,7); # w is dominant
	
	julia> dominant_transition(w)
	Any[]
	
	julia> max_transition(w)
	Any[]
	
	
	julia> u = Permutation(2,8,3,6,1,4,5,7); # u is vexillary, but not dominant
	
	julia> dominant_transition(u)
	3-element Vector{Any}:
	 (4, 7)
	 Permutation(2, 8, 3, 5, 1, 4, 6, 7)
	 Permutation[Permutation(2, 8, 5, 3, 1, 4, 6, 7)]
	
	julia> max_transition(u)
	Any[]
	
	
	julia> v = Permutation(2,1,5,4,3,8,6,7); # v is not vexillary

	julia> dominant_transition(v)
	3-element Vector{Any}:
	 (6, 8)
	 Permutation(2, 1, 5, 4, 3, 7, 6)
	 Permutation[Permutation(2, 1, 7, 4, 3, 5, 6), Permutation(2, 1, 5, 7, 3, 4, 6), Permutation(2, 1, 5, 4, 7, 3, 6)]

	julia> max_transition(v) == dominant_transition(v)
	true
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =# 



