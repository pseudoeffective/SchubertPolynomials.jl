
#==============================================================================================================================+
	
	subsets.jl
	Hugh Dennin, 1 August 2024.

	This file contains the definition of a iterator for k-element subsets of an n-element set {1,2,...,n}.
	n and k can take any Int8 value

+==============================================================================================================================#

export kSubsets

#==============================================================================================================================#

# iterator for k-element subsets of the set {1,2,...,n}
# subsets are iterated in rev-lex order
struct kSubsets
	n::Int8
	k::Int8
end

function Base.iterate(iter::kSubsets, state=Int8[iter.n-i+1 for i in iter.k:-1:1])
	# k<=0 cases
	if isempty(state)
		iter.k==0 && return ([],[0])
		return nothing
	end
	
	# end of iteration
    state[1] <= 0 && return nothing
	
    # next state
	new_state = copy(state)
	i = 1
	while i != iter.k && new_state[iter.k-i+1] == new_state[iter.k-i]+1
		new_state[iter.k-i+1] = iter.n-i+1
		i += 1
	end
	new_state[iter.k-i+1] -= 1
	
	return (state,new_state)
end

Base.length(iter::kSubsets) = iter.n >= iter.k ? binomial(iter.n,iter.k) : 0
Base.eltype(::Type{kSubsets}) = Vector{Int8}

#==============================================================================================================================+

	--Example--
	
	julia> collect(kSubsets(4,2))
	6-element Vector{Vector{Int8}}:
	 [3, 4]
	 [2, 4]
	 [2, 3]
	 [1, 4]
	 [1, 3]
	 [1, 2]
	
	julia> length(kSubsets(5,3))
	10
	
	julia> collect(kSubsets(3,0))
	1-element Vector{Vector{Int8}}
	 []
	
	julia> collect(kSubsets(3,-1))
	Vector{Int8}[]
	
+==============================================================================================================================#

