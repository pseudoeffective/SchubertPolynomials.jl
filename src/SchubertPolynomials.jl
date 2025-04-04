################################################################################
# SchubertPolynomials.jl
#
# A Schubert Calculus package for Julia.
#
# Copyright (C) 2025 Dave Anderson, pseudoeffective.github.io
################################################################################


module SchubertPolynomials

################################################################################
# Import
################################################################################

# Base
import Base:
	*, transpose

# Combinatorics
import Combinatorics:
	nthperm

# AbstractAlgebra
import AbstractAlgebra:
	base_ring, gen, gens, parent_type, nvars, polynomial_ring, MPolyBuildCtx, push_term!, finish

# Nemo
import Nemo:
	ZZ, QQ, libflint, ZZMPolyRing, ZZMPolyRingElem, evaluate, vars, coefficients

# LinearAlgebra for determinant
import LinearAlgebra: 
	det

# Memoization
import Memoization: 
	@memoize

# Plots
import Plots:
	plot, plot!, annotate!, savefig



################################################################################
# Export (more exports are in the source files)
################################################################################

export
	ZZ, QQ, PolyRing, ZZMPolyRing, ZZMPolyRingElem, QQMPolyRing, QQMPolyRingElem,
	
	base_ring, polynomial_ring, gen, gens, nvars, vars, coefficients, evaluate



################################################################################
# source files
################################################################################

include("subsets.jl")
include("double_poly_ring.jl")
include("permtools.jl")
include("permutations.jl") # meant to replace `permtools.jl`, missing some functionality; must be included after `subsets.jl`
include("ssyt.jl") # must be included after `double_poly_ring.jl`
include("pds.jl") # must be included after `permutations.jl`
include("bpds.jl")
include("drifts.jl")
include("draw_bpds.jl") # must be included after `bpds.jl` & `drifts.jl`
include("schub_polys.jl")
include("schubert_sum.jl")
include("mult_schub.jl")

include("dominant_sum.jl")
include("back_stable.jl")

include("specialize.jl")

end
