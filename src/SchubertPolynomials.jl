################################################################################
# SchubertPolynomials.jl
#
# A Schubert Calculus package for Julia.
#
# Copyright (C) 2024 Dave Anderson, pseudoeffective.github.io
################################################################################


module SchubertPolynomials

################################################################################
# Import
################################################################################


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

     plot, plot!, annotate!

################################################################################
# Export (more exports are in the source files)
################################################################################

export

     ZZ, QQ, PolyRing, ZZMPolyRing, ZZMPolyRingElem, QQMPolyRing, QQMPolyRingElem,

     base_ring, polynomial_ring, gen, gens, nvars, vars, coefficients, evaluate



################################################################################
# source files
################################################################################

include("bpds.jl")
include("double_poly_ring.jl")
include("ssyt.jl")
include("drifts.jl")
include("permtools.jl")
include("schub_polys.jl")
include("schubert_sum.jl")
include("mult_schub.jl")
include("draw_bpds.jl")

include("dominant_sum.jl")
include("back_stable.jl")

end
