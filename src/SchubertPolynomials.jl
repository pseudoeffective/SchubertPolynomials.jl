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


#=
# Plots
import Plots:
	plot, plot!, annotate!, savefig

# Printf 
import Printf:
	@sprintf
=#

# bpds
using BumplessPipeDreams

# Tableaux and Schur polynomials live in SemistandardTableaux now; import the
# names we re-export selectively (a blanket `using` would clash on `xy_ring`,
# `ZZ`, `QQ`, `schur_poly`, etc. — see the deprecation shims those packages keep).
import SemistandardTableaux
import SemistandardTableaux:
	Tableau, ssyt, ssyt_iterator, schur_poly, rsk, rsk_insert,
	edelman_greene, edelman_greene_insert

################################################################################
# Export (more exports are in the source files)
################################################################################

export
	ZZ, QQ, PolyRing, ZZMPolyRing, ZZMPolyRingElem, QQMPolyRing, QQMPolyRingElem,

	base_ring, polynomial_ring, gen, gens, nvars, vars, coefficients, evaluate,

	BPD, all_bpds, all_Kbpds, flat_bpds, Rothe, is_asm, bpd2asm, asm2bpd, bpd2word,

	bpd2perm, isreduced, isflat,

	# re-exported from SemistandardTableaux
	Tableau, ssyt, ssyt_iterator, schur_poly, rsk, rsk_insert,
	edelman_greene, edelman_greene_insert



################################################################################
# source files
################################################################################

include("subsets.jl")
include("double_poly_ring.jl")
include("permtools.jl")
include("permutations.jl") # meant to replace `permtools.jl`, missing some functionality; must be included after `subsets.jl`
# Tableaux/Schur moved to SemistandardTableaux.jl (imported above); ssyt.jl is
# kept on disk for one version but no longer included.
#include("ssyt.jl") # must be included after `double_poly_ring.jl`
include("pds.jl") # must be included after `permutations.jl`
#include("bpds.jl")
# Drift methods have moved to DriftPolynomials.jl.  The drift source files are
# kept for one version but no longer included; `compat.jl` provides migration
# shims for the names that used to live here.
#include("drifts.jl")
#include("drift_polys.jl") # must be included after `drifts.jl`
#include("draw_bpds.jl") # must be included after `bpds.jl` & `drifts.jl`
include("schub_polys.jl")
include("schubert_sum.jl")
include("mult_schub.jl")

include("dominant_sum.jl")
include("back_stable.jl")

include("specialize.jl")

include("compat.jl")

end
