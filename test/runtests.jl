using SchubertPolynomials
using Test


if isempty(ARGS)
    tests = [ 
		"subsets.jl",
		"permutations.jl",
        "permtools.jl", 
        "double_poly_ring.jl", 
        "ssyt.jl",
        "bpds.jl", 
        "drifts.jl",
        "schub_polys.jl",
        "schubert_sum.jl",
        "mult_schub.jl"
	]
else
    tests = ARGS
end

for test in tests
    include(test)
end

#=
@testset "SchubertPolynomials.jl" begin
    # Write your tests here.
end
=#