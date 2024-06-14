using SchubertPolynomials
using Test


if isempty(ARGS)
   tests = [ 
      "bpds.jl", 
      "drifts.jl", 
      "permtools.jl", 
      "double_poly_ring.jl", 
      "ssyt.jl", 
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