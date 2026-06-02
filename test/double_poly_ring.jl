@testset "schub_ring" begin

# schub_ring returns a bare MPolyRing; variables recovered via extract_vars
R = schub_ring(4,3)

@test nvars(R)==7
@test isa(extract_vars(R; varname=:x), Vector{ZZMPolyRingElem})
@test isa(extract_vars(R; varname=:y), Vector{ZZMPolyRingElem})
@test length(extract_vars(R; varname=:x))==4
@test length(extract_vars(R; varname=:y))==3

R = schub_ring(4)

@test nvars(R)==4
@test length(extract_vars(R; varname=:x))==4
@test length(extract_vars(R; varname=:y))==0

# coefficient ring keyword
@test base_ring(schub_ring(3; coeff=QQ)) == QQ

# custom variable names
Rab = schub_ring(2,2; xname=:a, yname=:b)
@test length(extract_vars(Rab; varname=:a))==2
@test length(extract_vars(Rab; varname=:b))==2

# deprecated DoublePolyRing / xy_ring still functions for one version
R2 = xy_ring(4,3)[1]
@test isa(R2, DoublePolyRing)
@test nvars(R2.ring)==7
@test length(R2.x_vars)==4
@test length(R2.y_vars)==3

end
