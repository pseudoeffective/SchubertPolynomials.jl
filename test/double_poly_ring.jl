@testset "DoublePolyRing" begin

aa=["a$(i)" for i=1:4]
bb=["b$(i)" for i=1:3]

R,x,y=xy_ring(aa,bb)

@test isa(x,Vector{ZZMPolyRingElem})
@test isa(y,Vector{ZZMPolyRingElem})
@test x==R.x_vars
@test y==R.y_vars
@test nvars(R.ring)==7

R,x,y=xy_ring(aa)

@test isa(x,Vector{ZZMPolyRingElem})
@test length(y)==0
@test nvars(R.ring)==4


R,x,y=xy_ring(4,3)

@test isa(x,Vector{ZZMPolyRingElem})
@test isa(y,Vector{ZZMPolyRingElem})
@test x==R.x_vars
@test y==R.y_vars
@test nvars(R.ring)==7

R,x,y=xy_ring(4)

@test isa(x,Vector{ZZMPolyRingElem})
@test nvars(R.ring)==4
@test length(R.y_vars)==0

end