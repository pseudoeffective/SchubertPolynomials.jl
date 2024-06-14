@testset "SSYT" begin


R,x,y=xy_ring(5)

la=[2,1]
ff=2
p = schur_poly(la,ff,R)

@test p==x[1]^2*x[2]+x[1]*x[2]^2

p2 = schur_poly(la,ff,R; mu=[1])

@test p2==x[1]^2+2*x[1]*x[2]+x[2]^2

p3 = schur_poly(la,ff,R; mu=[1],rowmin=true)

@test p3==x[1]*x[2]+x[2]^2

R,x,y=xy_ring(2,2)
la=[2,1]
ff=2
p = schur_poly(la,ff,R)

@test p == x[1]^2*x[2] + x[1]^2*y[1] + x[1]*x[2]^2 + 2*x[1]*x[2]*y[1] + x[1]*x[2]*y[2] + x[1]*y[1]^2 + x[1]*y[1]*y[2] + x[2]^2*y[1] + x[2]*y[1]^2 + x[2]*y[1]*y[2] + y[1]^2*y[2]

la = [3,2,1]
ff = [2,4,5]
R = xy_ring(5)[1]
p = schur_poly(la,ff,R)

@test length( ssyt(la,ff) ) == 44
@test sum(coefficients(p))==44

la = [4,4,2,1]
ff = [[3,4,4,5],[3,4,5,6],[4,5],[6]]
R = xy_ring(6)[1]
p = schur_poly(la,ff,R);

@test sum(coefficients(p))==3693
@test length(ssyt(la,ff))==3693

end