@testset "Schubert Polynomials" begin


w=[2,1,4,3]
R,x,y=xy_ring(3)
p = schub_poly(w,R)
q = groth_poly(w,R)

@test p == x[1]*(x[1]+x[2]+x[3])
@test q == x[1]*(x[1]+x[2]+x[3]-x[1]*x[2]-x[1]*x[3]-x[2]*x[3]+x[1]*x[2]*x[3])

w = [1,4,2,3,7,6,5]
R=xy_ring(6,6)[1]
p1=schub_poly(w,R,method="bpd");
p2=schub_poly(w,R,method="drift");
p3=schub_poly(w,R,method="transition");

@test p1==p2
@test p2==p3

R=xy_ring(6)[1]
p1=schub_poly(w,R,method="bpd");
p2=schub_poly(w,R,method="dd");

@test p1==p2


w=[1,4,3,2,10,9,8,7,6,5];
@test nschub(w)==4424420

w=[1,3,2,8,7,6,5,4];
@test ngroth(w)==1711251

end