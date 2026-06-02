@testset "SSYT" begin


R = schub_ring(5)
x = extract_vars(R; varname=:x)

la=[2,1]
ff=2
p = schur_poly(la,ff; ring=R)

@test p==x[1]^2*x[2]+x[1]*x[2]^2

p2 = schur_poly(la,ff; ring=R, mu=[1])

@test p2==x[1]^2+2*x[1]*x[2]+x[2]^2

p3 = schur_poly(la,ff; ring=R, mu=[1], rowmin=true)

@test p3==x[1]*x[2]+x[2]^2

R = schub_ring(2,2)
x = extract_vars(R; varname=:x)
y = extract_vars(R; varname=:y)
la=[2,1]
ff=2

# single vs double is controlled by the `double` keyword (not the ring's shape)
@test schur_poly(la,ff; ring=R, double=false) == x[1]^2*x[2] + x[1]*x[2]^2
@test schur_poly(la,ff; ring=R) == x[1]^2*x[2] + x[1]*x[2]^2   # default is single

p = schur_poly(la,ff; ring=R, double=true)   # double/factorial Schur polynomial
@test p == x[1]^2*x[2] + x[1]^2*y[1] + x[1]*x[2]^2 + 2*x[1]*x[2]*y[1] + x[1]*x[2]*y[2] + x[1]*y[1]^2 + x[1]*y[1]*y[2] + x[2]^2*y[1] + x[2]*y[1]^2 + x[2]*y[1]*y[2] + y[1]^2*y[2]

la = [3,2,1]
ff = [2,4,5]
R = schub_ring(5)
p = schur_poly(la,ff; ring=R)

@test length( ssyt(la,ff) ) == 44
@test sum(coefficients(p))==44

la = [4,4,2,1]
ff = [[3,4,4,5],[3,4,5,6],[4,5],[6]]
R = schub_ring(6)
p = schur_poly(la,ff; ring=R);
tt = ssyt(la,ff);

@test sum(coefficients(p))==3693
@test length(tt)==3693

end
