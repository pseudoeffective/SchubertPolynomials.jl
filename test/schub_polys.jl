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
p2=schub_poly(w,R,method="drift");   # :drift is deprecated and now aliased to :bpd
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


# --- `double` keyword (Step A) ---

w = [2,1,4,3]
Rxy = xy_ring(3,3)[1]   # ring carrying y-variables

# double=false: single polynomial, even though the ring has y-variables (new capability)
p_single = schub_poly(w, Rxy; double=false)
@test p_single == schub_poly(w, Rxy; double=false, method="bpd")
@test p_single == schub_poly(w, Rxy; double=false, method="transition")
@test p_single == schub_poly(w, Rxy; double=false, method="dd")
# the y-variables of the ring are genuinely unused
@test all( v -> !startswith(string(v), "y"), vars(p_single) )

# double=true: factorial polynomial in x and y, methods agree
p_double = schub_poly(w, Rxy; double=true)
@test p_double == schub_poly(w, Rxy; double=true, method="bpd")
@test p_double == schub_poly(w, Rxy; double=true, method="transition")
@test p_double != p_single
@test any( v -> startswith(string(v), "y"), vars(p_double) )

# double=true on a y-free ring is an error
@test_throws ArgumentError schub_poly(w, xy_ring(3,0)[1]; double=true)
@test_throws ArgumentError groth_poly(w, xy_ring(3,0)[1]; double=true)

# no ring supplied: double=true builds a ring with y-variables on its own
@test any( v -> startswith(string(v), "y"), vars( schub_poly(w; double=true) ) )


# --- dominant permutations: transition base case is drift-free (Step C2) ---
# For a dominant w, max_transition(w) is empty, so :transition falls back to the
# single flat (Rothe) BPD; this must agree with :bpd in both single and double.

for wd in ( [3,2,1], [1,3,4,2] )
  Rd = xy_ring(4,4)[1]
  @test schub_poly(wd, Rd; double=false, method="transition") ==
        schub_poly(wd, Rd; double=false, method="bpd")
  @test schub_poly(wd, Rd; double=true,  method="transition") ==
        schub_poly(wd, Rd; double=true,  method="bpd")
  # memoized transition path hits the same fallback
  @test schub_poly(wd, Rd; double=true, memo=true) ==
        schub_poly(wd, Rd; double=true, method="bpd")
end

end