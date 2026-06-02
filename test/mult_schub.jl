@testset "MultiplySchuberts" begin


#export mult_2schub, lrc, mult_schub, expand_schub





ss = SchubertSum([2,3,1])+SchubertSum([3,1,2])
ss2 = mult_2schub([2,1], [1,3,2])

@test ss==ss2


R = schub_ring(5,5)
x = extract_vars(R; varname=:x)
y = extract_vars(R; varname=:y)

u = [2,4,3,1]
v = [4,1,2,3]
w = [4,2,3,1]

@test lrc( u,v,w,4; ring=R, double=true ) == (y[1]-y[4])*(y[3]-y[4])


f = evaluate( schub_poly([2,3,1]; ring=R, double=true), [y[1],y[2],y[3]],[-x[3],-x[2],-x[1]] )
ss = expand_schub(f,3; ring=R, double=true)

@test ss.schubs == [ [1], [1,3,2], [2,3,1] ]
@test ss.coeffs == [ y[1]*y[2] - y[1]*y[3] - y[2]*y[3] + y[3]^2, -2*y[1] + y[2] + y[3], 3 ]



end