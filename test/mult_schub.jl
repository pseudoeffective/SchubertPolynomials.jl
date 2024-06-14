@testset "MultiplySchuberts" begin


#export mult_2schub, lrc, mult_schub, expand_schub





ss = SchubertSum([2,3,1])+SchubertSum([3,1,2])
ss2 = mult_2schub([2,1], [1,3,2])

@test ss==ss2


R,x,y=xy_ring(5,5)

u = [2,4,3,1]
v = [4,1,2,3]
w = [4,2,3,1]

@test lrc( u,v,w,4,R ) == (y[1]-y[4])*(y[3]-y[4])


f = evaluate( schub_poly([2,3,1],R), [y[1],y[2],y[3]],[-x[3],-x[2],-x[1]] )
ss = expand_schub(f,3,R)

@test ss.schubs == [ [1], [1,3,2], [2,3,1] ]
@test ss.coeffs == [ y[1]*y[2] - y[1]*y[3] - y[2]*y[3] + y[3]^2, -2*y[1] + y[2] + y[3], 3 ]



end