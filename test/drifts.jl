@testset "Drifts" begin

# constructors

mtx = Matrix( [ 0 0 2 ; 0 2 1 ; 2 1 1 ] )
mtxd = Matrix( [ 0 0 8 ; 0 8 1 ; 8 1 1 ] )
b = BPD( mtx )
d = bpd2drift(b)

@test d.m == mtxd

mtx2 = Matrix( [ "O" "O" ""; "O" "" "+"; "" "+" "+" ] )
d2 = Drift( mtx2 )

@test d==d2



# iterators

w=[1,2,7,5,3,4,6]
bps=collect(all_bpds(w))
b=bps[57]
d=bpd2drift(b)
dc = collect(drift_class(d))

@test length(dc)==66

w=[1,3,2,9,8,7,6,5,4]
fbs=collect(flat_bpds(w))
b=fbs[21]
d=bpd2drift(b)
dc=collect(drift_class(d))

@test length(dc)==24453

la=[3,2,1]
ff=[2,4,5]
d=partition2drift(la,ff)
dc=collect(drift_class(d))

@test length(dc)==44


end
