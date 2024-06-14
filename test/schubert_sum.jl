@testset "SchubertSum" begin


ee=[1,2,3]
uu=[2,1,3]
ww=[2,3,1]

sse = SchubertSum(ee)
ssu = SchubertSum(uu)
ssw = SchubertSum(ww)

@test sse==SchubertSum([1])
@test ssu.schubs==[[2,1]]

ss = ssu+2*ssw
ss2 = 2*ssw+ssu

@test ss==ss2
@test ss-2*ssw==ssu

end

