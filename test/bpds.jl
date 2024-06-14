@testset "BPDs" begin

# constructors

mtx = Matrix( [ 0 0 2 ; 0 2 1 ; 2 1 1 ] )
b = BPD( mtx )

@test b.m == mtx

mtx2 = Matrix( [ "O" "O" "/"; "O" "/" "+"; "/" "+" "+" ] )
b2 = BPD( mtx2 )

@test b==b2

bw = Rothe([3,2,1])

@test b==bw


# iterators

w=[2,1,4,3]
bps=collect(all_bpds(w))
fbps=collect(flat_bpds(w))

@test bps==fbps

w=[1,3,2,9,8,7,6,5,4]
bps=collect(all_bpds(w))
fbps=collect(flat_bpds(w))

@test length(bps)==163592
@test length(fbps)==21

Kbps=collect(all_Kbpds(w));

@test length(Kbps)==437569

# ASMs

bs = collect(all_bpds([1,3,2]))
a = bpd2asm(bs[2])

@test a == Matrix{Int8}([ 0 1 0; 1 -1 1; 0 1 0 ])

bs = collect(all_bpds([3,1,5,2,4]))
b = bs[4]
a = bpd2asm(b)

@test b == asm2bpd(a)

end
