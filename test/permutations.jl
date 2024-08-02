
#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
    
	test\permutations.jl
    Hugh Dennin, 2 August 2024.
	
    Testing environment for `permutations.jl`.
    
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#


@testset "Permutations" begin

	@testset "Constructors" begin
		w = Permutation([1,3,5,2,4])
		@test w[1:5] == [1,3,5,2,4]
		u = Permutation([-1,0,4,-2,1,2,3])
		@test u[-2:5] == [-1,0,4,-2,1,2,3,5]
	
		w = Permutation(1,3,5,2,4)
		@test w[1:5] == [1,3,5,2,4]
		u = Permutation(-1,0,4,-2,1,2,3)
		u[-2:5] == [-1,0,4,-2,1,2,3,5]

		v = Permutation([(1,1),(2,4),(3,2),(4,3),(5,0),(0,5)])
		@test v[0:8] == [5,1,4,2,3,0,6,7,8]
		
		v = Permutation((1,1),(2,4),(3,2),(4,3),(5,0),(0,5))
		@test v[0:8] == [5,1,4,2,3,0,6,7,8]
		
		id = Permutation()
		@test id[-10:10] == collect(-10:10)
	end
	
	@testset "Equality" begin
		w = Permutation(4,3,6,5)
		u = Permutation(2,4,3,6,5,7,8,9)
		@test w==u
		
		v = Permutation([(3,4),(4,3),(6,5),(5,6)])
		@test w==v
		
		wvec = [4,3,6,5]
		@test w != wvec
	end
	
	@testset "Mutators" begin
		w = Permutation(1,2,5,8,6,4,3,7,9,10,11)
		id = Permutation()
		
		@test min!(w) == 3
		@test w == Permutation(1,2,5,8,6,4,3,7,9,10,11)
		@test min!(id) == 1
		
		@test max!(w) == 8
		@test w == Permutation(1,2,5,8,6,4,3,7,9,10,11)
		@test max!(id) == 1
	end
	
	@testset "Iterators" begin
		S3 = [
			Permutation(1),
			Permutation(2,1),
			Permutation(3,1,2),
			Permutation(1,3,2),
			Permutation(2,3,1),
			Permutation(3,2,1)
		]
		@test collect(Permutations(3)) == S3
		@test collect(Permutations(0)) == [Permutation()]
		@test isempty(collect(Permutations(-3)))
		@test length(Permutations(5)) == 120
		@test eltype(Permutations(5)) == Permutation
	end
	
	
	@testset "Special Permutations" begin
		c = cycle([1,4,3])
		@test c == Permutation(4,2,1,3)
		c = cycle(1,4,3)
		@test c == Permutation(4,2,1,3)
		c = cycle(20)
		@test c == Permutation()
		
		t = t_ij(2,5)
		@test t == Permutation(1,5,3,4,2)
		t = t_ij(5,2)
		@test t == Permutation(1,5,3,4,2)
		t = t_ij(3,3)
		@test t == Permutation()
		
		s = s_i(3)
		@test s == Permutation(1,2,4,3)
		
		w0 = longest_perm(6)
		@test w0 == Permutation(6,5,4,3,2,1)
	end
	
	@testset "Group Operations" begin
		w = Permutation(1,5,3,7,2,6,8,4)
		u = Permutation(2,4,6,8,1,3,5,7)
		
		@test w*u == Permutation(5,7,6,4,1,3,2)
		@test u*w == Permutation(2,1,6,5,4,3)
		@test w*w == Permutation(1,2,3,8,5,6,4,7)
		@test w*w*w*w*w*w == Permutation()
		
		@test inverse(w) == Permutation(1,5,3,8,2,6,4,7)
		@test conjugate(w,u) == Permutation(5,1,6,8,7,3,4,2)
		@test cycle_decomposition(w) == [[2,5],[4,7,8]]
		@test cycle_decomposition(u) == [[1,2,4,8,7,5],[3,6]]
		@test cycle_decomposition(conjugate(w,u)) == [[1,5,7,4,8,2],[3,6]]
		@test order(w) == 6
		@test order(u) == 6
		@test support(w) == [2,4,5,7,8]
		@test support(u) == [1,2,3,4,5,6,7,8]
	end
	
	@testset "Reduced Words" begin
		@test coxeter_product(3,3) == Permutation()
		@test coxeter_product(3,2,3,1,3,2) == Permutation(4,2,1,3)
		@test demazure_product(3,3) == Permutation(1,2,4,3)
		@test demazure_product(3,2,3,1,3,2) == Permutation(4,3,1,2)
		
		@test !is_reduced(3,2,3,1,3,2)
		@test is_reduced(3,2,3,1)
	end
	
	@testset "Permutation Statistics" begin
		w = Permutation(1,3,5,2,7,6,4)
		u = Permutation(2,1,4,5,3)
		
		@test inversions(w) == [(2, 4), (3, 4), (3, 7), (5, 6), (5, 7), (6, 7)]
		@test len(w) == 6
		@test is_even(w)
		@test len(u) == 3
		@test !is_even(u)
		
		@test code(w) == [0,1,2,0,2,1]
		@test code_to_perm([0,1,2,0,2,1]) == Permutation(1,3,5,2,7,6,4)
		
		@test descents(w) == [3,5,6]
		@test maj(w) == 14
	end
	
	@testset "Pattern Containment" begin
		w = Permutation(2,5,9,6,1,7,3,4,6,8,10)
		
		@test max_pattern(w,[1,4,3,2]) == [2,3,6,9]
		@test max_pattern(w,[4,3,2,1]) == []
		@test max_inversion(w) == [6,9]
		@test max_132(w) == [5,6,9]
		@test max_2143(w) == [4,5,6,9]
		
		w = Permutation(8,6,3,2,1,4,5,7) # dominant
		u = Permutation(2,8,3,6,1,4,5,7) # vexillary
		v = Permutation(2,1,5,4,3,8,6,7)
		
		@test is_dominant(w)
		@test !is_dominant(u)
		@test !is_dominant(v)
		@test is_vexillary(w)
		@test is_vexillary(u)
		@test !is_vexillary(v)
	end
	
	@testset "Transition" begin
		w = Permutation(8,6,3,2,1,4,5,7) # dominant
		u = Permutation(2,8,3,6,1,4,5,7) # vexillary
		v = Permutation(2,1,5,4,3,8,6,7)
		
		@test dominant_transition(w) == []
		@test dominant_transition(u) == [
			(4,7),
			Permutation(2,8,3,5,1,4,6,7),
			[Permutation(2,8,5,3,1,4,6,7)]
		]
		@test dominant_transition(v) == [
			(6,8),
			Permutation(2,1,5,4,3,7,6),
			[Permutation(2,1,7,4,3,5,6),Permutation(2,1,5,7,3,4,6),Permutation(2,1,5,4,7,3,6)]
		]
		
		@test max_transition(w) == []
		@test max_transition(u) == []
		@test max_transition(v) == [
			(6,8),
			Permutation(2,1,5,4,3,7,6),
			[Permutation(2,1,7,4,3,5,6),Permutation(2,1,5,7,3,4,6),Permutation(2,1,5,4,7,3,6)]
		]
	end
end
