
#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	
	test\subsets.jl
	Hugh Dennin, 2 August 2024.
	
	Testing environment for `subsets.jl`
	
+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#


@testset "Subsets" begin
	@test collect(kSubsets(4,2)) == [
		[3,4],
		[2,4],
		[2,3],
		[1,4],
		[1,3],
		[1,2]
	]
	@test length(kSubsets(5,3)) == 10
	
	@test collect(kSubsets(3,0)) == [Int8[]]
	@test collect(kSubsets(3,3)) == [[1,2,3]]
	
	@test isempty(collect(kSubsets(3,-1)))
end
