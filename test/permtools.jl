@testset "PermutationTools" begin


w=[1, 4, 2, 8, 7, 5, 3, 6]

@test len(w)==10
@test descents(w)==[2,4,5,6]
@test sij(w,2,5)==[1,7,2,8,4,5,3,6]

trn = dominant_transition(w)

@test trn[1]==(6,7)
@test trn[2]==[1, 4, 2, 8, 7, 3, 5, 6]
@test trn[3]==[[1, 4, 3, 8, 7, 2, 5, 6]]

trn = max_transition(w)

@test trn[1]==(5, 8)
@test trn[2]==[1, 4, 2, 8, 6, 5, 3, 7]
@test trn[3]==[[1, 6, 2, 8, 4, 5, 3, 7], [1, 4, 6, 8, 2, 5, 3, 7]]


@test trimw([1,3,2,4])==[1,3,2]
@test trimw([1,2,3])==[]

end