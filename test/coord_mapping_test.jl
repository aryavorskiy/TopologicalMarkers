using StaticArrays

@testset "Coordinate mapping" begin
    size = (12, 15)
    @test adjacent_sites([5, 5], size) == Set([[5,4], [5,6], [4,5], [6,5]])
    @test adjacent_sites([1, 1], size) == Set([[1,2], [2,1]])
    @test adjacent_sites([12, 15], size) == Set([[12,14], [11,15]])
    @test adjacent_sites([12, 15], (12, 16)) == Set([[12,14], [12,16], [11,15]])
    @test adjacent_sites([12, 15], (20, 20)) == Set([[12,14], [12,16], [11,15], [13,15]])
    @test adjacent_sites([1, 1], size, order = 2) == Set([[1,3], [3,1], [2,2]])

    # Type stability check
    @test adjacent_sites([1,1], size) isa Set{Vector{Int}}
    @test adjacent_sites(SA[1,1], size) isa Set{SVector{2, Int}}
end
