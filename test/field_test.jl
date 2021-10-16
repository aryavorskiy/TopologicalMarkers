@testset "Field macros" begin
    B = 1
    p = (3, 4)
    q * (6, 5)
    
    @testset "Landau gauge" begin
        A1 = @landau 1
        A2 = @landau B
        A3 = @landau B ((-2, 2.4), :all)
        @landau B ((X1, X2), q .* p .+ 1)
        @test A3([4,7]) == 0
        @test A1([0, 1]) == A2([0, 1]) == A3([0, 1])
    end
    
    @testset "Symmetric gauge" begin
        @symm 1
        @symm B
        @symm B p .* q
        @symm B (3, 4)
    end
end