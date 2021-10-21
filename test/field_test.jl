@testset "Field macros" begin
    B = 1
    p = (3, 4)
    q = (6, 5)
    x = [2,3,4]
    y = [6,2,5]
    
    @testset "Landau gauge" begin
        A1 = @landau 1
        A2 = @landau B
        @test A1(x) == A2(x)
        @test A1(y) == A2(y)
    end
    
    @testset "Symmetric gauge" begin
        A1 = @symm 1
        A2 = @symm B
        A3 = @symm B p .* q
        A4 = @symm B (18, 20)
        @test A1(x) == A2(x)
        @test A1(y) == A2(y)
        @test A3(x) == A4(x)
        @test A3(y) == A4(y)

    end

    @testset "Flux" begin
        A1 = @flux 1
        A2 = @flux B
        A3 = @flux B p .* q
        A4 = @flux B (18, 20)
        @test A1(x) == A2(x)
        @test A1(y) == A2(y)
        @test A3(x) == A4(x)
        @test A3(y) == A4(y)
    end
end