@testset "Hamiltonian test" begin
    ms = ones(3, 3)
    a = CoordinateRepr(ms)
    a = CoordinateRepr(ms, :c)
    a = CoordinateRepr(ms, :n)
    h = hamiltonian(a)
    field!(h, x -> [0, x[1], 0])
    zmap = ones(3, 3)
    zmap[2,2] = 2
    zones!(h, CoordinateRepr(zmap))
end

@testset "Currents test" begin
    p_b = filled_projector(ham_b)
    curr_b = currents(ham_b, p_b) / B
    @test all(abs(sum(curr_b[i, :])) < 1e-10 for i in 1:size(curr_b)[2])
end