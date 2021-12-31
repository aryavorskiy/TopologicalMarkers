@testset "Hamiltonian test" begin
    ms = ones(3, 3)
    a = CoordinateRepr(ms)
    a = CoordinateRepr(ms, :c)
    a = CoordinateRepr(ms, :n)
    h = hamiltonian(a)
    field!(h, x -> [0, x[1]])
    zmap = fill(:ext, 3, 3)
    zmap[2,2] = :Int
    domains!(h, CoordinateRepr(zmap))
    domains!(h, zmap, :n)

    siz = (3, 3)
    @test hamiltonian(2, siz) == hamiltonian(2ones(siz), :c) == hamiltonian(CoordinateRepr(2ones(siz)))
    @test hamiltonian(2, siz) == hamiltonian(2, siz, pbc = (false, false))
    @test hamiltonian(2, siz) != hamiltonian(2, siz, pbc = (true, false))
    @test hamiltonian(2, siz) != hamiltonian(2, siz, pbc = (false, true))
    @test hamiltonian(2, siz) != hamiltonian(2, siz, pbc = (true, true))
end

@testset "Currents test" begin
    B = 1e-8
    ms = CoordinateRepr(ones(15, 15) * -1)
    ham_b = hamiltonian(ms, field=@symm(B))
    p_b = filled_projector(ham_b)
    curr_b = currents(ham_b, p_b)
    @test [abs(sum(curr_b[i, :])) for i in 1:size(curr_b)[2]] |> maximum < 1e-12
end
