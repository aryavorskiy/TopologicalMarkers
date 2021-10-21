@testset "Hamiltonian test" begin
    ms = ones(3, 3)
    a = CoordinateRepr(ms)
    a = CoordinateRepr(ms, :c)
    a = CoordinateRepr(ms, :n)
    h = hamiltonian(a)
    field!(h, x -> [0, x[1], 0])
    zmap = fill(:ext, 3, 3)
    zmap[2,2] = :Int
    zones!(h, CoordinateRepr(zmap))
    zones!(h, zmap, :n)
end

@testset "Currents test" begin
    B = 1e-8
    ms = CoordinateRepr(ones(15, 15) * -1)
    ham_b = hamiltonian(ms, field=@symm(B))
    p_b = filled_projector(ham_b)
    curr_b = currents(ham_b, p_b)
    @test [abs(sum(curr_b[i, :])) for i in 1:size(curr_b)[2]] |> maximum < 1e-12
end