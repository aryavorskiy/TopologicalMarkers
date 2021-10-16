@testset "Hamiltonian tests" begin
    ms = ones(3, 3)
    a = CoordinateRepr(ms)
    a = CoordinateRepr(ms, :c)
    a = CoordinateRepr(ms, :n)
    h = hamiltonian(a)
    field!(h, x -> [0, x[1], 0])
    zmap = ones(3, 3)
    zmap[2,2] = 2
    zones!(h, zmap)
end