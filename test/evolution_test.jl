using LinearAlgebra

@testset "Evolution operator" begin
    H = hamiltonian(ones(11, 12), :c, field=@landau(1))
    t = 1e-2
    ev_ex = exp(im * t * H)
    ev_li = I + im * t * H
    @test evolution_operator(H, t) == ev_ex
    TopologicalMarkers._configure_exp!(false, cache=false)
    @test evolution_operator(H, t) == ev_ex
    TopologicalMarkers._configure_exp!(true, order=1, threshold=1, cache=false)
    @test evolution_operator(H, t) == ev_li
    TopologicalMarkers._configure_exp!(true, cache=true)
    @test evolution_operator(H, t) == ev_li
    TopologicalMarkers._configure_exp!(false)
    @test evolution_operator(H, t) == ev_ex
    TopologicalMarkers._configure_exp!(true, threshold=1e-4)
    @test evolution_operator(H, 1e-5) == I + im * 1e-5 * H
    @test evolution_operator(H, 1e-3) == exp(im * 1e-3 * H)
    TopologicalMarkers._configure_exp!(true, threshold=nothing)
    @test evolution_operator(H, 1e-5) == I + im * 1e-5 * H
    @test evolution_operator(H, 1e-3) == I + im * 1e-3 * H
    TopologicalMarkers._configure_exp!(true, order = 2)
    @test evolution_operator(H, 1e-5) == TopologicalMarkers.taylor_exp(im * 1e-5 * H, 2)
    @test evolution_operator(H, 1e-3) == TopologicalMarkers.taylor_exp(im * 1e-3 * H, 2)
    TopologicalMarkers._configure_exp!(true, threshold=1e-5)
    @test evolution_operator(H, 1e-3) == exp(im * 1e-3 * H)
end
