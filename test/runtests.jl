print("Importing test dependencies... ")

using Test
using TopologicalMarkers
using LinearAlgebra
using Plots

println("done.")

const ENABLE_PROFILE = false
const ENABLE_BENCHMARKS = true

if ENABLE_PROFILE
    include(joinpath("performance", "profile.jl"))
end
if ENABLE_BENCHMARKS
    include(joinpath("performance", "benchmark.jl"))
end

# Macro tests
@testset "Usability tests" begin
    @testset "Static LCM" begin
        global siz = (25, 25)
        global B = 1e-8
        ms = ones(siz) * -1
        ham = hamiltonian(ms, :c)
        ham_b = hamiltonian(CoordinateRepr(ms), field=@symm(B))
        P = filled_projector(ham)
        P_b = filled_projector(ham_b)
        curr_b = currents(ham_b, P_b) / B
        streda = (P - P_b) / B
        plot_auto("Streda" => heatmap_data(streda, siz) => curr_b)

        X, Y = coord_operators()
        plot_auto("Streda" => streda => curr_b,
        "LCM" => 4π * im * P * X * P * Y * P,
        cutaway_view=(8, 8))

        pl = plot()
        currs = @currents @J_treq ham_b P X Y
        plot_figure!(pl, hmap=4π * im * P * X * P * Y * P, currents=currs, xlims=(6,10), ylims=(8, 12))
    end

    @testset "Density currents" begin
        ms = CoordinateRepr(ones(15, 15))
        Bf = 0.01
        τ = 30

        H0 = hamiltonian(ms)
        P0 = filled_projector(H0)
        h(t) = hamiltonian(ms, field=@landau(Bf * t / τ))

        @evolution [
            :ham => h => H,
            P0 => h => P
        ] for t in 0:1:τ
            cur = currents(H, P)
            plot_auto("f" => P => cur * 10, clims=(0.98, 1.02))
        end
    end

    @testset "Adiabatic field on, simplified" begin
        time_domain = 0:0.1:30
        ms = CoordinateRepr(ones(15, 15))
        Bf = 0.01
        τ = 30

        H0 = hamiltonian(ms)
        P0 = filled_projector(H0)
        h(t) = hamiltonian(ms, field=@landau(Bf * t / τ))
        a = Animation()
        println("Unitary evolution timing:")
        TopologicalMarkers._configure_evolution!(false)
        @time @evolution [
            :ham => h => H,
            P0 => h => P
        ] for t in 0:1:τ
            cur = currents(H, P)
            plot_auto("f" => P => cur * 10, clims=(0.98, 1.02))
            frame(a)
        end
        println("Simplified unitary evolution timing:")
        TopologicalMarkers._configure_evolution!(true)
        @time @evolution [
            :ham => h => H,
            P0 => h => P
        ] for t in 0:1:τ
            cur = currents(H, P)
            plot_auto("f" => P => cur * 10, clims=(0.98, 1.02))
            frame(a)
        end
        TopologicalMarkers._configure_evolution!(false)
    end
end
include("operator_test.jl")
include("field_test.jl")
include("evolution_test.jl")
include("lcm_currents_test.jl")
