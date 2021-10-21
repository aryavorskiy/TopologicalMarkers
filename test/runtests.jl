print("Importing test dependencies... ")

using Test
using TopologicalMarkers
using LinearAlgebra
using Plots

println("done.")

# Macro tests
@testset "Usability tests" begin
    @testset "Static LCM" begin
        siz = (15, 15)
        global B = 1e-8
        ms = ones(siz) * -1
        global ham = hamiltonian(ms, :c)
        global ham_b = hamiltonian(CoordinateRepr(ms), field=@symm(B))
        # field!(ham_b, @symm(B))

        global P = filled_projector(ham)
        P_b = filled_projector(ham_b)
        curr_b = currents(ham_b, P_b) / B
        streda = (P - P_b) / B
        plot_auto("Streda" => heatmap_data(streda, siz) => curr_b)

        global X, Y = coord_operators()
        plot_auto("Streda" => streda => curr_b,
        "LCM" => 4π * im * P * X * P * Y * P,
        control_site=(8, 8), lattice_size=siz)
        print("Auto plot timing: ")
        @time plot_auto(streda => curr_b,
        "LCM" => 4π * im * P * X * P * Y * P,
        control_site=(8, 8), lattice_size=siz)

        pl = plot()
        currs = @currents @J_best ham_b P X Y
        plot_marker!(pl, hmap=4π * im * P * X * P * Y * P, currents=currs, xlims=(6,10), ylims=(8, 12))
    end
end
include("operator_test.jl")
include("field_test.jl")
include("lcm_currents_test.jl")