print("Importing test dependencies... ")

using Test
using TopologicalMarkers
using LinearAlgebra

println("done.")

# Macro tests
@testset "Static LCM" begin
    siz = (15, 15)
    global B = 1e-8
    ms = CoordinateRepr(ones(siz) * -1)
    global ham = hamiltonian(ms)
    global ham_b = hamiltonian(ms, field=@symm(B))
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
end

include("operator_test.jl")
include("field_test.jl")
include("lcm_currents_test.jl")