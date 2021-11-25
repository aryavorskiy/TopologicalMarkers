using BenchmarkTools
using StatProfilerHTML

print("-"^25, " BENCHMARKS ", "-"^25)

B = 1e-6

print("\nHamiltonian + field apply benchmark: ")
b = @benchmark hamiltonian(CoordinateRepr(ones(11, 12)), field=@symm(B))
display(b)

ham = hamiltonian(CoordinateRepr(ones(11, 12)))
ham_b = hamiltonian(CoordinateRepr(ones(11, 12)), field=@symm(B))
P = filled_projector(ham)
P_b = filled_projector(ham_b)
streda = (P_b - P) / B
curr_b = currents(ham_b, P_b) / B
X, Y = coord_operators()
chern = 4π * im * P * X * P * Y * P

print("\nHeatmap data timing: ")
b = @benchmark heatmap_data(P)
display(b)

print("\nQuiver data timing: ")
b = @benchmark currents(ham_b, P_b)
display(b)

print("\nAuto plot timing: ")
b = @benchmark plot_auto(streda => curr_b,
    "LCM" => chern,
    cutaway_view=(8, 8)) seconds=15
display(b)
