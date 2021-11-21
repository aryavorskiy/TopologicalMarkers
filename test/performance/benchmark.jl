using BenchmarkTools

print("-"^25, " BENCHMARKS ", "-"^25)

B = 1e-6

print("\nHamiltonian + field apply benchmark: ")
b = @benchmark hamiltonian(CoordinateRepr(ones(11, 12)), field=@symm(B))
show(stdout, "text/plain", b)

ham = hamiltonian(CoordinateRepr(ones(11, 12)))
ham_b = hamiltonian(CoordinateRepr(ones(11, 12)), field=@symm(B))
P = filled_projector(ham)
P_b = filled_projector(ham_b)
streda = (P_b - P) / B
curr_b = currents(ham_b, P_b) / B
X, Y = coord_operators()

print("\nAuto plot timing: ")
b = @benchmark plot_auto(streda => curr_b,
    "LCM" => 4π * im * P * X * P * Y * P,
    cutaway_view=(8, 8)) seconds=15
show(stdout, "text/plain", b)