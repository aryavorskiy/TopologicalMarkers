using Profile
using StatProfilerHTML

# Init lattice
expname = :adiabatic_move_check

τ = 100
time_domain = 0:10:2τ

TopologicalMarkers._set_lattice_size!(22, 21)
X, Y = coord_operators()

m0 = fill(-3, 22, 21)
m0[6:12, 7:13] .= -1

m1 = fill(-3, 22, 21)
m1[7:13, 7:13] .= -1

h(t) = hamiltonian(m0 + (m1 - m0) * min(1, t/τ), :c)

B = 1e-6

function hb(t)
    local H = h(t)
    field!(H, @symm(B))
    return H
end
P0 = filled_projector(h(0))
P0b = filled_projector(hb(0))

TopologicalMarkers._configure_evolution!(true)

@profilehtml @time @evolution [
    :ham => h => H,
    :ham => hb => Hb,
    P0 => h => P,
    P0b => hb => Pb,
] for t in time_domain
    cur = currents(H, P)
    curb = currents(Hb, Pb)
end
