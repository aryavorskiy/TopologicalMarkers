# TopologicalMarkers

A package that simplifies calculation of different topological markers.

'''julia
m_lattice = ones(25, 25)
m_lattice[6:10, 6:10] .= -1
B = 0.
H = hamiltonian(m_lattice, :c)
P = filled_projector(H)
c = currents(H, P)
X, Y = coord_operators()
ch = -4π * im * P * X * P * Y * P

plot_auto("LCM" => ch => c * 10, hmapclims=(-1.5, 1.5), currentscolor=:yellow, control_site=(8, 8))
'''
