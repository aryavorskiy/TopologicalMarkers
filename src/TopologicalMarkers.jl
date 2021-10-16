module TopologicalMarkers

# TODO finish bugfixes, write docs, tests
# TODO new features: plot_arranged keyword arguments for plots, a file with marker calculators

export hamiltonian, field!, zones!, # hamiltonian.jl
    filled_projector, coord_operators, currents, # operators.jl
    @evolution, evolution_operator, # evolution.jl
    CoordinateRepr, # utils.jl
    heatmap_data, quiver_data, quiver_currents!, plot_boundaries!, plot_arranged, # visual.jl
    @J_c, @J_m, @J, @J_eq, @J_inv, @J_best, @currents, # lcm_currents.jl
    @landau, @symm # field.jl

include("utils.jl")

include("hamiltonian.jl")
include("operators.jl")
include("evolution.jl")
include("visual.jl")
include("lcm_currents.jl")
include("field.jl")
end # module
