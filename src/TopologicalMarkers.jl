module TopologicalMarkers

# TODO finish bugfixes and tests, write docs
# TODO profile and optimize plot_auto

export hamiltonian, field!, zones!, # hamiltonian.jl
    filled_projector, coord_operators, currents, # operators.jl
    @evolution, evolution_operator, # evolution.jl
    CoordinateRepr, pair_to_index, adjacent_sites, # utils.jl
    heatmap_data, quiver_data, quiver_currents!, plot_boundaries!, optimal_layout, plot_marker!, plot_auto, # visual.jl
    @J_c, @J_m, @J, @J_m_inv, @J_eq, @J_inv, @J_best, @currents, # lcm_currents.jl
    @landau, @symm, # field.jl
    local_chern, streda_p # markers.jl

include("utils.jl")
include("visual.jl")
include("operators.jl")
include("evolution.jl")
include("plot.jl")
include("hamiltonian.jl")
include("field.jl")
include(joinpath("formulas", "lcm_currents.jl"))
include(joinpath("formulas", "markers.jl"))
end # module
