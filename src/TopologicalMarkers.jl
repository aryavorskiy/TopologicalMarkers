module TopologicalMarkers

# TODO finish bugfixes and tests, write docs
# TODO profile and optimize plot_auto
# TODO new features for plots: allow title omitting
# TODO new attrs: axis labels,

export hamiltonian, field!, zones!, # hamiltonian.jl
    filled_projector, coord_operators, currents, # operators.jl
    @evolution, evolution_operator, # evolution.jl
    pair_to_index, adjacent_sites, # utils.jl
    CoordinateRepr, heatmap_data, quiver_data, # visual.jl
    plot_boundaries!, optimal_layout, plot_marker!, plot_auto, # plot.jl
    @J_c, @J_m, @J, @J_m_inv, @J_eq, @J_inv, @J_best, @currents, # lcm_currents.jl
    @landau, @symm, @flux, # field.jl
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
