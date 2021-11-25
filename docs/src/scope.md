# TopologicalMarkers.jl

```@index
Modules = [TopologicalMarkers]
```

## Hamiltonian matrix

```@docs
hamiltonian
field!
domains!
@landau
@symm
@flux
```

## Other linear operators

```@docs
coord_operators
filled_projector
currents
```

## [Bianca-Resta currents](@id bianca_resta_currents)

```@autodocs
Modules = [TopologicalMarkers]
Pages = [joinpath("formulas", "lcm_currents.jl")]
```

## Unitary evolution

```@docs
evolution_operator
@evolution
```

## Data visualization

```@docs
CoordinateRepr
heatmap_data
quiver_data
```

## Plotting shorthands

These functions may come in handy when you want to plot many heatmaps and electric current diagrams in one line.

```@docs
plot_boundaries!
plot_figure!
optimal_layout
plot_auto
```