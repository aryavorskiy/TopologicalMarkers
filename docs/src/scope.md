# TopologicalMarkers.jl

```@index
Modules = [TopologicalMarkers]
Filter = f -> (n = String(nameof(f)); !(startswith(n, "@J")) && n != "@currents")
```

## Hamiltonian matrix

```@docs
hamiltonian
field!
zones!
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

## LCM currents

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

These functions may come in handy when you want to plot many marker heatmaps and electric current diagrams in one line.

```@docs
plot_boundaries!
plot_marker!
plot_auto
```