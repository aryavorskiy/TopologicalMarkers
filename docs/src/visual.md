# Visualization

## [Coordinate representation](@id coord_repr)

The [`CoordinateRepr`](@ref) struct is a workaround to deal with certain ambiguity of matrix graphical representations. 
We usually index its elements row-first, and print the rows in ascending order (i. e. the natural representation). 
However, if you plot a heatmap for a matrix, its first index is treated as $x$, the second - as $y$, 
and the $Ox$ axis is directed upward. That's where we need this struct.

The constructor takes two arguments - the matrix itself and a symbol - a _representation specifier_, 
which tells if the matrix should be treated as _natural_ (`:n` or `:natural`) or _coordinate_ (`:c` or `:coord`) representation.

Let us look at an example:

```@example repr_comp
mat = zeros(12, 12)
mat[1:5, 1:5] .= 1
mat[1:10, 6:10] .= -1
mat
```

You can see how this matrix looks when printed in a _natural_ way. Let us look at a heatmap of both representations:

```@example repr_comp
using TopologicalMarkers, Plots

p = plot(layout = (1, 2), size = (800, 300))
heatmap!(p[1], title="'Natural' representation", CoordinateRepr(mat, :n))
heatmap!(p[2], title="'Coordinate' representation", CoordinateRepr(mat, :c))
```

## Data processing

Most 'raw' data like operator or current matrices require additional processing before being plotted.

The [`heatmap_data`](@ref) function generates a `CoordinateRepr` object, given a linear operator matrix. It evaluates traces of diagonal elements of the operator in the coordinate representation.

```julia
# c is the local Chern operator
c = 4pi * im * P * X * (I - P) * Y * P

marker = heatmap_data(c)
heatmap(marker, color=:viridis)
```

The [`quiver_data`](@ref) function generates data to be added to a quiver plot. The input is a matrix with currents, the output is a tuple of vectors of the same length - one contains the origin points of the arrows, the other contains the arrow vectors.

```julia
# j is the electric currents matrix
j = currents(H, P)

ps, qs = quiver_data(j)
quiver(ps, quiver = qs, color=:blue)
```

This function supports additional arguments - suppose you need to plot only arrows that are located in some specific domain, 
currents between sites that are to far away from each other (e. g. more than 5) do not need to be displayed 
and arrows shorter than `0.1` must be omitted. It is possible then to specify these parameters as follows:

```julia
# j is the electric currents matrix
j = currents(H, P)
lims = (5, 10)

ps, qs = quiver_data(j; threshold = 0.1, dist_threshold = 5, xlims = lims, ylims = lims)
quiver(ps, quiver = qs, color=:blue)
```

## Automatic plotting

These functions may come in handy if you want to plot multiple figures simultaneously and do not want to set the layout up manually.

The [`plot_figure!`](@ref) function can be used to plot one figure. It has one argument - the `Plot` object to draw on. 
You can draw a heatmap, a quiver plot or bounding lines between domains if you need to. What you have to do is pass a keyword argument:
- `hmap` for heatmap
- `currents` for quiver
- `domain_mapping` for boundaries

The `xlims` and `ylims` keywords have meaning traditional for plots.
All other keyword arguments are automatically parsed and then passed to the plot function at different stages depending on the prefix of the key.
- `hmap` is passed to heatmap
- `currents` is passed to quiver
- `bounds` is passed to boundaries
For example, `hmapcolor = :viridis` means that the argument `color = :viridis` will be passed to the function when plotting the heatmap. 

Here is an example:

```@setup vis_test
using TopologicalMarkers, Plots
```

```@example vis_test
ms = ones(15, 15)
ms[6:10, 6:10] .= 3
B = 0.1

zs = fill(:ext, 15, 15)
zs[6:10, 6:10] .= :isolate
domains = CoordinateRepr(zs)

H = hamiltonian(ms, :c, field = @landau(B), domain_mapping = domains)
P = filled_projector(H)
cur = currents(H, P)
p = plot(title="State density in two-domain insulator")
plot_figure!(p, hmap = P, currents = cur / maximum(abs.(cur)), bounds = domains, 
    hmapcolor = :viridis, currentscolor = :blue, boundsstyle = :dashdotdot, 
    boundscolor = :red)
```

The [`plot_auto`](@ref) function generates an optimal layout for the given amount of figures and then invokes `plot_figure!` multiple times.

The arguments it takes are pairs like this:
```julia
title => heatmap => currents
```

Here the first and the third arguments can be omitted.

domain mapping can be set via the same keyword argument. Also you can create figures with cutaway views if you define cutaway viewpoints via `cutaway_views` or `cutaway_view` arguments.

Each cutaway viewpoint is a lattice site. For each one the function takes heatmap values for every site with the same y-coordinate and plots the heatmap values dependent on `x` on a separate plot. Each cutaway view is shown on the plot as a line, which will be later referred to as a splitline.

This function also supports keyword argument parsing. The following prefixes are supported:
- `hmap` is passed to heatmap
- `currents` is passed to quiver
- `bounds` is passed to boundaries
- `splitline` is passed to splitline
- `cutaway` is passed to cutaway views

A good example is the piece of code on the main page of this website. Let us modify it a bit - add another cutaway viewpoint, set the splitline color to brown, and cutaway view line style to `:dashdot`

```@example vis_test
m_lattice = ones(25, 25)
m_lattice[11:15, 11:15] .= -1
H = hamiltonian(m_lattice, :c)
P = filled_projector(H)
X, Y = coord_operators()
ch = -4π * im * P * X * P * Y * P

B = 0.01
Hb = hamiltonian(m_lattice, :c, field = @landau(B))
Pb = filled_projector(Hb)
str = (Pb - P) / B

plot_auto("LCM" => ch, "Streda" => str, domain_mapping = CoordinateRepr(m_lattice), 
    plot_size = (800, 600), hmapclims = (-1.5, 1.5), currentscolor = :yellow, 
    cutaway_view = (13, 13), splitlinecolor = :brown, cutawaystyle = :dashdot)
```