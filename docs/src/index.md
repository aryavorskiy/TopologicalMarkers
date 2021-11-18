# TopologicalMarkers.jl

A package that simplifies calculation of different topological markers.

## Installation

To install it, simply copy this line to julia's REPL and execute it:

```
]add https://github.com/aryavorskiy/TopologicalMarkers/
```

## Examples

### Different local markers in an equilibrium

Let us take a Chern insulator, set the $m$ parameter to $-1$ in the middle of the lattice and $1$ everywhere else,
and then evaluate the local Chern marker using both Bianca-Resta and Streda formulas (i. e. the linear response of the local density to the magnetic field):

```@example
using TopologicalMarkers
using Plots

# Set the m parameter to -1 on the whole lattice
# Let the phase be topological
m_lattice = ones(25, 25)
m_lattice[11:15, 11:15] .= -1
H = hamiltonian(m_lattice, :c)

P = filled_projector(H)             # Ground state density matrix
X, Y = coord_operators()            # Coordinate operator matrices
ch = -4π * im * P * X * P * Y * P   # Bianca-Resta local Chern operator

# Calculate the perturbation of the density matrix
B = 0.01
Hb = hamiltonian(m_lattice, :c, field = @landau(B))
Pb = filled_projector(Hb)
str = (Pb - P) / B

# Finally, plot everything
plot_auto("Bianca-Resta" => ch, "Streda" => str, plot_size = (800, 600),
    hmapclims = (-1.5, 1.5), currentscolor = :yellow, cutaway_view = (13, 13), 
    markercolor = :brown)
```

See [Hamiltonian generation](@ref) and [Visualization](@ref) for detailed explanation.

### Adiabatic flux quantum application

To understand better how the particle density and electric currents react to adiabatic flux quantum appearance, let us create an animation of it.
The Chern insulator is in topological phase, with $m = 1$.

```@example
using TopologicalMarkers
using Plots

# Set the m parameter to 1 on the whole lattice
# Let the phase be topological again
ms = CoordinateRepr(ones(15, 15))
Bf = 0.01
τ = 30
time_domain = 0:0.5:τ

# Define how the hamiltonian depends on time
H0 = hamiltonian(ms)
P0 = filled_projector(H0)
h(t) = hamiltonian(ms, field = @flux(Bf * t / τ))   # This is a function definition

a = Animation()

# The code in the square brackets defines `H` and `P` variables in the loop body
# The comma can be omitted
@evolution [
    :ham => h => H,
    P0 => h => P
] for t in time_domain
    cur = currents(H, P)
    plot_auto("Local density" => P => 100cur, plot_size = (800, 600), hmapclims = (0.9, 1.1))
    frame(a)
end

# Export the animation
gif(a, "example_animation.gif", fps = 10)
```

See [Unitary evolution](@ref unitary_evolution) for detailed explanation.

!!! tip
    There are more examples of different use cases in this notebook: [github.com/aryavorskiy/TopologicalMarkers.jl/supplementary.ipynb](https://www.youtube.com/watch?v=E4WlUXrJgy4)
    All pictures from the article were created using this notebook.