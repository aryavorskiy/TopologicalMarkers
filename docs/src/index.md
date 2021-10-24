# TopologicalMarkers.jl

A package that simplifies calculation of different topological markers.

To install it, simply copy this line to julia's REPL and execute it:

```
]add https://github.com/aryavorskiy/TopologicalMarkers/
```

## Examples

Let us take a Chern insulator, set the $m$ parameter to $-1$ in the middle of the lattice and $1$ everywhere else,
and then evaluate the local Chern marker using both traditional and Streda formulas:

```@example
using TopologicalMarkers
using Plots

m_lattice = ones(25, 25)
m_lattice[6:10, 6:10] .= -1
H = hamiltonian(m_lattice, :c)
P = filled_projector(H)
X, Y = coord_operators()
ch = -4π * im * P * X * P * Y * P

B = 0.01
Hb = hamiltonian(m_lattice, :c, field=@landau(B))
Pb = filled_projector(Hb)
str = (Pb - P) / B

plot_auto("LCM" => ch, "Streda" => str, 
    hmapclims=(-1.5, 1.5), currentscolor=:yellow, control_site=(8, 8), markercolor=:brown)
savefig("example.png"); nothing # hide
```

The code here will produce the following graph:

![](example.png)

_TODO add more examples_