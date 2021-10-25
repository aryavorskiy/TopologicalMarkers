# TopologicalMarkers.jl

A package that simplifies calculation of different topological markers.

## Installation

To install it, simply copy this line to julia's REPL and execute it:

```
]add https://github.com/aryavorskiy/TopologicalMarkers/
```

## Examples

Let us take a Chern insulator, set the $m$ parameter to $-1$ in the middle of the lattice and $1$ everywhere else,
and then evaluate the local Chern marker using both traditional and Streda formulas (i. e. the linear response of the local density to the magnetic field):

```@example
using TopologicalMarkers
using Plots

m_lattice = ones(25, 25)
m_lattice[11:15, 11:15] .= -1
H = hamiltonian(m_lattice, :c)
P = filled_projector(H)
X, Y = coord_operators()
ch = -4π * im * P * X * P * Y * P

B = 0.01
Hb = hamiltonian(m_lattice, :c, field=@landau(B))
Pb = filled_projector(Hb)
str = (Pb - P) / B

plot_auto("LCM" => ch, "Streda" => str, 
    hmapclims=(-1.5, 1.5), currentscolor=:yellow, control_site=(13, 13), markercolor=:brown)
```

See (Hamiltonian generation)[@ref] and (Visualization)[@ref] for detailed explanation.

Another good example is a problem where unitary evolution is used. 
Here we create an animation of the local density changing during adiabatic magnetic field-on:

```@example
using TopologicalMarkers
using Plots

ms = CoordinateRepr(ones(15, 15))
Bf = 0.01
τ = 30
time_domain = 0:0.5:τ

H0 = hamiltonian(ms)
P0 = filled_projector(H0)
h(t) = hamiltonian(ms, field=@symm(Bf * t / τ))
a = Animation()
@evolution [
    :ham => h => H,
    P0 => h => P
] for t in time_domain
    cur = currents(H, P)
    plot_auto("Local density" => P => cur * 100, hmapclims=(0.9, 1.1))
    frame(a)
end

gif(a, "example_animation.gif", fps=10)
```

See (Unitary evolution)[@ref] for detailed explanation.