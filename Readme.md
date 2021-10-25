# TopologicalMarkers

This package contains software used in the article <РЫБА>. 

## Installation

This package can be installed using Julia's package manager via REPL:

```julia-repl
(@v1.6) pkg> add https://github.com/aryavorskiy/TopologicalMarkers
```

## Usage

### Hamiltonian matrix generation

```julia
using TopologicalMarkers

# connection parameter setup
sz = (15, 15)
m_lattice = ones(sz)
m_lattice[6:10, :] = 3

# zone mapping setup
zones = fill(:A, sz)
zones[6:10, :] = :B

# generate hamiltonian matrix
H = hamiltonian(m_lattice, :c)
zones!(H, zones, :c)

# magnetic field
B = 0.1
Hb = copy(H)
field!(H, @landau(B))
```

### Unitary evolution

```julia
using TopologicalMarkers, Plots

sz = (15, 15)
H1 = hamiltonian(ones(sz), :c)
H2 = hamiltonian(ones(sz) * 3, :c)

P_0 = filled_projector(H1) # generate ground state density matrix
X, Y = coord_operators()   # generate coordinate operators
h(t) = H2 # set hamiltonian time dependence

plot(title="LCM quench")

@evolution [
    :ham => h => H,
    P_0 => h => P
] for t in 0:10 
    ch = -4pi * im * P * X * P * Y * P
    lcm = heatmap_data(ch)
    plot!(lcm[8, :])
end
```

### Visualization tools

```julia
using TopologicalMarkers, Plots

sz = (15, 15)
H1 = hamiltonian(ones(sz), :c)
H2 = hamiltonian(ones(sz) * 3, :c)

P_0 = filled_projector(H1) # generate ground state density matrix
X, Y = coord_operators()   # generate coordinate operators
h(t) = H2 # set hamiltonian time dependence

plot(title="LCM quench")
anim = Animation()

@evolution [
    :ham => h => H,
    P_0 => h => P
] for t in 0:10 
    ch = -4pi * im * P * X * P * Y * P
    plot_auto("LCM" => ch, control_site = (8, 8))
    frame(anim)
end

gif(anim, "animation.gif", fps = 4)
```

See docs for more details.