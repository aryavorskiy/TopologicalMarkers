# Linear operators

## Hamiltonian generation

The [`hamiltonian`](@ref) function can be used to generate a Chern insulator hamiltonian. 
The only essential parameter is `m_lattice` - it is a matrix with the size of the Chern insulator lattice which contains the value of the `m` parameter for each site. Its type must be `CoordinateRepr`.
There is an alternative way to set the `m_lattice` parameter - if its type is a `Matrix`, an additional argument of type `Symbol` is required - a representation specifier
(see [Coordinate Representation](visual.md#Coordinate-representation) for more detail).

Here is an example of both usages:

```@setup ham_test
using TopologicalMarkers
```

```@example ham_test
m_lattice = ones(15, 15)
m_lattice[6:10, 7:11] .= -1
m_repr = CoordinateRepr(m_lattice, :n)
H1 = hamiltonian(m_repr)
H2 = hamiltonian(m_lattice, :n)

H1 == H2 # These are equal ways to specify the m_lattice parameter
```

There is an ability to isolate some sites from other ones. To do that, a `zone_mapping` matrix is required. 
It is a `CoordinateRepr{Symbol}` object, which maps lattice sites to symbols, where different symbols mean different zones, 
and hoppings between sites in different zones are set to zero.
You can define the mapping in the `hamiltonian` function call, or apply it to an existing hamiltonian matrix using the `zones!` function:

```@example ham_test
m_lattice = ones(15, 15)

mapping = CoordinateRepr(fill(:zone1, 15, 15))
mapping[6:10, 6:10] .= :zone2

H1 = hamiltonian(m_lattice, :n, zone_mapping = mapping)
H2 = hamiltonian(m_lattice, :n)
zones!(H2, mapping)

H1 == H2 # These are equal ways to specify the zone mapping
```

You also can apply some magnetic field to the Chern insulator. 
This can be done by passing a function that takes the position vector and returns the gauge vector potential.
You can write your own function  use one of the magnetic field macros defined: 
- `@landau` for Landau gauge field
- `@symm` for symmetric gauge field
- `@flux` for flux quantum
You can define the field in the `hamiltonian` function or apply it to an existing hamiltonian matrix using the `field!` function:

```@example ham_test
m_lattice = ones(15, 15)
B = 0.1

H0 = hamiltonian(m_lattice, :n, field = (x -> [0, x[1] * B, 0]))
H1 = hamiltonian(m_lattice, :n, field = @landau(B))
H2 = hamiltonian(m_lattice, :n)
field!(H2, @landau(B))

H0 == H1 == H2 # These are equal ways to specify the field
```

## Other linear operators

To calculate a density matrix for zero temperature, you can use the [`filled_projector`](@ref) function. 
It accepts the hamiltonian matrix and returns a density matrix. You can also specify the Fermi level (which is zero by default).

You can also generate the coordinate operators using the [`coord_operators`](@ref) function. It takes an optional parameter `symm` which defines if the operators are defined symmetrically (i. e. the point with `(0, 0)` coordinates is in the center of the lattice) or not (i. e. the point with `(0, 0)` coordinates is in the bottom-left angle of the lattice).

Here is an example of usage of both functions to calculate the local Chern marker:

```@example ham_test
using Plots, LinearAlgebra

m_lattice = ones(15, 15)
m_lattice[6:10, 7:11] .= -1
H = hamiltonian(m_lattice, :n)

P = filled_projector(H)
X, Y = coord_operators()
c = 4pi * im * P * X * (I - P) * Y * P

heatmap(heatmap_data(c), color=:viridis)
```

## Electric current

You can calculate electric currents using the [`currents`](@ref) function. 
It accepts the hamiltonian and the density matrix to produce a matrix with currents.
Here is an example of electric currents in a magnetic field:

```@example ham_test
using Plots

B = 0.1
H = hamiltonian(ones(15, 15), :n, field=@symm(B))
P = filled_projector(H)
ps, qs = quiver_data(currents(H, P) * 10 / B)
quiver(ps, quiver = qs)
```

## LCM current

There are several formulas which allow us to calculate currents of the local Chern marker. These formulas comply to the continuity equation, but, unlike electric currents, they are not localized. In this program, there are several macros that can interest you:

```@autodocs
Modules = [TopologicalMarkers]
Pages = [joinpath("formulas", "lcm_currents.jl")]
Filter = f -> nameof(f) != Symbol("@currents")
```

Each one takes matrices of the hamiltonian, the density and coordinate operators and returns a function that calculates a current given indices of 2 sites.

To generate a matrix with currents, use the @currents macro:

```@docs
@currents
```
