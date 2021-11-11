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

There is an ability to isolate some sites from other ones. To do that, a `domain_mapping` matrix is required. 
It is a `CoordinateRepr{Symbol}` object, which maps lattice sites to symbols, where different symbols mean different domains, 
and hoppings between sites in different domains are set to zero.
You can define the mapping in the `hamiltonian` function call, or apply it to an existing hamiltonian matrix using the `domains!` function:

```@example ham_test
m_lattice = ones(15, 15)

mapping = CoordinateRepr(fill(:domain1, 15, 15))
mapping[6:10, 6:10] .= :domain2

H1 = hamiltonian(m_lattice, :n, domain_mapping = mapping)
H2 = hamiltonian(m_lattice, :n)
domains!(H2, mapping)

H1 == H2 # These are equal ways to specify the domain mapping
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

## Density and coordinate operators

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

heatmap(heatmap_data(c), color = :viridis)
```

## Electric current

You can calculate electric currents using the [`currents`](@ref) function. 
It accepts the hamiltonian and the density matrix to produce a matrix with currents.
Here is an example of electric currents in a magnetic field:

```@example ham_test
using Plots

B = 0.1
H = hamiltonian(ones(15, 15), :n, field = @symm(B))
P = filled_projector(H)
ps, qs = quiver_data(currents(H, P) * 10 / B)
quiver(ps, quiver = qs)
```

## LCM current

There are several formulas which allow us to calculate currents of the local Chern marker. These formulas comply to the continuity equation, but, unlike electric currents, they are not localized. In this program, there are several macros that can interest you, all can be found [in the corresponding section](scope.md#LCM-Currents).

Each one takes matrices of the hamiltonian, the density and coordinate operators and returns a function that calculates a current given indices of 2 sites.

To generate a matrix with currents, use the [`@currents`](@ref) macro.

Here is an example animation of LCM behavior under after hamiltonian quench

```@example ham_test
using Plots

H1 = hamiltonian(ones(15, 15), :c)
H2 = hamiltonian(ones(15, 15) * 3, :c)

h(t) = H2
P0 = filled_projector(H1)
X, Y = coord_operators()
a = Animation()
@evolution [
    :ham => h => H,
    P0 => h => P
] for t in 0:0.1:3
    c = -4pi * im * P * X * P * Y * P
    cur = @currents @J_treq H P X Y
    plot_auto("LCM currents" => c => cur / 5, hmapclims = (-4, 1))
    frame(a)
end

gif(a, "example_animation.gif", fps = 10)
```

!!! warning
    You may get an error entitled as `Please specify lattice size explicitly` while running one of these functions. 
    This is caused by the fact that all these function require an additional argument, the lattice size, 
    but by default they use the lattice size for the last hamiltonian generated.

    It is quite unlikely that you encounter this error, but if you did, you can set the lattice size manually by using this function:

    `TopologicalMarkers._set_lattice_size!(lattice_size)`

    This function accepts two integers or a tuple of two integers as arguments.

    Note that if no hamiltonian was generated since module import, you will definitely get this error if you try to evaluate some other operator.
    Do not do it. Please.