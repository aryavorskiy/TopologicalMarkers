# Linear operators

## Hamiltonian generation

!!! info "Ambiguous coordinate representation"
    Trying to map lattice sites to some values (for example, value of the $m$ parameter) can be somewhat ambiguous.

    Usually the first index of the matrix stands for the row number and the second one for the column number.
    And when you write out the matrix values, row number 1 is on the upper side.

    However, if we treat the first index as the $x$-coordinate and the second one as the $y$-coordinate, the placement of the lattice values on the coordinate plane will change.

    To deal with this ambiguity, in all places where sites are mapped to values, use the `CoordinateRepr` struct.
    See [Coordinate Representation](@ref coord_repr) section for more detail.


The [`hamiltonian`](@ref) function can be used to generate a Chern insulator hamiltonian using the following formula:

$\hat{H} = 
\sum_i m_i c^\dagger_i \sigma_z c_i + 
\sum_{x-links} c^\dagger_i \frac{\sigma_z - i \sigma_x}{2} c_j + 
\sum_{y-links} c^\dagger_i \frac{\sigma_z - i \sigma_y}{2} c_j + 
h. c.$

The only essential parameter is `m_lattice` - it is a `CoordinateRepr` that defines the $m$ parameter for each site.
There is an alternative way to set the `m_lattice` - with two arguments, a `Matrix` and a `Symbol` representation specifier.

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

H0 = hamiltonian(m_lattice, :n, field = (r -> [0, r[1] * B, 0]))
H1 = hamiltonian(m_lattice, :n, field = @landau(B))
H2 = hamiltonian(m_lattice, :n)
field!(H2, @landau(B))

H0 == H1 == H2 # These are equal ways to specify the field
```

## Density and coordinate operators

To calculate a density matrix for zero temperature, you can use the [`filled_projector`](@ref) function. 
It accepts the hamiltonian matrix and returns a density matrix. You can also specify the Fermi level (which is zero by default).

You can also generate the coordinate operators using the [`coord_operators`](@ref) function. 
It takes an optional parameter `symm` which defines if the operators are defined symmetrically 
(i. e. the point with `(0, 0)` coordinates is in the center of the lattice) or not 
(i. e. the point with `(0, 0)` coordinates is in the bottom-left angle of the lattice).

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

## Bianca-Resta currents

There are several ways to calculate currents of the local Chern marker. 
For example, if we define it using the Bianca-Resta formula $c(r) = 4\pi i \langle r | PXQYP | r \rangle$, 
we can define the current $J(r, r')$ as some formula that complies the following rules:

$$
\begin{cases}
    d_t c(r) = \sum_{r'} J(r, r') \\
    J(r, r') = - J(r', r)
\end{cases}, \hspace{0.5cm}
d_t c(r) = 4\pi i \langle r | i[H, PXQYP] | r \rangle
= 4\pi i \langle r | i([H, P]XQYP + PX[H, Q]YP + PXQY[H, P]) | r \rangle
$$

Let us call them **Bianca-Resta currents**.

There are four formulas like that in this package. Each is a macro that takes matrices of the hamiltonian, the density and coordinate operators
to generate a function that calculates the current between two sites (given their indices as parameters).

The difference between them is that some of them are translationally invariant - if we shift the coordinate operators (e. g. redefine $\hat{X}$ as $\hat{X} + X_0$), the values of currents will not change. Also some of them are _stable_ - all currents are equal to zero if $P$ is a density matrix of a stationary state of hamiltonian $\hat{H}$ that does not depend on time.

These traits of macros in this package are in the following table:

|Macro|Translational invariance|Stability|
|---|---|---|
|`J_b`| No | No |
|`J_tr`| Yes | No |
|`J_eq`| No | Yes |
|`J_treq`| Yes | Yes |

To find out more about Bianca-Resta currents, [check out this section](scope.md#LCM-Currents).

Each one takes  and returns a function that calculates a current given indices of 2 sites.

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

    `TopologicalMarkers._set_lattice_size!(needed_lattice_size)`

    This function accepts two integers or a tuple of two integers as arguments.

    Note that if no hamiltonian was generated since module import, you will definitely get this error if you try to evaluate some other operator.
    Do not do it. Please.