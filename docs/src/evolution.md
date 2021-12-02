# [Unitary evolution](@id unitary_evolution)

Suppose we want to study the behavior of some quantum system in time-dependent conditions. We can use the unitary evolution operator to describe how the density matrix depends on time:

$\mathcal{U}(t) = T\left\{ e^{\frac{1}{i\hbar} \int_{t_0}^t \hat{H}(\tau) d\tau} \right\},\hspace{0.5cm}
\mathcal{P}(t) = \mathcal{U}(t) \mathcal{P}_0 \mathcal{U}^\dagger (t)$

## The static evolution function

The [`evolution_operator`](@ref) function calculates the evolution operator for a time-independent hamiltonian. 
Its input parameters are the hamiltonian matrix `H` and the time interval `t`.

## [Matrix exponent optimization](@id mexp_opt)

The matrix exponent is the heaviest linear algebra operation used in this project. 
To speed up calculations in some times, you can set up the matrix exponent to be calculated in a simpler way with this function:

```julia
TopologicalMarkers._configure_evolution!(simplify::Bool; <keyword arguments...>)
```

If the `simplify` parameter is set to to `true`, the matrix exponent is evaluated using the Taylor series.
You can set other parameters with following keywords:

- `order`: The quantity of members of the Taylor expansion to be calculated
- `threshold`: To avoid precision loss, if the `t` parameter is greater than `threshold`, the exponent will be evaluated using the `exp` function. Set to `nothing` to always use Taylor expansion.

## The evolution macro

This macro can be quite useful if your hamiltonian depends on time or if there are multiple hamiltonians in your experiment.
Let us define a function that takes the time and returns the hamiltonian:

```julia
h(t) = hamiltonian(ms, field = @symm(Bf * min(t, τ) / τ))
```

Here `h(t)` describes the magnetic field being adiabatically turned on.

Take a look at [this example](@ref example_flux) from the Home/Examples section.
The `@evolution` macro creates a block where the hamiltonian and density matrices are evaluated for the given time interval. 
It takes two arguments - a list/vector with evolution specifiers and a for-loop that iterates over the time interval:

```julia
a = Animation()

@evolution [
    :ham => h => H,
    P0 => h => P
] for time in time_domain
    cur = currents(H, P)
    plot_auto("Local density" => P => cur * 40, hmapclims = (0.98, 1.02))
    frame(a)
end
```

Let us make it clear what an evolution specifier is: it is a pair chain with three arguments:

```julia
arg1 => arg2 => arg3
```

There are two options what these arguments can be:

- The first argument is a density matrix at the moment $t = 0$. The second is the function describing the hamiltonian time dependence, and the third one is the name for the density matrix at the $t = \text{time}$ moment - a variable with such a name will be defined inside the loop body.
- The first argument is the `:ham` symbol. Then the second one is still the hamiltonian function, and the third one is the name for the hamiltonian matrix at the $t = \text{time}$ moment.

So, in the example before in the for-loop body `H` stands for `h(time)`, and `P` is the evolved `P0` density matrix.

!!! note
    It is important that `arg2` must be a _function name_ - therefore you cannot use a lambda expression in such context. 
    This is done with purpose to make the resulting code more readable.

    Note that if the hamiltonian matrix does not depend on time, you still need to define a function. 
    The reason is that it is nearly impossible to detect if it is a function name or a variable name at compile-time.

## Performance

To speed up calculations, one might use libraries such as [CUDA.jl](https://juliagpu.gitlab.io/CUDA.jl/) that provide an alternative linear algebra interface. 
To make them compatible with the `@evolution` macro, you should define these operatots and functions for the new matrix type:

- Equality operator: `==`
- Basic arithmetic functions: `+`, `-`, `*`, `adjoint`
- The matrix exponent `exp(A)`
    - If it is not possible to implement this function (like in `CUDA.jl`), you can define the `one(A)` function which returns the identity matrix with the size same as `A`. 
    Then you need to [configure the evolution operator function to enable the Taylor expansion formula](@ref mexp_opt)

!!! tip
    The `hamiltonian` function is slow enough - this is the price we have to pay for its flexibility.
    The most time-consuming tasks are allocating a new matrix and calling lambda-functions to evaluate the Peierls substitution.

    Avoid calling it when evaluating the `h(t)` for the evolution operator when possible.