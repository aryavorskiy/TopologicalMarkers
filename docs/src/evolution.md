# Unitary evolution

Suppose we want to study the behavior of some quantum system in time-dependent conditions. We can use the unitary evolution operator to describe how the density matrix depends on time:

$\mathcal{U}(t) = T\left\{ e^{\frac{1}{i\hbar} \int_{t_0}^t \hat{H}(\tau) d\tau} \right\},\hspace{0.5cm}
\mathcal{P}(t) = \mathcal{U}(t) \mathcal{P}_0 \mathcal{U}^\dagger (t)$

## The static evolution function

The [`evolution_operator`](@ref) function calculates the evolution operator for a time-independent hamiltonian. 
It also stores the output data, so if you call it multiple times with the same input, the output will be calculated only once.

## The evolution macro

This macro can be quite useful if your hamiltonian depends on time or if there are multiple hamiltonians in your experiment.
What you need is a function that takes the time and returns the matrix of the hamiltonian at that moment of time:

```julia
h(t) = hamiltonian(ms, field=@symm(Bf * t / τ))
```

Here `h(t)` describes the magnetic field being adiabatically turned on.

Let's take a look at the example 2 from the [Home](index.md) section.
The `@evolution` macro creates a block where the hamiltonian and density matrices are evaluated for the given time interval. 
It takes two arguments - a list/vector with evolution specifiers and a for-loop that iterates over the time interval:

```julia
a = Animation()

@evolution [
    :ham => h => H,
    P0 => h => P
] for time in time_domain
    cur = currents(H, P)
    plot_auto("Local density" => P => cur * 40, hmapclims=(0.98, 1.02))
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