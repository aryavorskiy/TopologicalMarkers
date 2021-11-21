using Logging
import LinearAlgebra: eigen

function _set_hopping!(H::Matrix{ComplexF64}, i::Int, j::Int, hop::Matrix)
    H[2 * i - 1:2 * i, 2 * j - 1:2 * j] = hop
    H[2 * j - 1:2 * j, 2 * i - 1:2 * i] = hop'
end

function _hamiltonian(m_repr::CoordinateRepr, pbc::Tuple{Bool, Bool})::Matrix{ComplexF64}
    lattice_size = size(m_repr)
    _set_lattice_size!(lattice_size)
    H = zeros(ComplexF64, (prod(lattice_size) * 2, prod(lattice_size) * 2))

    # Generate diagonal elements
    for i in 1:prod(lattice_size)
        H[2 * i - 1:2 * i, 2 * i - 1: 2 * i] = σ[3] * m_repr[i]
    end

    # Generate hoppings
    for i in 1:lattice_size[1], j in 1:lattice_size[2]
        if i != lattice_size[1] || pbc[1]
            i1 = pair_to_index(lattice_size, i, j)
            i2 = pair_to_index(lattice_size, i+1, j)
            _set_hopping!(H, i1, i2, (σ[3] - im * σ[1]) / 2)
        end
        if j != lattice_size[2] || pbc[2]
            i1 = pair_to_index(lattice_size, i, j)
            i2 = pair_to_index(lattice_size, i, j+1)
            _set_hopping!(H, i1, i2, (σ[3] - im * σ[2]) / 2)
        end
    end
    return H
end

@doc raw"""
    field!(H, A[, lattice_size]; intervals = 10)

Applies magnetic field to specified hamiltonian. In other words, all hoppings get multiplied on a specific phase factor that can be calculated using Peierls substitution:

$\varphi_{ij} = 2\pi \int_i^j A(r) \cdot dr$

This integral is calculated explicitly for every hopping, using the `A` function.

# Arguments
- `H`: the hamiltonian matrix
- `A`: a function that takes a `Vector` representing a point and returns a `Vector` representing the vector potential in that point
- `lattice_size`: the size of the lattice the hamiltonian is defined for. If not provided, this function will use the value for the hamiltonian matrix that was created last
- `intervals`: the number of intervals to use when calculating the Peierls substitution phase factor
"""
function field!(H::Matrix{<:Complex}, A::Function, lattice_size::SizeType = nothing; intervals::Int=10)
    lattice_size = _try_get_lattice_size(lattice_size)     
    local function peierls(i, j)::ComplexF64
        phase::Float64 = 0
        r1 = index_to_pair(lattice_size, i)
        r2 = index_to_pair(lattice_size, j)
        dr = (r2 - r1) / intervals
        for i in 1:intervals
            phase += dr' * A(r1 + (i - 0.5) * dr)[1:2]
        end
        return -2π * im * phase |> exp
    end
    for i in 1:lattice_size[1], j in 1:lattice_size[2]
        i1 = pair_to_index(lattice_size, i, j)
        i2 = pair_to_index(lattice_size, i+1, j)
        phase_mod = peierls(i1, i2)
        H[2 * i1 - 1:2 * i1, 2 * i2 - 1:2 * i2] *= phase_mod
        H[2 * i2 - 1:2 * i2, 2 * i1 - 1:2 * i1] *= phase_mod'
        
        i2 = pair_to_index(lattice_size, i, j+1)
        phase_mod = peierls(i1, i2)
        H[2 * i1 - 1:2 * i1, 2 * i2 - 1:2 * i2] *= phase_mod
        H[2 * i2 - 1:2 * i2, 2 * i1 - 1:2 * i1] *= phase_mod'
    end
end

"""
    domains!(H, domain_mapping[, repr_spec])

Divides the Chern insulator hamiltonian into several unconnected domains. The hoppings between these domains are erased.

# Arguments
- `H`: the hamiltonian matrix
- `domain_mapping`: an `AbstractMatrix{Symbol}` or `CoordinateRepr{Symbol}`. Each site is mapped to a symbol, different symbols mean different domains
- `repr_spec`: if `domain_mapping` is an `AbstractMatrix{Symbol}`, this argument is a representation specifier (see `CoordinateRepr` docs for more information)
"""
function domains!(H::Matrix{ComplexF64}, domain_mapping::CoordinateRepr{Symbol})
    lattice_size = size(domain_mapping)
    E = zeros(2, 2)
    for i in 1:lattice_size[1]
        for j in 1:lattice_size[2]
            i1 = pair_to_index(lattice_size, i, j)
            i2 = pair_to_index(lattice_size, i+1, j)
            if domain_mapping[i1] != domain_mapping[i2]
                _set_hopping!(H, i1, i2, E)
            end
            i2 = pair_to_index(lattice_size, i, j+1)
            if domain_mapping[i1] != domain_mapping[i2]
                _set_hopping!(H, i1, i2, E)
            end
        end
    end
end

domains!(H::Matrix{ComplexF64}, domain_mapping::AbstractMatrix{Symbol}, repr_spec::Symbol) = 
    domains!(H, CoordinateRepr(domain_mapping, repr_spec))

@doc raw"""
    hamiltonian{T}((m_repr | m_lattice, repr_spec); <keyword arguments>)

Generates a Hamiltonian operator for a Chern insulator using the following formula

$\hat{H} = 
\sum_i m_i c^\dagger_i \sigma_z c_i + 
\sum_{x-links} c^\dagger_i \frac{\sigma_z - i \sigma_x}{2} c_j + 
\sum_{y-links} c^\dagger_i \frac{\sigma_z - i \sigma_y}{2} c_j + 
h. c.$

# Arguments
- `m_repr`: The value of `m` on different sites, in `CoordinateRepr` format
Alternatively, pass a matrix and a representation specifier (see `CoordinateRepr` for more information)
- `pbc`: Periodic boundary conditions. A `Tuple{Bool, Bool}`, each element sets boundary conditions for the horizontal and vertical edge respectively. Default is `(false, false)`
- `domains`: A matrix with elements of arbitrary type, which maps sites to isolated domains. The hopping members between different domains are erased. There are no isolated domains by default
- `field`: A function/lambda that takes two coordinates and returns the vector potential of the magnetic field. Used to calculate phase factors on hoppings. There is no magnetic field by default
"""
function hamiltonian(m_repr::CoordinateRepr{Float64}; kw...)
    arg_keys = Set(keys(kw))
    
    function process_kw(fun::Function, k::Symbol, type::Type)
        if k in arg_keys
            @assert kw[k] isa type "Unsupproted type of keyword \"$k\": expected '$type', got '$(typeof(kw[k]))'"
            fun(kw[k])
        end
    end

    pbc = (false, false)
    process_kw(:pbc, NTuple{2, Bool}) do _pbc
        pbc = _pbc
    end
    H = _hamiltonian(m_repr, pbc)

    process_kw(:domains, CoordinateRepr{Symbol}) do mapping
        domains!(H, mapping)
    end

    process_kw(:field, Function) do A_fun
        field!(H, A_fun, size(m_repr))
    end

    return H
end

hamiltonian(m_lattice::AbstractMatrix{<:Real}, repr_spec::Symbol; kw...) = 
    hamiltonian(CoordinateRepr(m_lattice, repr_spec); kw...)
