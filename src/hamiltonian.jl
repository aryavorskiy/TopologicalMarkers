using Logging
import LinearAlgebra: eigen

_current_lattice_size = nothing

function _try_get_lattice_size(lattice_size::N{NTuple{2, Integer}}) 
    if lattice_size !== nothing
        return lattice_size
    end
    if _current_lattice_size !== nothing
        return _current_lattice_size
    else
        error("Please specify lattice size explicitly")
    end
end

function set_hopping!(H::AbstractMatrix{Complex{T}}, i::Integer, j::Integer, hop::AbstractMatrix) where T<:Real
    H[2 * i - 1:2 * i, 2 * j - 1:2 * j] = hop
    H[2 * j - 1:2 * j, 2 * i - 1:2 * i] = hop'
end

function _hamiltonian(dtype::Type{T}, m_repr::CoordinateRepr, pbc::Tuple{Bool, Bool})::AbstractMatrix{Complex{T}} where T<:Real
    lattice_size = size(m_repr)
    global _current_lattice_size = lattice_size
    H = zeros(Complex{dtype}, (prod(lattice_size) * 2, prod(lattice_size) * 2))

    # Generate diagonal elements
    for i in 1:prod(lattice_size)
        H[2 * i - 1:2 * i, 2 * i - 1: 2 * i] = σ[3] * m_repr[i]
    end

    # Generate hoppings
    for i in 1:lattice_size[1], j in 1:lattice_size[2]
        if i != lattice_size[1] || pbc[1]
            i1 = pair_to_index(lattice_size, i, j)
            i2 = pair_to_index(lattice_size, i+1, j)
            set_hopping!(H, i1, i2, (σ[3] - im * σ[1]) / 2)
        end
        if j != lattice_size[2] || pbc[2]
            i1 = pair_to_index(lattice_size, i, j)
            i2 = pair_to_index(lattice_size, i, j+1)
            set_hopping!(H, i1, i2, (σ[3] - im * σ[2]) / 2)
        end
    end
    return H
end

@doc raw"""
    field!(H, A[, lattice_size]; intervals = 10)

Applies magnetic field to specified hamiltonian. In other words, all hoppings get multiplied on a specific phase factor that can be calculated using Peierls substitution:

$ \varphi_{ij} = \frac{1}{2\pi} \int_i^j A(r) \cdot dr $

This integral is calculated explicitly for every hopping, using the `A` function.

# Arguments
- `H`: the hamiltonian matrix.
- `A`: a function that takes a `Vector` representing a point and returns a `Vector` representing the vector potential in that point.
- `lattice_size`: the size of the lattice the hamiltonian is defined for. If not provided, this function will use the value for the hamiltonian matrix that was created last.
- `intervals`: the number of intervals to use when calculating the Peierls substitution phase factor
"""
function field!(H::AbstractMatrix{<:Complex}, A::Function, lattice_size::N{NTuple{2, Integer}} = nothing; intervals::Integer=10)
    lattice_size = _try_get_lattice_size(lattice_size)
    local function peierls(i, j)
        phase = 0.
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
        H[2 * i1 - 1:2 * i1, 2 * i2 - 1:2 * i2] *= peierls(i1, i2)
        H[2 * i2 - 1:2 * i2, 2 * i1 - 1:2 * i1] *= peierls(i2, i1)
        
        i2 = pair_to_index(lattice_size, i, j+1)
        H[2 * i1 - 1:2 * i1, 2 * i2 - 1:2 * i2] *= peierls(i1, i2)
        H[2 * i2 - 1:2 * i2, 2 * i1 - 1:2 * i1] *= peierls(i2, i1)
    end
end

"""
    zones!(H, zone_mapping[, repr_spec])

Divides the Chern insulator hamiltonian into several unconnected zones. The hoppings between these zones are erased.

# Arguments
- `H`: the hamiltonian matrix.
- `zone_mapping`: an `AbstractMatrix{Symbol}` or `CoordinateRepr{Symbol}`. Each site is mapped to a symbol, different symbols mean different zones.
- `repr_spec`: if `zone_mapping` is an `AbstractMatrix{Symbol}`, this argument is a representation specifier (see `CoordinateRepr` docs for more information).
"""
function zones!(H::AbstractMatrix{<:Complex{<:Real}}, zone_mapping::CoordinateRepr{Symbol})
    lattice_size = size(zone_mapping)
    E = zeros(2, 2)
    for i in 1:lattice_size[1]
        for j in 1:lattice_size[2]
            i1 = pair_to_index(lattice_size, i, j)
            i2 = pair_to_index(lattice_size, i+1, j)
            if zone_mapping[i1] != zone_mapping[i2]
                set_hopping!(H, i1, i2, E)
            end
            i2 = pair_to_index(lattice_size, i, j+1)
            if zone_mapping[i1] != zone_mapping[i2]
                set_hopping!(H, i1, i2, E)
            end
        end
    end
end

zones!(H::AbstractMatrix{<:Complex}, zone_mapping::AbstractMatrix{Symbol}, repr_spec::Symbol) = 
    zones!(H, CoordinateRepr(zone_mapping, repr_spec))

@doc raw"""
    hamiltonian((m_repr | m_lattice, repr_spec); <keyword arguments>)

Generates a Hamiltonian operator for a Chern insulator using the following formula

$\hat{H} = 
\sum_i m_i c^\dagger_i \sigma_z c_i + 
\sum_{x-links} c^\dagger_i \frac{\sigma_z - i \sigma_x}{2} c_j + 
\sum_{y-links} c^\dagger_i \frac{\sigma_z - i \sigma_y}{2} c_j + 
h. c.$

# Arguments
- `m_repr`: The value of `m` on different sites, in `CoordinateRepr` format. 
Alternatively, pass a matrix and a representation specifier (see `CoordinateRepr` for more information).
- `type`: Element type of the hamiltonian matrix. Default is `Float64`.
- `pbc`: Periodic boundary conditions. A `Tuple{Bool, Bool}`, each element sets boundary conditions for the horizontal and vertical edge respectively. Default is `(false, false)`.
- `zones`: A matrix with elements of arbitrary type, which maps sites to isolated zones. The hopping members between different zones are erased. There are no isolated zones by default.
- `field`: A function/lambda that takes two coordinates and returns the vector potential of the magnetic field. Used to calculate phase factors on hoppings. There is no magnetic field by default.
"""
function hamiltonian(m_repr::CoordinateRepr{<:Real}; kw...)
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

    type = Float64
    process_kw(:dtype, Type{<:Real}) do _type
        type = _type
    end
    H = _hamiltonian(type, m_repr, pbc)

    process_kw(:zones, CoordinateRepr) do mapping
        zones!(H, mapping)
    end

    process_kw(:field, Function) do A_fun
        field!(H, A_fun, size(m_repr))
    end

    return H
end

hamiltonian(m_lattice::AbstractMatrix{<:Real}, repr_spec::Symbol; kw...) = 
    hamiltonian(CoordinateRepr(m_lattice, repr_spec); kw...)