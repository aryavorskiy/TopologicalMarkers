using Logging
import LinearAlgebra: eigen

function set_hopping!(H::AbstractMatrix{Complex{T}}, i::Integer, j::Integer, hop::AbstractMatrix{Complex{T}}) where T<:Real
    H[2 * i - 1:2 * i, 2 * j - 1:2 * j] = hop
    H[2 * j - 1:2 * j, 2 * i - 1:2 * i] = hop'
end

function _hamiltonian(type::Type{T}, m_lattice::CoordinateRepr, pbc::Tuple{Bool, Bool})::AbstractMatrix{Complex{T}} where T<:Real
    sz = size(m_lattice)
    H = zeros(Complex{type}, (prod(sz) * 2, prod(sz) * 2))

    # Generate diagonal elements
    for i in 1:prod(sz)
        H[2 * i - 1:2 * i, 2 * i - 1: 2 * i] = σ[3] * m_lattice[i]
    end

    # Generate hoppings
    for i in 1:sz[1], j in 1:sz[2]
        if i != sz[1] || pbc[1]
            i1 = pair_to_index(sz, i, j)
            i2 = pair_to_index(sz, i+1, j)
            set_hopping!(H, i1, i2, (σ[3] - im * σ[1]) / 2)
        end
        if j != sz[2] || pbc[2]
            i1 = pair_to_index(sz, i, j)
            i2 = pair_to_index(sz, i, j+1)
            set_hopping!(H, i1, i2, (σ[3] - im * σ[2]) / 2)
        end
    end
    return H
end

function field!(H::AbstractMatrix{Complex{<:Real}}, sz::NTuple{2, Integer}, A::Function, n_interval::Integer=10)
    local function peierls(i, j)
        phase = 0.
        r1 = index_to_pair(sz, i)
        r2 = index_to_pair(sz, j)
        dr = (r2 - r1) / n_interval
        for i in 1:n_interval
            phase += dr' * A(r1 + (i - 0.5) / n * r2)[1:2]
        end
        return 2π * im * phase |> exp
    end
    for i in 1:sz[1]
        for j in 1:sz[2]
            i1 = pair_to_index(sz, i, j)
            i2 = pair_to_index(sz, i+1, j)
            H[2 * i1 - 1:2 * i1, 2 * i2 - 1:2 * i2] .*= peierls(i1, i2)
            H[2 * i2 - 1:2 * i2, 2 * i1 - 1:2 * i1] .*= peierls(i2, i1)
            
            i2 = pair_to_index(sz, i, j+1)
            H[2 * i1 - 1:2 * i1, 2 * i2 - 1:2 * i2] .*= peierls(i1, i2)
            H[2 * i2 - 1:2 * i2, 2 * i1 - 1:2 * i1] .*= peierls(i2, i1)
        end
    end
end

function zones!(H::AbstractMatrix{Complex{<:Real}}, zone_mapping::CoordinateRepr)
    sz = size(zone_mapping)
    E = zeros(2, 2)
    for i in 1:sz[1]
        for j in 1:sz[2]
            i1 = pair_to_index(sz, i, j)
            i2 = pair_to_index(sz, i+1, j)
            if zone_mapping[i1] != zone_mapping[i2]
                set_hopping!(H, i1, i2, E)
            end
            i2 = pair_to_index(sz, i, j+1)
            if zone_mapping[i1] != zone_mapping[i2]
                set_hopping!(H, i1, i2, E)
            end
        end
    end
end

@doc raw"""
    hamiltonian(m_lattice; <keyword arguments>)

Generates a Hamiltonian operator for a Chern insulator using the following formula

$\hat{H} = 
\sum_i m_i c^\dagger_i \sigma_z c_i + 
\sum_{x-links} c^\dagger_i \frac{\sigma_z - i \sigma_x}{2} c_j + 
\sum_{y-links} c^\dagger_i \frac{\sigma_z - i \sigma_y}{2} c_j + 
h. c.$

# Arguments
- `m_lattice`: The value of `m` on different sites.
- `type`: Element type of the hamiltonian matrix. Default is `Float64`.
- `pbc`: Periodic boundary conditions. A `Tuple{Bool, Bool}`, each element sets boundary conditions for the horizontal and vertical edge respectively. Default is `(false, false)`.
- `zones`: A matrix with elements of arbitrary type, which maps sites to isolated zones. The hopping members between different zones are erased. There are no isolated zones by default.
- `field`: A function/lambda that takes two coordinates and returns the vector potential of the magnetic field. Used to calculate phase factors on hoppings. There is no magnetic field by default.
"""
function hamiltonian(m_lattice::CoordinateRepr{<:Real}; kwargs...)::AbstractMatrix{Complex{<:Real}}
    arg_keys = Set(keys(kwargs))
    
    function process_kw(fun::Function, k::Symbol, type::Type)
        if k in arg_keys
            @assert kwargs[k] isa type "Unsupproted type of keyword \"$k\": expected '$type', got '$(typeof(kwargs[k]))'"
            fun(kwargs[k])
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
    H = _hamiltonian(type, m_lattice, pbc)

    process_kw(:zones, CoordinateRepr) do mapping
        zones!(H, mapping)
    end

    process_kw(:field, Function) do A_fun
        field!(H, size(m_lattice),A_fun)
    end

    return H
end