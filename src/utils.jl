using RecipesBase

N{T} = Union{T,Nothing}
SizeType = N{NTuple{2,Int}}

# Lattice size stored here

CURRENT_LATTICE_SIZE = nothing

function _try_get_lattice_size(lattice_size::SizeType) 
    if lattice_size !== nothing
        return lattice_size
    end
    if CURRENT_LATTICE_SIZE !== nothing
        return CURRENT_LATTICE_SIZE
    else
        error("Please specify lattice size explicitly")
    end
end

function _set_lattice_size!(lattice_size::NTuple{2, Int})
    global CURRENT_LATTICE_SIZE = lattice_size
end

_set_lattice_size!(a::Int, b::Int) = _set_lattice_size!((a, b))

# Pauli matrices
const σ_x = [0 1; 1 0]
const σ_y = [0 -1im; 1im 0]
const σ_z = [1 0; 0 -1]
const E = [1 0; 0 1]
const σ = [σ_x, σ_y, σ_z]

# Pair-index conversion

pair_to_index(lattice_size::NTuple{2,Int}, a::Int, b::Int) = (b - 1) % lattice_size[2] * lattice_size[1] + (a - 1) % lattice_size[1] + 1

pair_to_index(lattice_size::NTuple{2,Int}, pair::Vector{Int}) = pair_to_index(lattice_size, pair[1], pair[2])

function index_to_pair(lattice_size::NTuple{2,Int}, i::Int)::Vector{Int}
    a = i % lattice_size[1] == 0 ? lattice_size[1] : i % lattice_size[1]
    return [a, round(Int, (i - a) / lattice_size[1]) + 1]
end

dist(lattice_size, i, j) = norm(index_to_pair(lattice_size, i) - index_to_pair(lattice_size, j))

function adjacent_sites(lattice_size::NTuple{2,Int}, site::Vector{Int}, order::Int)
    adj = Set()
    for i in 0:order, sig1 in (-1, 1), sig2 in (-1, 1)
        new_site = @. site + [i * sig1, order - i] * sig2
        if all(@. 0 < new_site < lattice_size)
            push!(adj, new_site)
        end
    end
    return adj
end