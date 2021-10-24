using RecipesBase

N{T} = Union{T,Nothing}
SizeType = N{NTuple{2,Int}}

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

# Progressbar
let previous_pct = -1.
    global function pbar(progress::Real, len::Int=60; fill::Char='#', empty::Char='-')
        fill_len = round(Int, progress * len)
        pb = fill^fill_len * empty^(len - fill_len)
        pct = round(progress * 100, digits=1)
        if pct != previous_pct
            print("\r [$pb] $pct%")
            previous_pct = pct
        end
    end
end