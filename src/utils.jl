using RecipesBase

N{T} = Union{T,Nothing}

# Pauli matrices
const σ_x = [0 1; 1 0]
const σ_y = [0 -1im; 1im 0]
const σ_z = [1 0; 0 -1]
const E = [1 0; 0 1]
const σ = [σ_x, σ_y, σ_z]

# Pair-index conversion

pair_to_index(sz::NTuple{2,Integer}, a::Integer, b::Integer) = (b - 1) % sz[2] * sz[1] + (a - 1) % sz[1] + 1

pair_to_index(sz::NTuple{2,Integer}, pair::Vector{Integer}) = pair_to_index(sz, pair[1], pair[2])

function index_to_pair(sz::NTuple{2,Integer}, i::Integer)::Vector{Integer}
    a = i % sz[1] == 0 ? sz[1] : i % sz[1]
    return [a, round(Integer, (i - a) / sz[1]) + 1]
end

dist(sz, i, j) = norm(index_to_pair(sz, i) - index_to_pair(sz, j))

function adjacent_sites(sz::NTuple{2,Integer}, site::Vector{Integer}, order::Integer)
    adj = Set()
    for i in 0:order, sig1 in (-1, 1), sig2 in (-1, 1)
        new_site = @. site + [i * sig1, order - i] * sig2
        if all(@. 0 < new_site < sz)
            push!(adj, new_site)
        end
    end
    return adj
end

# Progressbar
function pbar(p::Real, len::Integer=60; fill::Char='#', empty::Char='-')
    fill_len = round(Int, p * len)
    return fill^fill_len * empty^(len - fill_len)
end