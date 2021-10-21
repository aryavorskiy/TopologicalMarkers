using LinearAlgebra

function filled_projector(H::Matrix{<:Complex{<:Real}}, energy_thr::Real=0)
    val, vec = eigen(Hermitian(H))
    d = (val .≤ energy_thr) .|> Int |> Diagonal
    return Hermitian(vec * d * vec')
end

function coord_operators(lattice_size::N{NTuple{2,Integer}}; 
    symmetric::Bool=true, type::Type{T}=Float64)::NTuple{2,Matrix{T}} where T <: Real
    lattice_size = _try_get_lattice_size(lattice_size)
    len = prod(lattice_size)
    operators = []
    for axis in 1:2
        op_diagonal = fill(zero(type), len * 2)
        for i in 1:lattice_size[1], j in 1:lattice_size[2]
            site = pair_to_index(lattice_size, i, j)
            op_diagonal[2 * site - 1:2 * site] .= ((i, j)[3 - axis] - symmetric * (lattice_size[3 - axis] + 1) / 2)
        end
        push!(operators, Diagonal(op_diagonal))
    end
    return Tuple(operators)
end

coord_operators(; kw...) = coord_operators(nothing; kw...)

function currents(H::AbstractMatrix{Complex{T}}, P::AbstractMatrix{Complex{T}}, lattice_size::N{NTuple{2,Integer}}=nothing) where T <: Real
    lattice_size = _try_get_lattice_size(lattice_size)
    currents_mat = zeros(prod(lattice_size), prod(lattice_size))
    for i in 1:prod(lattice_size), j in 1:(i - 1)
        if dist(lattice_size, i, j) == 1
            hop = H[2 * i - 1:2 * i, 2 * j - 1:2 * j]
            curr = 2 * tr(im * hop * P[2 * j - 1:2 * j, 2 * i - 1:2 * i]) |> real
            currents_mat[i, j] = curr
            currents_mat[j, i] = -curr
        end 
    end
    return currents_mat
end