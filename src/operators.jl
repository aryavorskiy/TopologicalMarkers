using LinearAlgebra

function filled_projector(H::Matrix{<:Complex{<:Real}}, energy_thr::Real=0)
    val, vec = eigen(Hermitian(H))
    d = (val .≤ energy_thr) .|> Int |> Diagonal
    return Hermitian(vec * d * vec')
end

function coord_operators(sz::N{NTuple{2,Integer}}=nothing; 
    symmetric::Bool=true, type::Type{T}=Float64)::NTuple{2,Matrix{T}} where T <: Real
    sz = _try_get_sz(sz)
    len = prod(sz)
    operators = []
    for axis in 1:2
        op_diagonal = fill(zero(type), len * 2)
        for i in 1:sz[1], j in 1:sz[2]
            site = pair_to_index(sz, i, j)
            op_diagonal[2 * site - 1:2 * site] .= ((i, j)[3 - axis] - symmetric * (sz[3 - axis] + 1) / 2)
        end
        push!(operators, Diagonal(op_diagonal))
    end
    return Tuple(operators)
end

function currents(H::AbstractMatrix{Complex{T}}, P::AbstractMatrix{Complex{T}}, sz::N{NTuple{2,Integer}}=nothing) where T <: Real
    sz = _try_get_sz(sz)
    currents_mat = zeros(prod(sz), prod(sz))
    for i in 1:prod(sz), j in 1:(i - 1)
        if dist(sz, i, j) == 1
            hop = H[2 * i - 1:2 * i, 2 * j - 1:2 * j]
            curr = 2 * tr(im * hop * P[2 * j - 1:2 * j, 2 * i - 1:2 * i]) |> real
            currents_mat[i, j] = curr
            currents_mat[j, i] = -curr
        end 
    end
    return currents_mat
end