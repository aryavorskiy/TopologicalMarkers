using LinearAlgebra: norm, normalize
import Base:getindex, setindex!, size, +, -, *, /

# Coordinate representation

mutable struct CoordinateRepr{T}
    _inner_mat::Matrix{T}
end

function CoordinateRepr(lattice::Matrix{T}, repr_spec::Symbol) where T
    if repr_spec ∈ (:c, :coord)
        CoordinateRepr(copy(lattice))
    elseif repr_spec ∈ (:n, :natural)
        CoordinateRepr(permutedims(lattice[end:-1:1, :], (2, 1)) |> Matrix)
    else
        error("Unsupported indexing type $repr_spec")
    end
end

size(coord::CoordinateRepr) = size(coord._inner_mat)

getindex(coord::CoordinateRepr, args...) = getindex(coord._inner_mat, args...)
setindex!(coord::CoordinateRepr, args...) = setindex!(coord._inner_mat, args...)

@recipe f(::Type{T}, val::T) where T <: CoordinateRepr = val._inner_mat |> transpose

+(c1::CoordinateRepr, c2::CoordinateRepr) = CoordinateRepr(c1._inner_mat + c2._inner_mat)
-(c1::CoordinateRepr, c2::CoordinateRepr) = CoordinateRepr(c1._inner_mat - c2._inner_mat)
*(c1::CoordinateRepr, a::Number) = CoordinateRepr(c1._inner_mat * a)
/(c1::CoordinateRepr, a::Number) = CoordinateRepr(c1._inner_mat / a)

function heatmap_data(op::AbstractMatrix{<:Complex{T}}, lattice_size::N{NTuple{2,Integer}})::CoordinateRepr{T} where T <: Real
    lattice_size = _try_get_lattice_size(lattice_size)
    markers = zeros(T, lattice_size)
    for i in 1:prod(lattice_size)
        markers[i] = tr(op[2 * i - 1:2 * i, 2 * i - 1:2 * i]) |> real
    end
    return CoordinateRepr(markers, :c)
end

heatmap_data(op) = heatmap_data(op, nothing)

# Arrow visualization

function _arrow_data(lattice_size::NTuple{2,Integer}, cur::Real, i::Integer, j::Integer)
    if cur < 0
        i, j = j, i
        cur = -cur
    end
    local site = index_to_pair(lattice_size, i)
    local other = index_to_pair(lattice_size, j)
    local vec = normalize(other - site) * cur
    return Tuple(site), Tuple(vec)
end

function quiver_data(currents_mat::AbstractMatrix{<:Real}, lattice_size::N{NTuple{2,Integer}}=nothing; threshold::Real=0.1, dist_threshold::Real=Inf,
    xlims::NTuple{2, <:Real}=(-Inf, Inf), ylims::NTuple{2, <:Real}=(-Inf, Inf))
    lattice_size = _try_get_lattice_size(lattice_size)
    ps = Vector{NTuple{2,<:Integer}}()
    qs = Vector{NTuple{2,<:Real}}()
    for i in 1:prod(lattice_size), j in 1:(i - 1)
        p, q = _arrow_data(lattice_size, currents_mat[i, j], i, j)
        if norm(q) > threshold && dist(lattice_size, i, j) < dist_threshold && xlims[1] < p[1] < xlims[2] && ylims[1] < p[2] < ylims[2]
            push!(ps, p)
            push!(qs, q)
        end 
    end
    return ps, qs
end

