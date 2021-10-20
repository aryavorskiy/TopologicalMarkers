using LinearAlgebra: norm, normalize
import Base:getindex, setindex!, size, +, -, *, /

# Coordinate representation

mutable struct CoordinateRepr{T}
    _inner_mat::Matrix{T}
end

function CoordinateRepr(mat::Matrix{T}, indexing_type::Symbol) where T <: Real
    if indexing_type ∈ (:c, :coord)
        CoordinateRepr(copy(mat))
    elseif indexing_type ∈ (:n, :natural)
        CoordinateRepr(mat[end:-1:1, :] |> transpose |> Matrix)
    else
        error("Unsupported indexing type $indexing_type")
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

function heatmap_data(op::AbstractMatrix{<:Complex{T}}, sz::N{NTuple{2,Integer}}=nothing)::CoordinateRepr{T} where T <: Real
    sz = _try_get_sz(sz)
    markers = zeros(T, sz)
    for i in 1:prod(sz)
        markers[i] = tr(op[2 * i - 1:2 * i, 2 * i - 1:2 * i]) |> real
    end
    return CoordinateRepr(markers, :c)
end

# Arrow visualization

function _arrow_data(sz::NTuple{2,Integer}, cur::Real, i::Integer, j::Integer)
    if cur < 0
        i, j = j, i
        cur = -cur
    end
    local site = index_to_pair(sz, i)
    local other = index_to_pair(sz, j)
    local vec = normalize(other - site) * cur
    return Tuple(site), Tuple(vec)
end

function quiver_data(currents_mat::AbstractMatrix{<:Real}, sz::N{NTuple{2,Integer}}=nothing; threshold::Real=0.1, dist_threshold::Real=Inf)
    sz = _try_get_sz(sz)
    ps = Vector{NTuple{2,<:Integer}}()
    qs = Vector{NTuple{2,<:Real}}()
    for i in 1:prod(sz), j in 1:(i - 1)
        p, q = _arrow_data(sz, currents_mat[i, j], i, j)
        if norm(q) > threshold && dist(sz, i, j) < dist_threshold
            push!(ps, p)
            push!(qs, q)
        end 
    end
    return ps, qs
end

