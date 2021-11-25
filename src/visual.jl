using LinearAlgebra: norm, normalize
using StaticArrays: SVector
import Base:getindex, setindex!, size, +, -, *, /

# Coordinate representation

mutable struct CoordinateRepr{T}
    _inner_mat::Matrix{T}
end

"""
    CoordinateRepr{T}

A wrapper struct for matrices to be conveniently plotted.

# Arguments
- `lattice`: A matrix representing some quantity defined on the lattice (e. g. the LCM)
- `repr_spec`: A symbol defining the way how the lattice sites match the matrix values:
    - `:c` or `:coord`: The value for `(x, y)` site is `A[x, y]`. This is the default value
    - `:n` of `:natural`: If you print out the matrix as you usually do, and then imagine
    a coordinate system with its center in the bottom-left corner, this will be the mapping between
    sites and matrix values
"""
function CoordinateRepr(lattice::AbstractArray{T}, repr_spec::Symbol) where T
    if repr_spec ∈ (:c, :coord)
        CoordinateRepr(copy(lattice |> Array))
    elseif repr_spec ∈ (:n, :natural)
        CoordinateRepr(permutedims(Array(lattice)[end:-1:1, :], (2, 1)))
    else
        error("Unsupported indexing type '$repr_spec'")
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

@doc raw"""
    heatmap_data(op[, lattice_size])

Generates a CoordinateRepr for coordinate representation traces:
$\langle r | \hat{\mathcal{O}} | r \rangle$.

# Arguments
- `op`: the operator to find values for
- `lattice_size`: the size of the lattice
"""
function heatmap_data(op::AbstractMatrix{Complex{T}},
    lattice_size::SizeType)::CoordinateRepr{T} where T <: Real
    lattice_size = _try_get_lattice_size(lattice_size)
    markers = zeros(T, lattice_size)
    for i in 1:prod(lattice_size)
        markers[i] = tr(op[2 * i - 1:2 * i, 2 * i - 1:2 * i]) |> real
    end
    return CoordinateRepr(markers, :c)
end

heatmap_data(op) = heatmap_data(op, nothing)

# Arrow visualization

function _arrow_data(lattice_size::NTuple{2,Int}, cur::Real, i::Int, j::Int)
    if cur < 0
        i, j = j, i
        cur = -cur
    end
    local site::SVector = index_to_pair(lattice_size, i)
    local other::SVector = index_to_pair(lattice_size, j)
    local vec = normalize(other - site) * cur
    return site.data, vec.data
end

"""
    quiver_data(currents_mat[, lattice_size]; <keyword arguments>)

Generates data for a quiver plot using a matrix with currents.
The output is a tuple of two vectors with equal length: one contains arrow origins, the other one - arrow vectors.

# Arguments
- `currents_mat`: a matrix with currents
- `lattice_size`: the size of the lattice

# Keyword arguments
- `threshold`: minimum value of the current to be put to output. Default is `0.1`
- `dist_threshold`: maximum distance between sites for which the current will be evaluated. Infinite by default
- `xlims` and `ylims`: limit the area in which the currents will be evaluated. Infinite by default
"""
function quiver_data(currents_mat::AbstractMatrix{<:Real}, lattice_size::SizeType=nothing;
    threshold::Real=0.1, dist_threshold::Real=Inf, xlims::NTuple{2, <:Real}=(-Inf, Inf),
    ylims::NTuple{2, <:Real}=(-Inf, Inf))
    lattice_size = _try_get_lattice_size(lattice_size)
    ps = NTuple{2,<:Int}[]
    qs = NTuple{2,<:Real}[]
    for i in 1:prod(lattice_size), j in 1:(i - 1)
        if abs(currents_mat[i, j]) < threshold || dist(lattice_size, i, j) > dist_threshold
            continue
        end
        p, q = _arrow_data(lattice_size, currents_mat[i, j], i, j)
        if xlims[1] < p[1] < xlims[2] && ylims[1] < p[2] < ylims[2]
            push!(ps, p)
            push!(qs, q)
        end
    end
    return ps, qs
end
