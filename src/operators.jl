using LinearAlgebra

@doc raw"""
    coord_operators([lattice_size]; <keyword arguments>)

Returns a tuple of coordinate operators (e. g. $\hat{X}$ and $\hat{Y}$).

# Arguments
- `lattice_size`: the size of the lattice
- `symmetric`: defines if the operators are symmetricall (e. g. the center of the lattice corresponds to `(0, 0)`)
"""
function coord_operators(lattice_size::SizeType; 
    symmetric::Bool=true)::NTuple{2,Matrix{Float64}}
    lattice_size = _try_get_lattice_size(lattice_size)
    len = prod(lattice_size)
    operators = []
    for axis in 1:2
        op_diagonal = fill(0., len * 2)
        for i in 1:lattice_size[1], j in 1:lattice_size[2]
            site = pair_to_index(lattice_size, i, j)
            op_diagonal[2 * site - 1:2 * site] .= ((i, j)[3 - axis] - symmetric * (lattice_size[3 - axis] + 1) / 2)
        end
        push!(operators, Diagonal(op_diagonal))
    end
    return Tuple(operators)
end

coord_operators(; kw...) = coord_operators(nothing; kw...)

"""
    filled_projector(H[, energy_thr=0])

Returns a projector onto the filled states (in other words, a ground state density matrix).

# Arguments
- `H`: the hamiltonian matrix
- `energy_thr`: the Fermi energy level
"""
function filled_projector(H::Matrix{<:Complex{<:Real}}, energy_thr::Real=0)
    val, vec = eigen(Hermitian(H))
    d = (val .≤ energy_thr) .|> Int |> Diagonal
    return Hermitian(vec * d * vec')
end

"""
    currents(H, P[, lattice_size])

Returns a skew-symmetric matrix of electric currents between sites.

# Arguments
- `H`: the hamiltonian matrix.
- `P`: the density matrix.
- `lattice_size`: the size of the lattice the hamiltonian is defined for.
"""
function currents(H::AbstractMatrix{Complex{T}}, P::AbstractMatrix{Complex{T}}, lattice_size::SizeType=nothing) where T <: Real
    lattice_size = _try_get_lattice_size(lattice_size)
    currents_mat = zeros(prod(lattice_size), prod(lattice_size))
    for i in 1:prod(lattice_size), j in 1:(i - 1)
        if dist(lattice_size, i, j) == 1
            hop = H[2 * i - 1:2 * i, 2 * j - 1:2 * j]
            curr = 2 * tr(-im * hop * P[2 * j - 1:2 * j, 2 * i - 1:2 * i]) |> real
            currents_mat[i, j] = curr
            currents_mat[j, i] = -curr
        end 
    end
    return currents_mat
end