using LinearAlgebra

@doc raw"""
    coord_operators([lattice_size]; symmetric = true)

Returns a tuple of coordinate operators (i. e. $\hat{X}$ and $\hat{Y}$).

# Arguments
- `lattice_size`: the size of the lattice

# Keyword arguments
- `symmetric`: true if the operator is symmetrically defined (in other words, the central site of the lattice corresponds to `(0, 0)`), false otherwise. True by default
"""
function coord_operators(lattice_size::SizeType;
    symmetric::Bool = true)::NTuple{2,Matrix{Float64}}
    lattice_size = _current_lattice_size(lattice_size)
    len = prod(lattice_size)
    operators = []
    for axis in 1:2
        op_diagonal = fill(0., 2len)
        for i in 1:lattice_size[1], j in 1:lattice_size[2]
            site = pair_to_index(lattice_size, i, j)
            op_diagonal[2site - 1:2site] .= ((i, j)[3 - axis] - symmetric * (lattice_size[3 - axis] + 1) / 2)
        end
        push!(operators, Diagonal(op_diagonal))
    end
    return Tuple(operators)
end

coord_operators(; kw...) = coord_operators(nothing; kw...)

@doc raw"""
    filled_projector(H; <keyword arguments>)

Returns a projector onto the filled states (in other words - a density matrix of the ground state).

# Arguments
- `H`: the hamiltonian matrix

# Keyword arguments
- `fermi_energy`: the Fermi energy level. 0 by default
- `T`: temperature of the system, in energy units. 0 by default.
If value is non-zero, the state density will correspond to the Fermi-Dirac distribution
"""
function filled_projector(H::Matrix{<:Complex{<:Real}}; fermi_energy::Real = 0, T::Real = 0)
    val, vec = eigen(Hermitian(H))
    if T == 0
        d = @. (val ≤ fermi_energy) |> Int
    else
        d = @. 1 / (exp((val - fermi_energy) / T) + 1)
    end
    return Hermitian(vec * Diagonal(d) * vec')
end

"""
    currents(H, P[, lattice_size])

Returns a skew-symmetric matrix of electric currents between sites.

# Arguments
- `H`: the hamiltonian matrix
- `P`: the density matrix
- `lattice_size`: the size of the lattice the hamiltonian is defined for
"""
function currents(H::AbstractMatrix{Complex{T}}, P::AbstractMatrix{Complex{T}},
    lattice_size::SizeType=nothing) where T <: Real
    lattice_size = _current_lattice_size(lattice_size)
    currents_mat = zeros(prod(lattice_size), prod(lattice_size))
    for i in 1:prod(lattice_size), j in 1:(i - 1)
        if dist(lattice_size, i, j) == 1
            hop = H[2i - 1:2i, 2j - 1:2j]
            curr = 2tr(-im * hop * P[2j - 1:2j, 2i - 1:2i]) |> real
            currents_mat[i, j] = curr
            currents_mat[j, i] = -curr
        end
    end
    return currents_mat
end
