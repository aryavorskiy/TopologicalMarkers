using LinearAlgebra: normalize

"""
    @landau(B)
Generates a function that returns the Landau gauge vector potential. This can be used as an argument for the `field!` function.

# Arguments
- `B`: the magnetic field value
"""
macro landau(B)
    quote
        @inline A(r::Vector{Float64})::Vector{Float64} = [0, $(esc(B)) * r[1]]
    end
end

"""
    @symm(B[, center])
Generates a function that returns the symmetric gauge vector potential. This can be used as an argument for the `field!` function.

# Arguments
- `B`: the value of magnetic field
- `center`: the center of symmetry of the vector potential
"""
macro symm(B, center=nothing)
    if center !== nothing
        return quote
            local c = [$(esc(center))...]
            @inline A(r::Vector{Float64})::Vector{Float64} = [-r[2] + c[2], r[1] - c[1]] * $(esc(B)) / 2
        end
    else 
        return quote
            local c = [(_current_lattice_size .- 1)...] / 2 .+ 1
            @inline A(r::Vector{Float64})::Vector{Float64} = [-r[2] + c[2], r[1] - c[1]] * $(esc(B)) / 2
        end
    end
end

"""
    flux(Φ[, center])
Generates a function that returns the vector potential for a flux quantum. This can be used as an argument for the `field!` function.

# Arguments
- `B`: the value of magnetic field
- `center`: the point where the flux is located
"""
macro flux(Φ, center=nothing)
    if center !== nothing
        return quote
            local c = [$(esc(center))...]
            @inline A(r::Vector{Float64})::Vector{Float64} = normalize([-r[2] + c[2], r[1] - c[1]]) * $(esc(Φ))
        end
    else 
        return quote
            local c = [(_current_lattice_size .- 1)...] / 2 .+ 1
            @inline A(r::Vector{Float64}) = normalize([-r[2] + c[2], r[1] - c[1]]) * $(esc(Φ))
        end
    end
end