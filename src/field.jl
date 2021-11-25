using LinearAlgebra: normalize
using StaticArrays

"""
    @landau(B)
Generates a function that returns the Landau gauge vector potential. This can be used as an argument for the `field!` function.

# Arguments
- `B`: the magnetic field value
"""
macro landau(B)
    quote
        @inline A(r::SVector{2, Float64}) = SA[0, $(esc(B)) * r[1], 0]
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
            @inline A(r::SVector{2, Float64}) = SA[-r[2] + c[2], r[1] - c[1], 0] * $(esc(B)) / 2
        end
    else 
        return quote
            @inline A(r::SVector{2, Float64}) = (
                local c = [(_try_get_lattice_size(nothing) .- 1)...] / 2 .+ 1; 
                SA[-r[2] + c[2], r[1] - c[1], 0] * $(esc(B)) / 2
            )
        end
    end
end

"""
    flux(Φ[, point])
Generates a function that returns the vector potential for a flux quantum. This can be used as an argument for the `field!` function.

# Arguments
- `Φ`: the value of magnetic field
- `point`: the point where the flux is located
"""
macro flux(Φ, point=nothing)
    if point !== nothing
        return quote
            local c = [$(esc(point))...]
            @inline A(r::SVector{2, Float64}) = normalize(SA[-r[2] + c[2], r[1] - c[1], 0]) * $(esc(Φ))
        end
    else 
        return quote
            @inline A(r::SVector{2, Float64}) = (
                local c = [(_try_get_lattice_size(nothing) .- 1)...] / 2 .+ 1;
                normalize(SA[-r[2] + c[2], r[1] - c[1], 0]) * $(esc(Φ))
            )
        end
    end
end