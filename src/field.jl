using LinearAlgebra: normalize

"""
    @landau(B)
Generates a function that returns the Landau gauge vector potential.

# Arguments
- `B`: the magnetic field value
"""
macro landau(B)
    quote
        A(r) = [0, $(esc(B)) * r[1], 0]
    end
end

"""
    @symm(B[, center])
Generates a function that returns the symmetric gauge vector potential.

# Arguments
- `B`: the value of magnetic field
- `center`: the center of symmetry of the vector potential
"""
macro symm(B, center=nothing)
    if center !== nothing
        return quote
            local c = [$(esc(center))...] / 2
            A(r) = [-r[2] + c[2], r[1] - c[1], 0] * $(esc(B)) / 2
        end
    else 
        return quote
            local c = [_current_lattice_size...] / 2
            A(r) = [-r[2] + c[2], r[1] - c[1], 0] * $(esc(B)) / 2
        end
    end
end

"""
    flux(Φ[, center])
Generates a function that returns the vector potential for a flux quantum.

# Arguments
- `B`: the value of magnetic field
- `center`: the point where the flux is located
"""
macro flux(Φ, center=nothing)
    if center !== nothing
        return quote
            local c = [$(esc(center))...] / 2
            A(r) = normalize([-r[2] + c[2], r[1] - c[1], 0]) * $(esc(Φ)) / 2
        end
    else 
        return quote
            local c = [_current_lattice_size...] / 2
            A(r) = normalize([-r[2] + c[2], r[1] - c[1], 0]) * $(esc(Φ)) / 2
        end
    end
end