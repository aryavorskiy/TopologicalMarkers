"""
    @landau(B[, zone])
Generates a function that returns the Landau gauge vector potential.

# Arguments
- `B`: the magnetic field value
- `zone`: a `NTuple{2, Any}` which describes the part of the space where the field is present.
"""
macro landau(B, zone=:(:all, :all))
    zone_check = :(true)
    if zone isa Expr && zone.head != :tuple
        error("zone descriptor should be either a symbol or a tuple expression, $(zone.head) expression found")
    end
    for i in 1:2
        if !(zone isa Expr)
            zone_check = Expr(:call, :&, zone_check, :($(esc(zone))[$i][1] < r[$i] < $(esc(zone))[$i][2]))
        elseif zone.args[i] isa Symbol || (zone.args[i] isa Expr && zone.args[i].head != :tuple)
            tup = zone.args[i]
            zone_check = Expr(:call, :&, zone_check, :($(esc(tup))[1] < r[$i] < $(esc(tup))[2]))
        elseif zone.args[i] != :(:all)
            x1, x2 = zone.args[i].args
            zone_check = Expr(:call, :&, zone_check, :($(esc(x1)) < r[$i] < $(esc(x2))))
        end
    end
    quote
        A(r) = [0, $(esc(B)) * r[1], 0] * $zone_check
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