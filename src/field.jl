macro landau(B, zone=:(:all, :all))
    zone_check = :(true)
    if zone isa Expr && zone.head != :tuple
        error("zone descriptor should be either a symbol or a tuple expression, $(zone.head) expression found")
    end
    for i in 1:2
        if !(zone isa Expr)
            zone_check = Expr(:call, :&&, zone_check, :($zone[$i][1] < r[$i] < $zone[$i][2]))
        elseif zone.args[i] isa Union{Symbol,Expr}
            tup = zone.args[i]
            zone_check = Expr(:call, :&&, zone_check, :($tup[1] < r[$i] < $tup[2]))
        elseif zone.args[i] != :(:all)
            x1, x2 = zone.args[i].args
            zone_check = Expr(:call, :&&, zone_check, :($x1 < r[$i] < $x2))
        end
    end
    quote
        A(x) = [0, $(esc(B)) * x[1], 0] * $(esc(zone_check))
    end
end

# TODO write macro
macro symm(B, center=nothing)
    if center !== nothing
        return quote
            local c = [$(esc(center))...] / 2
            A(x) = [-x[2] + c[2], x[1] - c[1], 0] * $(esc(B)) / 2
        end
    else 
        return quote
            local c = [_current_size...] / 2
            A(x) = [-x[2] + c[2], x[1] - c[1], 0] * $(esc(B)) / 2
        end
    end
end