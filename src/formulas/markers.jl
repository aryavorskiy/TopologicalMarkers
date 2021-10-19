using Logging:@warn
using LinearAlgebra:eigen, Hermitian

local_chern(P, X, Y) = -4π * im * P * X * P * Y * P

function streda_p(H, Hb, P, B)
    local Q = I - P
    local Hc = Q * H * Q + P * H * P
    local dH = Hb - Hc
    local vec, val = eigen(Hc)
    local dP = zeros(ComplexF64, P |> size)
    local rhs = vec' * (dH * P - P * dH) * vec
    len = round(Int64, size(P)[1] / 2)

    # Check gap
    gap = abs(val[len] - val[len + 1])
    if gap < 1e-3
        @warn("small gap size: $gap")
    end
    
    for i in 1:len, j in len + 1:2len
        dP[i, j] = rhs[i, j] / (val[j] - val[i])
        dP[j, i] = rhs[j, i] / (val[i] - val[j])
    end
    return vec * Hermitian(dP) * vec' / B
end