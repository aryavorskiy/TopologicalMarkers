get_curr(operator, i, j) = operator[2 * i - 1:2 * i, 2 * j - 1:2 * j]

# Bianca-Resta current macros

@doc raw"""
    @J_c(H, P, X, Y)

Calculates the Bianca-Resta current part using the following formula:

$J_c(r, r') = \\ 4\pi (
\langle r | PXPYP | r' \rangle \langle r' | H | r \rangle -
\langle r | H | r' \rangle \langle r' | PXPYP | r \rangle
)$
"""
macro J_c(H, P, X, Y)
    quote
        global J_type = "J_c"
        local C = 4π * $(esc(:($P * $X * $P * $Y * $P)))
        curr_c(i::Int, j::Int) = - real(tr(
            get_curr($(esc(H)), i, j) * get_curr(C, j, i) - get_curr(C, i, j) * get_curr($(esc(H)), j, i)
        ))
    end
end

@doc raw"""
    @J_m_tr(H, P, X, Y)

Calculates the Bianca-Resta current part using the following formula:

$J_{m_{tr}}(r, r') = \\ 4\pi (
\langle r | P | r' \rangle \langle r' | [H,X] PYP | r \rangle -
\langle r | [H,X] PYP | r' \rangle \langle r' | P | r \rangle \\ -
\langle r | PYP | r' \rangle \langle r' | [H,X] P | r \rangle +
\langle r | [H,X] P | r' \rangle \langle r' | PYP | r \rangle \\ -
\langle r | P | r' \rangle \langle r' | [H,Y] PXP | r \rangle +
\langle r | [H,Y] PXP | r' \rangle \langle r' | P | r \rangle \\ +
\langle r | PXP | r' \rangle \langle r' | [H,Y] P | r \rangle -
\langle r | [H,Y] P | r' \rangle \langle r' | PXP | r \rangle
)$
"""
macro J_m_tr(H, P, X, Y)
    quote
        global J_type = "J_m_tr"
        local P = $(esc(P))
        local Jx =  $(esc(:($H * $X - $X * $H)))
        local Jy =  $(esc(:($H * $Y - $Y * $H)))
        local Xp = $(esc(:($P * $X * $P)))
        local Yp = $(esc(:($P * $Y * $P)))

        local Jxp = Jx * P
        local Jyp = Jy * P
        local Jyxp = Jy * Xp
        local Jxyp = Jx * Yp
        curr_m_tr(i::Int, j::Int) = 2π * real(tr(
            get_curr(P, i, j) * get_curr(Jxyp, j, i) - get_curr(Jxyp, i, j) * get_curr(P, j, i)
            - get_curr(Yp, i, j) * get_curr(Jxp, j, i) + get_curr(Jxp, i, j) * get_curr(Yp, j, i)
            - get_curr(P, i, j) * get_curr(Jyxp, j, i) + get_curr(Jyxp, i, j) * get_curr(P, j, i)
            + get_curr(Xp, i, j) * get_curr(Jyp, j, i) - get_curr(Jyp, i, j) * get_curr(Xp, j, i)
        ))
    end
end

@doc raw"""
    @J_eq(H, P, X, Y)

Calculates the Bianca-Resta current using the following formula:

$J_{eq}(r, r') = \\ 4\pi (
\langle r | [H,P] P | r' \rangle \langle r' | PXPYP | r \rangle -
\langle r | PXPYP | r' \rangle \langle r' | [H,P] P | r \rangle \\ -
\langle r | [H,P] P | r' \rangle \langle r' | PYPXP | r \rangle +
\langle r | PYPXP | r' \rangle \langle r' | [H,P] P | r \rangle \\ +
\langle r | P | r' \rangle \langle r' | [H,P] XPYP | r \rangle -
\langle r | [H,P] XPYP | r' \rangle \langle r' | P | r \rangle \\ +
\langle r | P [H,P] XP | r' \rangle \langle r' | PYP | r \rangle -
\langle r | PYP | r' \rangle \langle r' | P [H,P] XP | r \rangle \\ -
\langle r | P | r' \rangle \langle r' | [H,P] YPXP | r \rangle +
\langle r | [H,P] YPXP | r' \rangle \langle r' | P | r \rangle \\ -
\langle r | P [H,P] YP | r' \rangle \langle r' | PXP | r \rangle +
\langle r | PXP | r' \rangle \langle r' | P [H,P] YP | r \rangle
)$
"""
macro J_eq(H, P, X, Y)
    quote
        global J_type = "J_eq"
        local Pd = $(esc(:(($H * $P - $P * $H))))
        local P = $(esc(P))
        local X = $(esc(X))
        local Y = $(esc(Y))
        local pdp = Pd * P
        local pxpyp = P * X * P * Y * P
        local pxp = P * X * P
        local pyp = P * Y * P
        local pdxp = P * Pd * X * P
        local pdyp = P * Pd * Y * P
        local pdxpyp = pdxp * pyp
        local pdypxp = pdyp * pxp
        curr_eq(i::Int, j::Int) = - 4π * real(tr(
            get_curr(pdp, i, j) * get_curr(pxpyp, j, i) - get_curr(pxpyp, i, j) * get_curr(pdp, j, i)
            - get_curr(pdp, i, j) * get_curr(pxpyp', j, i) + get_curr(pxpyp', i, j) * get_curr(pdp, j, i)
            + get_curr(P, i, j) * get_curr(pdxpyp, j, i) - get_curr(pdxpyp, i, j) * get_curr(P, j, i)
            + get_curr(pdxp, i, j) * get_curr(pyp, j, i) - get_curr(pyp, i, j) * get_curr(pdxp, j, i)
            - get_curr(P, i, j) * get_curr(pdypxp, j, i) + get_curr(pdypxp, i, j) * get_curr(P, j, i)
            - get_curr(pdyp, i, j) * get_curr(pxp, j, i) + get_curr(pxp, i, j) * get_curr(pdyp, j, i)
        ))
    end
end

"""
    @J_loc(H, P, X, Y)

Calculates the sum of `@J_c` and `@J_m_tr` currents.
"""
macro J_loc(H, P, X, Y)
    quote
        jc = $(esc(:(@J_c($H, $P, $X, $Y))))
        jminv = $(esc(:(@J_m_inv($H, $P, $X, $Y))))
        global J_type = "J_loc"
        curr_inv(i, j) = jc(i, j) + jminv(i, j)
    end
end

"""
    @currents(currents_lambda)

Generates a matrix with currents, given a lambda/macro call that takes lattice site indices and returns the current value.

**Example usage:**

`currents_mat = @currents @J_b H P X Y`
"""
macro currents(call)
    if !(call isa Expr) || call.head != :macrocall
        return quote
            J = $(esc(call))
            [J(i, j) for i in 1:prod(CURRENT_LATTICE_SIZE), j in 1:prod(CURRENT_LATTICE_SIZE)]
        end
    end
    _size = :(prod(CURRENT_LATTICE_SIZE))

    for arg in (call.args[2:end])
        if arg isa Union{Symbol,Expr}
            _size = :(Base.round(Int, Base.size($arg)[1] / 2))
            break
        end
    end

    quote
        J = $(esc(call))
        [J(i, j) for i in 1:$(esc(_size)), j in 1:$(esc(_size))]
    end
end
