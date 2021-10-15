get_curr(operator, i, j) = operator[2 * i - 1:2 * i, 2 * j - 1:2 * j]

# LCM current macros

macro J_c(H, P, X, Y)
    quote 
        global J_type = "J_c"
        local C = 4π * $(esc(:($P * $X * $P * $Y * $P)))
        curr_c(i::Int, j::Int) = - real(tr(
            get_curr($(esc(H)), i, j) * get_curr(C, j, i) - get_curr(C, i, j) * get_curr($(esc(H)), j, i) 
        ))
    end
end

macro J_m(H, P, X, Y)
    quote
        global J_type = "J_m"
        Jpx = $(esc(:($P * $H * $X * $P)))
        Jpy = $(esc(:($P * $H * $Y * $P)))
        Xp = $(esc(:($P * $X * $P)))
        Yp = $(esc(:($P * $Y * $P)))
        curr_m(i::Int, j::Int) = 4π * real(tr(
            get_curr(Jpx, i, j) * get_curr(Yp, j, i) - get_curr(Yp, i, j) * get_curr(Jpx, j, i)
            - get_curr(Jpy, i, j) * get_curr(Xp, j, i) + get_curr(Xp, i, j) * get_curr(Jpy, j, i)
        ))
    end
end

macro J_eq(H, P, X, Y)
    quote
        global J_type = "J_eq"
        local Jpx = $(esc(:(($H * $P - $P * $H) * $X * $P)))
        local Jpy = $(esc(:(($H * $P - $P * $H) * $Y * $P)))
        local Xp = $(esc(:($P * $X * $P)))
        local Yp = $(esc(:($P * $Y * $P)))
        curr_nu(i::Int, j::Int) = - 4π * real(tr(
            get_curr(Jpx, i, j) * get_curr(Yp, j, i) - get_curr(Yp, i, j) * get_curr(Jpx, j, i)
            - get_curr(Jpy, i, j) * get_curr(Xp, j, i) + get_curr(Xp, i, j) * get_curr(Jpy, j, i)
        ))
    end
end

macro J_m_inv(H, P, X, Y)
    quote
        global J_type = "J_m_inv"
        local P = $(esc(P))
        local Jx =  $(esc(:($H * $X - $X * $H)))
        local Jy =  $(esc(:($H * $Y - $Y * $H)))
        local Xp = $(esc(:($P * $X * $P)))
        local Yp = $(esc(:($P * $Y * $P)))

        local Jxp = Jx * P
        local Jyp = Jy * P
        local Jyxp = Jy * Xp
        local Jxyp = Jx * Yp
        curr_m_inv(i::Int, j::Int) = 2π * real(tr(
            get_curr(P, i, j) * get_curr(Jxyp, j, i) - get_curr(Jxyp, i, j) * get_curr(P, j, i)
            - get_curr(Yp, i, j) * get_curr(Jxp, j, i) + get_curr(Jxp, i, j) * get_curr(Yp, j, i)
            - get_curr(P, i, j) * get_curr(Jyxp, j, i) + get_curr(Jyxp, i, j) * get_curr(P, j, i)
            + get_curr(Xp, i, j) * get_curr(Jyp, j, i) - get_curr(Jyp, i, j) * get_curr(Xp, j, i)
        ))
    end
end

macro J_best(H, P, X, Y)
    quote        
        global J_type = "J_best"
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
        curr_best(i::Int, j::Int) = - 4π * real(tr(
            get_curr(pdp, i, j) * get_curr(pxpyp, j, i) - get_curr(pxpyp, i, j) * get_curr(pdp, j, i)
            - get_curr(pdp, i, j) * get_curr(pxpyp', j, i) + get_curr(pxpyp', i, j) * get_curr(pdp, j, i)
            + get_curr(P, i, j) * get_curr(pdxpyp, j, i) - get_curr(pdxpyp, i, j) * get_curr(P, j, i)
            + get_curr(pdxp, i, j) * get_curr(pyp, j, i) - get_curr(pyp, i, j) * get_curr(pdxp, j, i)
            - get_curr(P, i, j) * get_curr(pdypxp, j, i) + get_curr(pdypxp, i, j) * get_curr(P, j, i)
            - get_curr(pdyp, i, j) * get_curr(pxp, j, i) + get_curr(pxp, i, j) * get_curr(pdyp, j, i)
        ))
    end
end

macro J_inv(H, P, X, Y)
    quote
        jc = $(esc(:(@J_c($H, $P, $X, $Y))))
        jminv = $(esc(:(@J_m_inv($H, $P, $X, $Y))))
        global J_type = "J_inv"
        curr_inv(i, j) = jc(i, j) + jminv(i, j)
    end
end

macro J(H, P, X, Y)
    quote
        jc = $(esc(:(@J_c($H, $P, $X, $Y))))
        jm = $(esc(:(@J_m($H, $P, $X, $Y))))
        global J_type = "J_mn"
        curr_inv(i, j) = jc(i, j) + jm(i, j)
    end
end