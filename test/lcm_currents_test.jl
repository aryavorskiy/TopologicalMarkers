@testset "LCM current macros" begin
    js_lambda::Function
    get_ham(_) = ham_b
    @evolution [
        :ham => get_ham => current_ham,
        P => get_ham => P_ev
    ] for t in time_domain
        js_lambda = @J_best current_ham, P_ev, X, Y
        js = [js_lambda(i, j) for i in 1:prod(siz), j in 1:prod(siz)]
        js2 = @currents js_lambda
        @test maximum(js + js') < 1e-10
        @test maximum(js - js2) < 1e-10
        global deriv_c = zeros(Float64, siz)
        for i in 1:prod(siz)
            local site = index_to_pair(siz, i)
            deriv_c[site[1], site[2]] = sum(js[:, i])
        end

        Jpx =  (P_ev * current_ham - current_ham * P_ev) * X * P_ev
        Jpy =  (P_ev * current_ham - current_ham * P_ev) * Y * P_ev
        Xp = P_ev * X * P_ev
        Yp = P_ev * Y * P_ev

        deriv = - 4π * (Jpx * Yp - Yp * Jpx + Xp * Jpy - Jpy * Xp)
        deriv_f = (deriv, siz)

        @test maximum(deriv_c - deriv_f) < 1e-10
    end
end