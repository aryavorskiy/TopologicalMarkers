function index_to_pair(lattice_size::NTuple{2,Integer}, i::Integer)::Vector{Integer}
    a = i % lattice_size[1] == 0 ? lattice_size[1] : i % lattice_size[1]
    return [a, round(Integer, (i - a) / lattice_size[1]) + 1]
end

@testset "LCM current macros" begin
    siz = (15, 15)
    ms = CoordinateRepr(ones(siz))
    ms2 = CoordinateRepr(ones(siz) * 3)
    ham = hamiltonian(ms)
    ham2 = hamiltonian(ms2)
    get_ham(_) = ham2
    @evolution [
        :ham => get_ham => current_ham,
        P => get_ham => P_ev
    ] for t in 0:5:30
        js_lambda = @J_best current_ham P_ev X Y
        js = [js_lambda(i, j) for i in 1:prod(siz), j in 1:prod(siz)]
        js2 = @currents js_lambda
        @test maximum(js + js') < 1e-12
        @test maximum(js - js2) < 1e-12
        global deriv_c = zeros(Float64, siz)
        for i in 1:prod(siz)
            local site = index_to_pair(siz, i)
            deriv_c[site...] = sum(js[:, i])
        end

        Jpx =  (P_ev * current_ham - current_ham * P_ev) * X * P_ev
        Jpy =  (P_ev * current_ham - current_ham * P_ev) * Y * P_ev
        Xp = P_ev * X * P_ev
        Yp = P_ev * Y * P_ev

        deriv = - 4π * (Jpx * Yp - Yp * Jpx + Xp * Jpy - Jpy * Xp)
        deriv_f = heatmap_data(deriv)

        @test maximum(deriv_c - deriv_f._inner_mat) < 1e-10
    end
end