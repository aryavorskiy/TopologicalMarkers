using LinearAlgebra
using ProgressMeter

SIMPLIFY_EVOLUTION = false
TAYLOR_ORDER = 1
THRESHOLD = nothing
EVOL_CACHE = true

function _configure_exp!(simplify::Bool; kw...)
    global SIMPLIFY_EVOLUTION = simplify
    reset_cache()
    ks = Set(keys(kw))
    if :order in ks
        global TAYLOR_ORDER = kw[:order]::Real
        pop!(ks, :order)
    end
    
    if :threshold in ks
        global THRESHOLD = kw[:threshold]::Nullable{Real}
        pop!(ks, :threshold)
    end
    
    if :cache in ks
        global EVOL_CACHE = kw[:cache]::Bool
        pop!(ks, :cache)
    end
    
    if !isempty(ks)
        error("Unsupproted keywords: $(join(ks, ", "))")
    end
end

function taylor_exp(A::AbstractMatrix, k::Int)
    B = one(A) + A
    if k == 1
        return B
    end
    M = copy(A)
    for _ in 2:k
        M *= A / k
        B += M
    end
    return B
end

simple_exp(A::AbstractMatrix) = taylor_exp(A, TAYLOR_ORDER)

let cached_h::Any = nothing
    cached_t::Real = 0.1
    cached_e::Any = nothing
@doc raw"""
    evolution_operator(H, t)

Calculates the unitary evolution operator using the formula

$ \mathcal{U}(t) = e^{-\frac{1}{i\hbar} \hat{H} t} $

# Arguments
- `H`: the hamiltonian matrix
- `t`: the evolution time
"""
    global function evolution_operator(H::AbstractMatrix{ComplexF64}, t::Real)
        if t == cached_t && H == cached_h
            return cached_e
        end
        local evol = (SIMPLIFY_EVOLUTION && (THRESHOLD === nothing || t < THRESHOLD)) ? simple_exp((im * t) * H) : exp((im * t) * H)
        if EVOL_CACHE
            cached_h = H
            cached_t = t
            cached_e = evol
        end
        return evol
    end

    global function reset_cache()
        cached_h = nothing
        cached_t = 0.1
        cached_e = nothing
    end
end

function _expand(chain::Expr)
    @assert chain.head == :call
    @assert chain.args[1] == :(=>)
    @assert chain.args[3].head == :call
    @assert chain.args[3].args[1] == :(=>)
    local members = chain.args[2], chain.args[3].args[2], chain.args[3].args[3]
    @assert all((typeof(m) in (QuoteNode, Symbol)) for m in members)
    return members
end

"""
    @evolution [rules...] for_loop

Generates an environment with defined hamiltonian and density matrices that evolve by certain laws.
See [Unitary evolution](evolution.md) for more details.
"""
macro evolution(rules, loop)
    if typeof(loop) != Expr || loop.head != :for
        error("Expression type should be 'for', not '$(loop.head)'")
    end
    loop_iter, loop_body = loop.args
    loop_var, loop_range = loop_iter.args
    
    main_block = quote end
    p_evolvers = quote end
    h_evaluated = quote end
    
    if typeof(rules) != Expr || rules.head ∉ (:vect, :vcat, :hcat, :braces, :bracescat)
        error("Evolution specifier should be iterable, not '$(loop.head)'")
    end
    
    for statement in rules.args
        if typeof(statement) == LineNumberNode
            continue
        elseif typeof(statement) == Expr
            local a, b, c = _expand(statement)
            if a == :(:ham)
                local ham_fun, ham_var = b, c
                local h_init = :(local $(esc(ham_var)) = zero($(esc(ham_fun))(0)))
                local h_eval = :($(esc(ham_var)) = $(esc(ham_fun))($(esc(loop_var))))
                push!(main_block.args, h_init)
                push!(h_evaluated.args, h_eval)
            else
                local p_initial, ham_fun, p_target = a, b, c
                local p_init = :(local $(esc(p_target)) = copy($(esc(p_initial))))
                local p_evol = quote
                    local ev = evolution_operator($(esc(ham_fun))($(esc(loop_var))), dt)
                    $(esc(p_target)) = ev * $(esc(p_target)) * ev'
                end
                push!(main_block.args, p_init)
                append!(p_evolvers.args, p_evol.args)
            end
        end
    end

    push!(h_evaluated.args, :($(esc(loop_body))))
    local new_loop = quote
        local t_inner = 0.
        local len = length($(esc(loop_range)))
        p = Progress(len, dt=0.5, barglyphs = BarGlyphs("[=> ]"), showspeed = true)
        for $(esc(loop_var)) in $(esc(loop_range))
            local dt = $(esc(loop_var)) - t_inner
            $p_evolvers
            $h_evaluated
            ProgressMeter.next!(p)
            t_inner = $(esc(loop_var))
        end
    end
    append!(main_block.args, new_loop.args)
    return main_block
end