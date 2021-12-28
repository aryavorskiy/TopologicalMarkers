using LinearAlgebra
using ProgressMeter

SIMPLIFY_EVOLUTION = false
TAYLOR_ORDER = 1
THRESHOLD = nothing

function _configure_evolution!(simplify::Bool; kw...)
    global SIMPLIFY_EVOLUTION = simplify
    ks = Set(keys(kw))
    if :order in ks
        global TAYLOR_ORDER = kw[:order]::Real
        pop!(ks, :order)
    end

    if :threshold in ks
        global THRESHOLD = kw[:threshold]::Nullable{Real}
        pop!(ks, :threshold)
    end

    if !isempty(ks)
        error("Unsupproted keyword(s): $(join(ks, ", "))")
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

@doc raw"""
    evolution_operator(H, t)

Calculates the unitary evolution operator using the formula

$ \mathcal{U}(t) = e^{-\frac{1}{i\hbar} \hat{H} t} $

# Arguments
- `H`: the hamiltonian matrix
- `t`: the evolution time
"""
evolution_operator(H::AbstractMatrix, t::Real) =
    (SIMPLIFY_EVOLUTION && (THRESHOLD === nothing || t < THRESHOLD)) ?
    simple_exp((im * t) * H) : exp((im * t) * H)

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

    hamiltonian_functions = Set()
    for statement in rules.args
        if typeof(statement) == LineNumberNode
            continue
        elseif typeof(statement) == Expr
            push!(hamiltonian_functions, _expand(statement)[2])
        end
    end
    for h in hamiltonian_functions
        h_eval = Symbol(string(h) * "_eval")
        h_eval_new = Symbol(string(h) * "_eval_new")
        local h_init = :(local $(esc(h_eval)) = nothing)
        local h_eval = :(local $(esc(h_eval_new)) = $(esc(h))($(esc(loop_var))))
        push!(main_block.args, h_init)
        push!(h_evaluated.args, h_eval)
    end

    for statement in rules.args
        if typeof(statement) == LineNumberNode
            continue
        elseif typeof(statement) == Expr
            local a, b, c = _expand(statement)
            if a == :(:ham)
                local ham_fun, ham_var = b, c
                h_eval_new = Symbol(string(ham_fun) * "_eval_new")
                local h_eval = :($(esc(ham_var)) = $(esc(h_eval_new)))
                push!(h_evaluated.args, h_eval)
            else
                local p_initial, ham_fun, p_target = a, b, c
                h_eval = Symbol(string(ham_fun) * "_eval")
                h_eval_new = Symbol(string(ham_fun) * "_eval_new")
                p_target_ev = Symbol(string(p_target) * "_evolutor")
                local p_init = quote
                    local $(esc(p_target)) = copy($(esc(p_initial)));
                    local $(esc(p_target_ev)) = nothing
                end
                local p_evol = quote
                    if $(esc(h_eval_new)) != $(esc(h_eval)) || dt != dt_old
                        $(esc(p_target_ev)) = evolution_operator($(esc(h_eval_new)), dt)
                    end
                    $(esc(h_eval)) = $(esc(h_eval_new))
                    $(esc(p_target)) = $(esc(p_target_ev)) * $(esc(p_target)) * adjoint($(esc(p_target_ev)))
                end
                append!(main_block.args, p_init.args)
                append!(p_evolvers.args, p_evol.args)
            end
        end
    end

    new_loop_body = quote end
    append!(new_loop_body.args, h_evaluated.args)
    append!(new_loop_body.args, p_evolvers.args)
    push!(new_loop_body.args, :($(esc(loop_body))))
    local new_loop = quote
        local t_inner = 0.
        local dt_old = 0.
        local len = length($(esc(loop_range)))
        p = Progress(len, dt=0.5, barglyphs = BarGlyphs("[=> ]"), showspeed = true)
        for $(esc(loop_var)) in $(esc(loop_range))
            local dt = $(esc(loop_var)) - t_inner
            $new_loop_body
            ProgressMeter.next!(p)
            t_inner = $(esc(loop_var))
            dt_old = dt
        end
    end
    append!(main_block.args, new_loop.args)
    return main_block
end
