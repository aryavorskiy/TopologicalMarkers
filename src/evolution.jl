using LinearAlgebra

let store = Dict()
    global function evolution_operator(H::AbstractMatrix{<:Complex{<:Real}}, t::Real)
        if length(store) > 200
            store = Dict()
        end
        if (H, t) ∉ keys(store)
            store[(H, t)] = exp(im * H * t)
        end
        return store[(H, t)]
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

macro evolution(evolvers, loop)
    if typeof(loop) != Expr || loop.head != :for
        error("Expression type should be 'for', not '$(loop.head)'")
    end
    loop_iter, loop_body = loop.args
    loop_var, loop_range = loop_iter.args
    
    main_block = quote end
    p_evolvers = quote end
    h_evaluated = quote end
    
    if typeof(evolvers) != Expr || evolvers.head ∉ (:vect, :vcat, :hcat, :braces, :bracescat)
        error("Evolution specifier should be iterable, not '$(loop.head)'")
    end
    
    for statement in evolvers.args
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

    # Cleanup this

    push!(h_evaluated.args, :($(esc(loop_body))))
    local new_loop = quote
        local t_inner = 0.
        local len = length($(esc(loop_range)))
        local l = 60
        local counter = 0
        for $(esc(loop_var)) in $(esc(loop_range))
            pbar(counter / len, l)
            local dt = $(esc(loop_var)) - t_inner
            $p_evolvers
            $h_evaluated

            t_inner = $(esc(loop_var))
            counter += 1
        end
        print("\rProcessing complete!" * " "^l * "\r")
    end
    append!(main_block.args, new_loop.args)
    return main_block
end