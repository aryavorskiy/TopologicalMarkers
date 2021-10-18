using Plots: plot, plot!, Plot, heatmap!, quiver!, xlims!, ylims!, hline!, vline!, current, @layout, grid, EmptyLayout, pct
using LinearAlgebra: norm, normalize
using Base:getindex

function heatmap_data(op::AbstractMatrix{<:Complex{T}}, sz::N{NTuple{2,Integer}}=nothing)::CoordinateRepr{T} where T <: Real
    sz = _try_get_sz(sz)
    markers = zeros(T, sz)
    for i in 1:prod(sz)
        markers[i] = tr(op[2 * i - 1:2 * i, 2 * i - 1:2 * i]) |> real
    end
    return CoordinateRepr(markers, :c)
end

# Arrow visualization

function _arrow_data(sz::NTuple{2,Integer}, cur::Real, i::Integer, j::Integer)
    if cur < 0
        i, j = j, i
        cur = -cur
    end
    local site = index_to_pair(sz, i)
    local other = index_to_pair(sz, j)
    local vec = normalize(other - site) * cur
    return Tuple(site), Tuple(vec)
end

function quiver_data(currents_mat::AbstractMatrix{<:Real}, sz::N{NTuple{2,Integer}}=nothing; threshold::Real=0.1, dist_threshold::Real=Inf)
    sz = _try_get_sz(sz)
    ps = Vector{NTuple{2,<:Integer}}()
    qs = Vector{NTuple{2,<:Real}}()
    for i in 1:prod(sz), j in 1:(i - 1)
        p, q = _arrow_data(sz, currents_mat[i, j], i, j)
        if norm(q) > threshold && dist(sz, i, j) < dist_threshold
            push!(ps, p)
            push!(qs, q)
        end 
    end
    return ps, qs
end

function quiver_currents!(pl::Plot, currents_mat::AbstractMatrix{<:Real}, sz::N{NTuple{2,Integer}}=nothing;
    threshold::Real=0.1, dist_threshold::Real=Inf, scale::Real=0, kw...)::Plots.Plot
    sz = _try_get_sz(sz)
    ps, qs = quiver_data(currents_mat, sz; threshold=threshold, dist_threshold=dist_threshold)
    if scale != 0
        qs .|> (arrow -> @. arrow * scale)
    else
        mx = qs .|> norm |> maximum
        preferred_length = max(1, √prod(sz) / 10)
        qs .|> (arrow -> @. arrow * preferred_length / mx)
    end
    return quiver!(pl, ps, quiver=qs, kw...)
end

quiver_currents!(currents_mat::AbstractMatrix{<:Real}, sz::NTuple{2,Integer}; kw...) = 
    quiver_currents!(current(), currents_mat, sz; kw...)

# Boundary visualization

function plot_boundaries!(pl::Plot, zone_mapping::CoordinateRepr; kw...)::Plots.plot
    local sz = size(zone_mapping)
    for i in 1:sz[1] - 1, j in 1:sz[2] - 1
        if zone_mapping[i, j] != zone_mapping[i + 1, j]
            plot!(pl, [(i + 0.5, j - 0.5), (i + 0.5, j + 0.5)], lab="", kw...)
        end
        if zone_mapping[i, j] != zone_mapping[i, j + 1]
            plot!(pl, [(i - 0.5, j + 0.5), (i + 0.5, j + 0.5)], lab="", kw...)
        end
    end
    return pl
end

plot_boundaries!(zone_mapping::AbstractMatrix; kw...)::Plots.Plot =
    plot_boundaries!(current(), zone_mapping, kw...)

function _expand_arg(arg, sz)
    mat = nothing
    tit = ""
    cur = nothing
    mat_type = Union{AbstractMatrix,CoordinateRepr}
    obtain_mat(obj) = obj isa AbstractMatrix ? heatmap_data(obj, sz) : obj
    
    if arg isa mat_type
        mat = obtain_mat(arg)
    elseif arg isa Pair
        tit = arg.first
        if arg.second isa mat_type
            mat = obtain_mat(arg.second)
        elseif arg.second isa Pair
            mat = obtain_mat(arg.second.first)
            cur = arg.second.second
        else
            error("could not expand arg $arg:\nInvalid operand#2 of pair, must be Pair or CoordinateRepr, not $(typeof(arg.second))")
        end
    else
        error("could not expand arg $arg:\nMust be Pair or CoordinateRepr, not $(typeof(arg))")
    end
    return mat, tit, cur, NamedTuple()
end

# Optimal layout

function _optimal_grid_size(plots_total, plot_aspect_ratio, cell_aspect_ratio)
    rows = round(Int, √(plots_total * plot_aspect_ratio[2] / plot_aspect_ratio[1] * cell_aspect_ratio[1] / cell_aspect_ratio[2]))
    rows = min(max(rows, 1), plots_total)
    cols = ceil(Int64, plots_total / rows)
    return rows, cols
end

function _optimal_size(plots_total; cell_size::NTuple{2,Integer}=(450, 350))
    rows, cols = _optimal_grid_size(plots_total, (1, 1), cell_size)
    return cell_size[1] * cols, cell_size[2] * rows
end


function optimal_layout(plots_total::Integer; plot_size::N{NTuple{2,Integer}}=nothing, cell_size::NTuple{2,Integer}=(450, 350))
    rows, cols = _optimal_grid_size(plots_total, plot_size === nothing ? (1, 1) : plot_size, cell_size)

    col_remainder = 1 - (rows * cols - plots_total) / cols

    if rows * cols == plots_total
        return grid(rows, cols)
    else
        return @layout [grid(rows - 1, cols);
        _ grid(1, plots_total - cols * (rows - 1), width=col_remainder * pct, height=1 / rows * pct) _]
    end
end

function _keys_by_prefix(dct, prefix::AbstractString)
    out = Dict()
    pl = length(prefix)
    for k in keys(dct)
        sk = String(k)
        if startswith(sk, prefix) && sk != prefix
            out[Symbol(sk[pl + 1:end])] = dct[k]
        end
    end
    return out
end

function _marker!(p1, p2, site; kw...)
    hline!(p1, [site[2]]; lab=nothing, kw...)
    plot!(p1, [site]; st=:scatter, lab=nothing, kw...)
    vline!(p2, [site[1]]; lab=nothing, kw...)
end

function plot_marker!(p, arg; lattice_size=nothing, clims=:auto, kw...)
    lattice_size = _try_get_sz(lattice_size)
    
    hmap_kw = _keys_by_prefix(kw, "hmap")
    currs_kw = _keys_by_prefix(kw, "currents")
    repr, tit, cur, p_args = _expand_arg(arg, lattice_size)
    heatmap!(p, repr; title=tit, clims=clims, cbar=:right, hmap_kw...)
    if cur !== nothing
        ps, qs = quiver_data(cur, lattice_size)
        quiver!(p, ps; quiver=qs, currs_kw...)
    end
    xlims!(p, (0, lattice_size[1] + 1))
    ylims!(p, (0, lattice_size[2] + 1))
end

"""
Plots multiple heatmaps, slices or currents simultaneously.

The subplots are automatically arranged into an optimal layout.

# Arguments

Each argument can be either a `CoordinateRepr` object or a chain of pairs.
"""
function plot_auto(args...; layout=nothing, plot_size=nothing, 
    zone_mapping::N{CoordinateRepr}=nothing, title="", control_site=nothing, clims=:auto, lattice_size=nothing, kw...)
    plots_total = length(args) + (control_site !== nothing)
    
    # Generate plot
    if plot_size === nothing
        plot_size = _optimal_size(plots_total)
        if layout !== nothing
            @warn "Plot size was not specified, fallback to optimal plot layout. Specify plot size to avoid this"
        end
        layout = optimal_layout(plots_total)
        println("Autogenerated size")
        
    elseif layout === nothing
        layout = optimal_layout(plots_total; plot_size=plot_size)
    end
    println(layout.grid)
    p = plot(layout=layout, size=plot_size)

    lattice_size = _try_get_sz(lattice_size)

    marker_kw = _keys_by_prefix(kw, "marker")
    bounds_kw = _keys_by_prefix(kw, "bounds")

    # Process args
    for i in 1:length(args)
        plot_marker!(p[i], args[1]; lattice_size=lattice_size, clims=clims, kw...)
        if zone_mapping !== nothing
            plot_boundaries!(p, zone_mapping, color=:black; bounds_kw...)
        end
        if control_site !== nothing
            _marker!(p[i], p[plots_total], control_site; marker_kw...)

            plot!(p[plots_total], repr[:, control_site[2]], lab=tit, ylim=(clims isa NTuple{2,<:Real} ? clims .* 2 : clims))
            hline!(p[plots_total], [-1, 1]; style=:dot, lab=nothing)
        end
        xlims!(p[i], (0, lattice_size[1] + 1))
        ylims!(p[i], (0, lattice_size[2] + 1))
    end
return plot!(p, plot_title=title)
end