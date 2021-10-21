using Plots: plot, plot!, Subplot, Plot, heatmap!, quiver!, xlims!, ylims!, hline!, vline!, current, @layout, grid, EmptyLayout, pct

AbstractPlot = Union{Subplot,Plot}

# Boundary visualization

function plot_boundaries!(pl::AbstractPlot, zone_mapping::CoordinateRepr; color=:black, kw...)::Plots.plot
    local lattice_size = size(zone_mapping)
    for i in 1:lattice_size[1] - 1, j in 1:lattice_size[2] - 1
        if zone_mapping[i, j] != zone_mapping[i + 1, j]
            plot!(pl, [(i + 0.5, j - 0.5), (i + 0.5, j + 0.5)]; lab="", color=color, kw...)
        end
        if zone_mapping[i, j] != zone_mapping[i, j + 1]
            plot!(pl, [(i - 0.5, j + 0.5), (i + 0.5, j + 0.5)]; lab="", color=color, kw...)
        end
    end
    return pl
end

plot_boundaries!(zone_mapping::AbstractMatrix; kw...)::Plots.Plot =
    plot_boundaries!(current(), zone_mapping; kw...)

_unchain_arg(arg) = 
    arg isa Pair ? Any[_unchain_arg(arg.first)..., _unchain_arg(arg.second)...] : Any[arg]

function _expand_arg(arg, lattice_size)
    mat_type = Union{AbstractMatrix,CoordinateRepr}
    _obtain_repr(obj) = obj isa CoordinateRepr ? obj : heatmap_data(obj, lattice_size)
    mat::N{CoordinateRepr} = nothing
    tit::AbstractString = ""
    cur::N{AbstractMatrix} = nothing

    arg_list = _unchain_arg(arg)
    if arg_list[1] isa AbstractString
        tit = popfirst!(arg_list)
    end
    if isempty(arg_list)
        error("Missing plot data in argument")
    end
    if length(arg_list) > 2
        error("Too long argument")
    end
    mat = _obtain_repr(arg_list[1])
    if length(arg_list) > 1
        cur = arg_list[2]
    end
    return mat, tit, cur
end

# Optimal layout

function _optimal_grid_size(plots_total::Integer, plot_aspect_ratio::NTuple{2,Integer}, subplot_aspect_ratio::NTuple{2,Integer})
    rows = round(Int, √(plots_total * plot_aspect_ratio[2] / plot_aspect_ratio[1] * subplot_aspect_ratio[1] / subplot_aspect_ratio[2]))
    rows = min(max(rows, 1), plots_total)
    cols = ceil(Int64, plots_total / rows)
return rows, cols
end

function _optimal_size(plots_total::Integer; subplot_size::NTuple{2,Integer}=(450, 350))
    rows, cols = _optimal_grid_size(plots_total, subplot_size, subplot_size)
    return subplot_size[1] * cols, subplot_size[2] * rows
end


function optimal_layout(plots_total::Integer; plot_size::N{NTuple{2,Integer}}=nothing, subplot_size::NTuple{2,Integer}=(450, 350))
    rows, cols = _optimal_grid_size(plots_total, plot_size === nothing ? subplot_size : plot_size, subplot_size)
    col_remainder = 1 - (rows * cols - plots_total) / cols

        if rows * cols == plots_total
        return grid(rows, cols)
    else
        return @layout [grid(rows - 1, cols, height=(rows - 1) / rows * pct);
        _ grid(1, plots_total - cols * (rows - 1), width=col_remainder * pct) _]
    end
end

function _keys_by_prefix(dct::Iterators.Pairs, prefix::AbstractString)
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

function plot_marker!(pl::AbstractPlot; hmap=nothing, currents=nothing, zone_mapping=nothing, xlims=nothing, ylims=nothing,
    lattice_size::N{NTuple{2,Integer}}=nothing, kw...)
    lattice_size = _try_get_lattice_size(lattice_size)
    hmap_kw = _keys_by_prefix(kw, "hmap")
    currs_kw = _keys_by_prefix(kw, "currents")
    bounds_kw = _keys_by_prefix(kw, "bounds")
    xlims = xlims !== nothing ? xlims : (0, lattice_size[1] + 1)
    ylims = ylims !== nothing ? ylims : (0, lattice_size[2] + 1)
    if hmap != nothing
        obtain_repr(obj) = obj isa CoordinateRepr ? obj : heatmap_data(obj, lattice_size)
        heatmap!(pl, obtain_repr(hmap); hmap_kw...)
    end
    if zone_mapping !== nothing
        plot_boundaries!(pl, zone_mapping; bounds_kw...)
    end
    if currents !== nothing
        ps, qs = quiver_data(currents, lattice_size, xlims=xlims, ylims=ylims)
        quiver!(pl, ps; quiver=qs, currs_kw...)
    end

    xlims!(pl, xlims)
    ylims!(pl, ylims)
end

"""
Plots multiple heatmaps, slices or currents simultaneously.

The subplots are automatically arranged into an optimal layout.

# Arguments

Each argument can be either a `CoordinateRepr` object or a chain of pairs.
"""
function plot_auto(args...; layout=nothing, plot_size=nothing, zone_mapping::N{CoordinateRepr}=nothing, title="",
     control_site=nothing, control_sites::AbstractVector{NTuple{2,Integer}}=Vector{NTuple{2,Integer}}(), lattice_size=nothing, kw...)
    lattice_size = _try_get_lattice_size(lattice_size)
    sites = []
    if control_site !== nothing
        push!(sites, control_site)
    end
    append!(sites, control_sites)
    plots_total = length(args) + length(sites)
    
    # Generate plot
    if plot_size === nothing
        plot_size = _optimal_size(plots_total)
        if layout !== nothing
            @warn "Plot size was not specified, fallback to optimal plot layout. Specify plot size to avoid this"
        end
        layout = optimal_layout(plots_total)
    elseif layout === nothing
        layout = optimal_layout(plots_total; plot_size=plot_size)
    end
    p = plot(layout=layout, size=plot_size)

    marker_kw = _keys_by_prefix(kw, "marker")
    slice_kw = _keys_by_prefix(kw, "slice")

    # Process args
    for i in 1:length(args)
        repr, tit, cur = _expand_arg(args[i], lattice_size)
        plot!(p[i], title=tit)
        plot_marker!(p[i]; hmap=repr, currents=cur, lattice_size=lattice_size, zone_mapping=zone_mapping, kw...)
        for site_idx = 1:length(sites)
            site = sites[site_idx]
            p_slice = p[site_idx + length(args)];
            plot!(p_slice, repr[:, site[2]]; title="@ $site", lab=tit, slice_kw...)
            vline!(p_slice, )
            hline!(p[i], [site[2]]; lab=nothing, marker_kw...)
            plot!(p[i], [site]; st=:scatter, lab=nothing, marker_kw...)
        end
        xlims!(p[i], (0, lattice_size[1] + 1))
        ylims!(p[i], (0, lattice_size[2] + 1))
    end
    return plot!(p, plot_title=title)
end