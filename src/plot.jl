using Plots: plot, plot!, Subplot, Plot, heatmap!, quiver!, xlims!, ylims!, hline!, vline!, current, @layout, grid, EmptyLayout, pct

AbstractPlot = Union{Subplot,Plot}

# Boundary visualization

"""
    plot_boundaries!([pl, ]domain_mapping; <keyword arguments>)
Draws boundaries between different domains described by the `domain_mapping` matrix.

# Arguments
- `pl`: a `Plots.Plot` object to visualize data on
- `domain_mapping`: a `CoordinateRepr` object that represents domain mapping. The boundaries between different domains will be drawn

All keyword arguments will be passed to the `plot!` function used for drawing - this can be used to change the line thickness or style, for example.
"""
function plot_boundaries!(pl::AbstractPlot, domain_mapping::CoordinateRepr; color=:black, kw...)
    local lattice_size = size(domain_mapping)
    for i in 1:lattice_size[1] - 1, j in 1:lattice_size[2] - 1
        if domain_mapping[i, j] != domain_mapping[i + 1, j]
            plot!(pl, [(i + 0.5, j - 0.5), (i + 0.5, j + 0.5)]; lab="", color=color, kw...)
        end
        if domain_mapping[i, j] != domain_mapping[i, j + 1]
            plot!(pl, [(i - 0.5, j + 0.5), (i + 0.5, j + 0.5)]; lab="", color=color, kw...)
        end
    end
    return pl
end

plot_boundaries!(domain_mapping::AbstractMatrix; kw...) =
    plot_boundaries!(current(), domain_mapping; kw...)

_unchain_arg(arg) = 
    arg isa Pair ? Any[_unchain_arg(arg.first)..., _unchain_arg(arg.second)...] : Any[arg]

_obtain_repr(obj, lattice_size::SizeType)::CoordinateRepr = obj isa CoordinateRepr ? obj : heatmap_data(obj, lattice_size)

function _expand_arg(arg, lattice_size)
    mat_type = Union{AbstractMatrix,CoordinateRepr}
    mat::Nullable{CoordinateRepr} = nothing
    tit::AbstractString = ""
    cur::Nullable{AbstractMatrix} = nothing

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
    mat = _obtain_repr(arg_list[1], lattice_size)
    if length(arg_list) > 1
        cur = arg_list[2]
    end
    return mat, tit, cur
end

# Optimal layout

function _optimal_grid_size(plots_total::Int, plot_aspect_ratio::NTuple{2,Int}, subplot_aspect_ratio::NTuple{2,Int})::NTuple{2, Int}
    rows = round(Int, √(plots_total * plot_aspect_ratio[2] / plot_aspect_ratio[1] * subplot_aspect_ratio[1] / subplot_aspect_ratio[2]))
    rows = min(max(rows, 1), plots_total)
    cols = ceil(Int64, plots_total / rows)
return rows, cols
end

function _optimal_size(plots_total::Int; subplot_size::NTuple{2,Int}=(450, 350))::NTuple{2, Int}
    rows, cols = _optimal_grid_size(plots_total, subplot_size, subplot_size)
    return subplot_size[1] * cols, subplot_size[2] * rows
end


function optimal_layout(plots_total::Int; plot_size::SizeType=nothing, subplot_size::NTuple{2,Int}=(450, 350))
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

"""
    plot_figure!(pl; <keyword arguments>)

Plots complicated splitline data series (heatmap, boundaries, quiver) on a single figure.

# Arguments
- `pl`: a `Plots.Plot` object to visualize data on
- `hmap`: data to be visualized on a heatmap. It can be a `CoordinateRepr` object (then it will be plotted directly) 
or a linear operator matrix (then the `CoordinateRepr` will be generated automatically)
- `domain_mapping`: a `CoordinateRepr` object that represents domain mapping. The boundaries between different domains will be drawn.
- `currents`: a matrix containing currents between sites
- `xlims` and `ylims`: objects of type `Tuple{Int, Int}` that define the limits of the x- and y- axes respectively

All keyword arguments with different prefixes are passed to the `plot!` function:
- `hmap` for the heatmap
- `bounds` for the boundaries
- `currents` for the quiver

This can be used to style the plot:

`plot_figure!(..., hmapclims=(-3, 3), boundsstyle=:dot, :currentscolor=:green)`
"""
function plot_figure!(pl::AbstractPlot; hmap=nothing, currents=nothing, domain_mapping=nothing, xlims::SizeType=nothing, ylims::SizeType=nothing,
    lattice_size::SizeType=nothing, kw...)
    lattice_size = _try_get_lattice_size(lattice_size)
    hmap_kw = _keys_by_prefix(kw, "hmap")
    currs_kw = _keys_by_prefix(kw, "currents")
    bounds_kw = _keys_by_prefix(kw, "bounds")
    xlims = xlims !== nothing ? xlims : (0, lattice_size[1] + 1)
    ylims = ylims !== nothing ? ylims : (0, lattice_size[2] + 1)
    if hmap !== nothing
        heatmap!(pl, _obtain_repr(hmap, lattice_size); hmap_kw...)
    end
    if domain_mapping !== nothing
        plot_boundaries!(pl, domain_mapping; bounds_kw...)
    end
    if currents !== nothing
        ps, qs = quiver_data(currents, lattice_size, xlims=xlims, ylims=ylims)
        quiver!(pl, ps; quiver=qs, currs_kw...)
    end

    xlims!(pl, xlims)
    ylims!(pl, ylims)
end

"""
    plot_auto(<arguments>; <keyword arguments>)
Plots multiple heatmaps, cutaway views or currents simultaneously.

The subplots are automatically arranged into an optimal layout.

# Arguments

Each argument can be either a `CoordinateRepr` object or a chain of pairs.
"""
function plot_auto(args...; layout=nothing, plot_size=nothing, domain_mapping::Nullable{CoordinateRepr}=nothing, title="",
     cutaway_view=nothing, cutaway_views::AbstractVector{NTuple{2,Int}}=Vector{NTuple{2,Int}}(), lattice_size=nothing, kw...)
    lattice_size = _try_get_lattice_size(lattice_size)
    sites = []
    if cutaway_view !== nothing
        push!(sites, cutaway_view)
    end
    append!(sites, cutaway_views)
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

    splitline_kw = _keys_by_prefix(kw, "splitline")
    cutaway_kw = _keys_by_prefix(kw, "cutaway")

    # Process args
    for i in 1:length(args)
        repr, tit, cur = _expand_arg(args[i], lattice_size)
        plot!(p[i], title=tit)
        plot_figure!(p[i]; hmap=repr, currents=cur, lattice_size=lattice_size, domain_mapping=domain_mapping, kw...)
        for site_idx = 1:length(sites)
            site = sites[site_idx]
            p_cutaway = p[site_idx + length(args)];
            plot!(p_cutaway, repr[:, site[2]]; title="@ $site", lab=tit, cutaway_kw...)
            vline!(p_cutaway, )
            hline!(p[i], [site[2]]; lab=nothing, splitline_kw...)
            plot!(p[i], [site]; st=:scatter, lab=nothing, splitline_kw...)
        end
        xlims!(p[i], (0, lattice_size[1] + 1))
        ylims!(p[i], (0, lattice_size[2] + 1))
    end
    return plot!(p, plot_title=title)
end