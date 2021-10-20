using Plots: plot, plot!, Subplot, Plot, heatmap!, quiver!, xlims!, ylims!, hline!, vline!, current, @layout, grid, EmptyLayout, pct

AbstractPlot = Union{Subplot,Plot}

function quiver_currents!(pl::AbstractPlot, currents_mat::AbstractMatrix{<:Real}, sz::N{NTuple{2,Integer}}=nothing;
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

function plot_boundaries!(pl::AbstractPlot, zone_mapping::CoordinateRepr; kw...)::Plots.plot
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

_unchain_arg(arg)::Vector{Any} = 
    arg isa Pair ? [_unchain_arg(arg.first)..., _unchain_arg(arg.second)...] : [arg]

function _expand_arg(arg, sz)
    mat_type = Union{AbstractMatrix,CoordinateRepr}
    obtain_mat(obj) = obj isa AbstractMatrix ? heatmap_data(obj, sz) : obj
    mat::N{mat_type} = nothing
    tit::AbstractString = ""
    cur::N{AbstractMatrix} = nothing

    arg_list = _unchain_arg(arg)
    if length(arg_list) == 1
        mat = obtain_mat(arg_list[1])
    else
        tit = arg_list[1]
        mat = obtain_mat(arg_list[2])
    end
    if length(arg_list) ≥ 3 && arg_list[3] isa AbstractMatrix
        cur = arg_list[3]
    end
    
    return mat, tit, cur
end

function _obtain_attr(arg)::NamedTuple
    if _unchain_arg(arg)[end] isa NamedTuple
        return _unchain_arg(arg)[end]
    else
        return NamedTuple()
    end
end

# Optimal layout

function _optimal_grid_size(plots_total::Integer, plot_aspect_ratio::NTuple{2, Integer}, cell_aspect_ratio::NTuple{2, Integer})
    rows = round(Int, √(plots_total * plot_aspect_ratio[2] / plot_aspect_ratio[1] * cell_aspect_ratio[1] / cell_aspect_ratio[2]))
    rows = min(max(rows, 1), plots_total)
    cols = ceil(Int64, plots_total / rows)
return rows, cols
end

function _optimal_size(plots_total::Integer; cell_size::NTuple{2,Integer}=(450, 350))
    rows, cols = _optimal_grid_size(plots_total, cell_size, cell_size)
    return cell_size[1] * cols, cell_size[2] * rows
end


function optimal_layout(plots_total::Integer; plot_size::N{NTuple{2,Integer}}=nothing, cell_size::NTuple{2,Integer}=(450, 350))
    rows, cols = _optimal_grid_size(plots_total, plot_size === nothing ? cell_size : plot_size, cell_size)

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

function plot_marker!(p::AbstractPlot, arg; lattice_size::N{NTuple{2,Integer}}=nothing, clims=:auto, kw...)
    lattice_size = _try_get_sz(lattice_size)
    hmap_kw = _keys_by_prefix(kw, "hmap")
    currs_kw = _keys_by_prefix(kw, "currents")
    repr, tit, cur = _expand_arg(arg, lattice_size)
    heatmap!(p, repr; title=tit, clims=clims, cbar=:right, hmap_kw...)
    if cur !== nothing
        ps, qs = quiver_data(cur, lattice_size)
        quiver!(p, ps; quiver=qs, currs_kw...)
    end

    xlims!(p, (0, lattice_size[1] + 1))
    ylims!(p, (0, lattice_size[2] + 1))
end

function _arg_subplot_mapping(args, additional_sites)
    arg_mapping = Dict()
    sites = []
    function index_of_site(site)
        if site ∉ sites
            push!(sites, site)
        end
        return indexin([site], sites)[1]
    end
    for arg in args
        attr = _obtain_attr(arg)
        arg_mapping[arg] = Vector{Int}()
        for s in additional_sites
            push!(arg_mapping[arg], index_of_site(s))
        end
        if :site ∈ keys(attr)
            push!(arg_mapping[arg], index_of_site(attr.site))
        elseif :sites ∈ keys(attr)
            for s in attr.sites
                push!(arg_mapping[arg], index_of_site(s))
            end
        end
    end
    return arg_mapping, sites
end

"""
Plots multiple heatmaps, slices or currents simultaneously.

The subplots are automatically arranged into an optimal layout.

# Arguments

Each argument can be either a `CoordinateRepr` object or a chain of pairs.
"""
function plot_auto(args...; layout=nothing, plot_size=nothing, zone_mapping::N{CoordinateRepr}=nothing, title="",
     control_site=nothing, control_sites::AbstractVector{NTuple{2,Integer}}=Vector{NTuple{2,Integer}}(), clims=:auto, lattice_size=nothing, kw...)
    lattice_size = _try_get_sz(lattice_size)
    if control_site !== nothing
        push!(control_sites, control_site)
    end
    arg_mapping, sites = _arg_subplot_mapping(args, Set(control_sites))

    plots_total = length(args) + (control_site !== nothing) + length(sites)
    
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

    bounds_kw = _keys_by_prefix(kw, "bounds")
    marker_kw = _keys_by_prefix(kw, "marker")

    # Process args
    for i in 1:length(args)
        if zone_mapping !== nothing
            plot_boundaries!(p, zone_mapping, color=:black; bounds_kw...)
        end
        plot_marker!(p[i], args[i]; lattice_size=lattice_size, clims=clims, kw...)
        repr, tit = _expand_arg(args[i], lattice_size)[1:2]
        for site_idx in arg_mapping[args[i]]
            site = sites[site_idx]
            p_slice = p[site_idx + length(args)];
            plot!(p_slice, repr[:, site[2]], title="@ $site", lab=tit, ylim=(clims isa NTuple{2,<:Real} ? clims .* 2 : clims))
            vline!(p_slice, )
            hline!(p[i], [site[2]]; lab=nothing, marker_kw...)
            plot!(p[i], [site]; st=:scatter, lab=nothing, marker_kw...)
        end
        xlims!(p[i], (0, lattice_size[1] + 1))
        ylims!(p[i], (0, lattice_size[2] + 1))
    end
    return plot!(p, plot_title=title)
end