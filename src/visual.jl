using Plots: plot, plot!, Plot, heatmap!, quiver!, xlims!, ylims!, 
    hline!, vline!, current, @layout, grid
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
    return mat, tit, cur
end

"""
Plots multiple heatmaps, slices or currents simultaneously.

The subplots are automatically arranged according to plot size (if provided) or into a square grid (by default).

# Arguments

Each argument can be either a `CoordinateRepr` object or a chain of pairs.
"""
function plot_arranged(args...; zone_mapping=nothing, title="", control_site=nothing, plot_size=nothing, cell_size=(450, 350), title_margin=50, clims=:auto, legend=:best, lattice_size=nothing)
    # Obtain optimal size & layout
    plots_total = length(args) + (control_site !== nothing)
    if plot_size === nothing
        rows = round(Int, √(plots_total * cell_size[1] / cell_size[2]))
    else
        rows = round(Int, √(plots_total * (plot_size[2] - title_margin) / plot_size[1] * cell_size[1] / cell_size[2]))
    end
    rows = min(max(rows, 1), plots_total)
    cols = ceil(Int64, plots_total / rows)
    if plot_size === nothing
        plot_size = (cell_size[1] * cols, title_margin + cell_size[2] * rows)
    end
    
    # Generate plot
    if rows == one(rows)
        l = grid(1, cols)
    else
        l = @layout [grid(rows - 1, cols); grid(1, plots_total - cols * (rows - 1))]
    end
    p = plot(layout=l, legend=legend, size=plot_size)

    lattice_size = _try_get_sz(lattice_size)

    # Process args
    for i in 1:length(args)
        repr, tit, cur = _expand_arg(args[i], lattice_size)
        heatmap!(p[i], repr, title=tit, clims=clims, cbar=:right)
        if cur !== nothing
            ps, qs = quiver_data(cur, lattice_size)
            quiver!(p[i], ps; quiver=qs, color="brown",)
        end
        if zone_mapping !== nothing
            plot_boundaries!(p[i], zone_mapping, color=:black)
        end
        if control_site !== nothing
            hline!(p[i], [control_site[2]]; color="brown", lab=nothing)
            plot!(p[i], [control_site]; st=:scatter, color="brown", lab=nothing)

            plot!(p[plots_total], ylim=(clims isa NTuple{2,<:Real} ? clims .* 2 : clims))
            plot!(p[plots_total], repr[:, control_site[2]], lab=tit)
            hline!(p[plots_total], [-1, 1]; style=:dot, lab=nothing)
            vline!(p[plots_total], [control_site[1]]; color="brown", lab=nothing)
        end
        xlims!(p[i], (0, lattice_size[1] + 1))
        ylims!(p[i], (0, lattice_size[2] + 1))
    end
    return plot!(p, plot_title=title)
end