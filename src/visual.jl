using Plots: plot, plot!, Plot, heatmap!, current, @layout
using LinearAlgebra: norm, normalize
using Base:getindex

function heatmap_data(op::AbstractMatrix{<:Complex{T}}, sz::NTuple{2,Integer})::CoordinateRepr{T} where T <: Real
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
    return site, vec |> Tuple
end

function quiver_data(sz::NTuple{2,Integer}, currents_mat::AbstractMatrix{Real}; threshold::Real=0.1, dist_threshold::Real=Inf)
    ps = []
    qs = []
    for i in 1:prod(sz), j in 1:(i - 1)
        p, q = _arrow_data(sz, currents_mat[i, j], i, j)
        if norm(q) > threshold && dist(sz, i, j) < dist_threshold
            push!(ps, p)
            push!(qs, q)
        end 
    end
    return ps, qs
end

function quiver_currents!(pl::Plot, sz::NTuple{2,Integer}, currents_mat::AbstractMatrix{Real};
    threshold::Real=0.1, dist_threshold::Real=Inf, scale::Real=1, kw...)::Plots.Plot
    ps, qs = quiver_data(sz, currents_mat, threshold=threshold, dist_threshold=dist_threshold)
    if scale != 1
        qs .|> (arrow -> @. arrow * scale)
    else
        mx = qs .|> norm |> maximum
        preferred_length = max(1, √prod(sz) / 10)
        qs .|> (arrow -> @. arrow * preferred_length / mx)
    end
    return quiver!(pl, ps, quiver=qs, kw...)
end

quiver_currents!(sz::NTuple{2,Integer}, currents_mat::AbstractMatrix{Real}; kw...) =
    quiver_currents!(current(), sz, currents_mat; kw...)

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

function _expand_arg(arg)
    mat = nothing
    tit = ""
    cur = nothing
    if arg isa CoordinateRepr
        mat = arg
    elseif arg isa Pair
        tit = arg.first
        if arg.second isa CoordinateRepr
            mat = arg.second
        elseif arg.second isa Pair
            mat = arg.second.first
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
function plot_arranged(args...; zone_mapping=nothing, title="", control_site=nothing, plot_size=nothing, cell_size=(400, 350), title_margin=50, clims=:auto, legend=:best)
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

    # Process args
    for i in 1:length(args)
        repr, tit, cur = _expand_arg(args[i])
        heatmap!(p[i], repr, title=tit, clims=clims, cbar=:right)
        if cur !== nothing
            xs, ys, qs = quiver_data(siz, cur)
            quiver!(p[i], xs, ys; quiver=qs, color="brown",)
        end
        if zone_mapping !== nothing
            plot_boundaries!(p[i], zone_mapping, color=:black)
        end
        if control_site !== nothing
            hline!(p[i], [control_site[2]]; color="brown")
            plot!(p[i], [control_site]; st=:scatter, color="brown")

            plot!(p[plots_total], ylim=(clims isa NTuple{2,<:Real} ? clims .* 2 : clims))
            plot!(p[plots_total], mat[:, control_site[2]], lab=tit, kw...)
            hline!(p[plots_total], [-1, 1]; style=:dot, lab=nothing)
            vline!(p[plots_total], [control_site[1]]; color="brown", lab=nothing)
        end
        xlims!(p[i], (0, size(mat)[1] + 1))
        ylims!(p[i], (0, size(mat)[2] + 1))
    end
    return plot!(p, plot_title=title)
end