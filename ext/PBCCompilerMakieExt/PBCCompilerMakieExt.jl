module PBCCompilerMakieExt

using Makie
using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, affectedqubits
using Graphs, SimpleWeightedGraphs, KaHyPar, SparseArrays, LinearAlgebra
using Moshi.Match: @match

import PBCCompiler: circuitplot, circuitplot!, circuitplot_axis, plot_histogram, plot_graph, plot_weight_histogram, plot_std_graph, plot_partition, plot_hypergraph, plot_hypergraph_partition, plot_hyperedge_frequency

# Define the recipe with attributes
Makie.@recipe(CircuitPlot, circuit) do scene
    Makie.Theme(;
        # Gate dimensions
        gatewidth = 0.8,
        qubitspacing = 1.0,
        # Wire appearance
        wirecolor = :black,
        wirelinewidth = 1.0,
        # Gate colors by variant
        paulicolor = Makie.RGB(0.2, 0.6, 0.2),           # green
        measurementcolor = Makie.RGB(0.8, 0.2, 0.2),     # red
        exphalfpicolor = Makie.RGB(0.2, 0.4, 0.8),       # blue
        expquatpicolor = Makie.RGB(0.5, 0.2, 0.8),       # purple
        expeighpicolor = Makie.RGB(0.8, 0.4, 0.8),       # magenta
        conditionalcolor = Makie.RGB(0.8, 0.6, 0.2),     # orange
        bitconditionalcolor = Makie.RGB(0.6, 0.6, 0.6),  # gray
        # Text appearance
        fontsize = 0.5,
        textcolor = :white,
    )
end

"""Get the color for a CircuitOp variant."""
function gate_color(op::CircuitOp.Type, plot::CircuitPlot)
    @match op begin
        CircuitOp.Measurement(_, _, _) => plot.measurementcolor[]
        CircuitOp.Pauli(_, _) => plot.paulicolor[]
        CircuitOp.ExpHalfPiPauli(_, _) => plot.exphalfpicolor[]
        CircuitOp.ExpQuatPiPauli(_, _) => plot.expquatpicolor[]
        CircuitOp.ExpEighPiPauli(_, _) => plot.expeighpicolor[]
        CircuitOp.PrepMagic(_, _) => :transparent  # Not visualized
        CircuitOp.PauliConditional(_, _, _, _) => plot.conditionalcolor[]
        CircuitOp.BitConditional(_, _) => plot.bitconditionalcolor[]
    end
end

"""Get a short label for a CircuitOp variant."""
function gate_label(op::CircuitOp.Type)
    @match op begin
        CircuitOp.Measurement(_, _, _) => "M"
        CircuitOp.Pauli(_, _) => "P"
        CircuitOp.ExpHalfPiPauli(_, _) => "S"
        CircuitOp.ExpQuatPiPauli(_, _) => "T"
        CircuitOp.ExpEighPiPauli(_, _) => "R"
        CircuitOp.PrepMagic(_, _) => ""
        CircuitOp.PauliConditional(_, _, _, _) => "CP"
        CircuitOp.BitConditional(inner, _) => gate_label(inner)
    end
end

"""Check if an operation is a PrepMagic (not visualized for now)."""
function is_prepmagic(op::CircuitOp.Type)
    @match op begin
        CircuitOp.PrepMagic(_, _) => true
        _ => false
    end
end

"""Get the classical bit index for a Measurement, or nothing."""
function measurement_bit(op::CircuitOp.Type)
    @match op begin
        CircuitOp.Measurement(_, bit, _) => bit
        CircuitOp.BitConditional(inner, _) => measurement_bit(inner)
        _ => nothing
    end
end

"""Get the conditioning bit index for a BitConditional, or nothing."""
function conditioning_bit(op::CircuitOp.Type)
    @match op begin
        CircuitOp.BitConditional(_, bit) => bit
        _ => nothing
    end
end

function Makie.plot!(plot::CircuitPlot)
    circuit = plot[:circuit][]

    if isempty(circuit)
        return plot
    end

    # Get all qubits in the circuit
    all_qubits = affectedqubits(circuit)
    if isempty(all_qubits)
        return plot
    end

    min_qubit = minimum(all_qubits)
    max_qubit = maximum(all_qubits)

    gw = plot.gatewidth[]
    qs = plot.qubitspacing[]

    # Draw qubit wires (horizontal lines)
    for q in min_qubit:max_qubit
        y = q * qs
        # Wire extends from before first gate to after last gate
        Makie.lines!(plot, [0.5, length(circuit) + 0.5], [y, y];
            color = plot.wirecolor[],
            linewidth = plot.wirelinewidth[]
        )
    end

    # Draw each gate
    for (idx, op) in enumerate(circuit)
        # Skip PrepMagic for now
        if is_prepmagic(op)
            continue
        end

        qubits = affectedqubits(op)
        if isempty(qubits)
            continue
        end

        # Calculate rectangle bounds
        x_center = idx
        x_left = x_center - gw / 2
        x_right = x_center + gw / 2

        y_min = minimum(qubits) * qs
        y_max = maximum(qubits) * qs
        y_min -= 0.3 * qs
        y_max += 0.3 * qs

        # Draw rectangle
        color = gate_color(op, plot)
        Makie.poly!(plot,
            Makie.Point2f[(x_left, y_min), (x_right, y_min), (x_right, y_max), (x_left, y_max)];
            color = color,
            strokecolor = :black,
            strokewidth = 1
        )

        # Draw gate label in center
        label = gate_label(op)
        if !isempty(label)
            Makie.text!(plot, x_center, (y_min + y_max) / 2;
                text = label,
                align = (:center, :center),
                fontsize = plot.fontsize[],
                color = plot.textcolor[],
                markerspace = :data
            )
        end

        # Draw measurement bit index (top-right corner)
        mbit = measurement_bit(op)
        if mbit !== nothing
            Makie.text!(plot, x_right - 0.05, y_min + 0.1;
                text = "$mbit",
                align = (:right, :bottom),
                fontsize = plot.fontsize[] * 0.7,
                color = :black,
                markerspace = :data
            )
        end

        # Draw conditioning bit index (bottom-left corner)
        cbit = conditioning_bit(op)
        if cbit !== nothing
            Makie.text!(plot, x_left + 0.05, y_min + 0.1;
                text = "$cbit",
                align = (:left, :bottom),
                fontsize = plot.fontsize[] * 0.7,
                color = :black,
                markerspace = :data
            )
        end
    end

    return plot
end

"""
    circuitplot_axis(subfig, circuit; kwargs...)

Create a complete Makie figure panel with a circuit plot and appropriate axis settings.

Returns a tuple of (subfig, axis, plot).
"""
function circuitplot_axis(subfig, circuit::Circuit; kwargs...)
    ax = Makie.Axis(subfig[1, 1])
    p = circuitplot!(ax, circuit; kwargs...)

    # Configure axis
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    ax.aspect = Makie.DataAspect()

    # Set axis limits with padding
    if !isempty(circuit)
        all_qubits = affectedqubits(circuit)
        if !isempty(all_qubits)
            min_q = minimum(all_qubits)
            max_q = maximum(all_qubits)
            Makie.xlims!(ax, 0, length(circuit) + 1)
            Makie.ylims!(ax, min_q - 0.5, max_q + 0.5)
        end
    end

    return (subfig, ax, p)
end

"""
    plot_histogram(data)

Plot a histogram of integer-valued `data` using CairoMakie.

- Each bar is centered over its integer value, with the x-axis label appearing
  directly below the bar (not at bar edges).
- X-axis tick labels are shown as bitstrings.
- The frequency count is labeled on top of each bar.
"""
function plot_histogram(data)
    freq_dict = Dict{Int,Int}()
    for v in data
        freq_dict[v] = get(freq_dict, v, 0) + 1
    end
    values = sort(collect(keys(freq_dict)))
    freqs  = [freq_dict[v] for v in values]

    fig = Figure()
    ax  = Axis(fig[1, 1])

    barplot!(ax, values, freqs; width = 1.0, gap = 0.0)

    # X-axis: one tick per unique value, labeled as compact bitstring
    # Use only as many bits as needed to represent the largest value
    nbits = values[end] == 0 ? 1 : Int(ceil(log2(values[end] + 1)))
    ax.xticks = (values, [string(v, base=2, pad=nbits) for v in values])
    ax.xticklabelrotation = π / 2
    ax.title = "Measurement Result Distribution"
    ax.ylabel = "Count"
    ax.xlabel = "Bitstring"

    # Frequency labels on top of each bar
    for (v, f) in zip(values, freqs)
        text!(ax, v, f; text = string(f), align = (:center, :bottom))
    end

    return fig
end

"""
    plot_graph(weights::AbstractMatrix{<:Real}) -> Figure

Draw a weighted graph from its adjacency matrix.

Vertices are placed on a circle. Edges with positive weight are drawn in gray
(opacity and thickness scale with weight); edges with negative weight are drawn
in red (opacity and thickness scale with absolute weight). Zero-weight entries
are skipped.

# Arguments
- `weights`: symmetric weight matrix of size n × n

# Returns
A `CairoMakie.Figure`.
"""
function plot_graph(weights::AbstractMatrix{<:Real})
    n = size(weights, 1)

    fig = Figure(size=(600, 600))
    ax = Axis(fig[1,1], title="Interaction Graph",
              aspect=DataAspect(), limits=(-1.5, 1.5, -1.5, 1.5))
    hidedecorations!(ax)
    hidespines!(ax)

    if n == 0 || all(iszero, weights)
        return fig
    end

    angles = [2π * (i-1) / n - π/2 for i in 1:n]
    xs = cos.(angles)
    ys = sin.(angles)

    pos_max = let ws = [weights[i,j] for i in 1:n for j in i+1:n if weights[i,j] > 0]
        isempty(ws) ? 1.0 : Float64(maximum(ws))
    end
    neg_max = let ws = [abs(weights[i,j]) for i in 1:n for j in i+1:n if weights[i,j] < 0]
        isempty(ws) ? 1.0 : Float64(maximum(ws))
    end

    for i in 1:n, j in i+1:n
        w = weights[i, j]
        iszero(w) && continue
        if w > 0
            frac = w / pos_max
            edge_color = (:gray20, 0.3 + 0.7 * frac)
        else
            frac = abs(w) / neg_max
            edge_color = (:crimson, 0.3 + 0.7 * frac)
        end
        lines!(ax, [xs[i], xs[j]], [ys[i], ys[j]];
               color=edge_color, linewidth=1 + 3 * frac)
        mx, my = (xs[i] + xs[j]) / 2, (ys[i] + ys[j]) / 2
        text!(ax, mx, my; text=string(w), fontsize=11, align=(:center, :center))
    end

    scatter!(ax, xs, ys; markersize=35, color=:steelblue)
    for i in 1:n
        text!(ax, xs[i], ys[i]; text=string(i), fontsize=14,
              align=(:center, :center), color=:white)
    end

    return fig
end

"""
    plot_weight_histogram(g::SimpleWeightedGraph; bins=nothing)

Plot a histogram of edge weights from `g`, assuming integer edge weights.

When `bins` is `nothing` (default), each distinct integer weight gets its own bar.
When `bins` is an integer, weights are grouped into approximately that many bins
with integer-aligned boundaries. Returns a `Makie.Figure`.
"""
function plot_weight_histogram(g::SimpleWeightedGraph; bins::Union{Int,Nothing}=nothing)
    ws = [e.weight for e in edges(g)]
    lo = round(Int, minimum(ws))
    hi = round(Int, maximum(ws))

    if isnothing(bins)
        bin_edges = (lo - 1 : hi) .+ 0.5   # half-integer edges, one bin per integer
        midpoints = Float64.(lo:hi)
        labels    = string.(lo:hi)
    else
        step      = max(1, ceil(Int, (hi - lo) / bins))
        start     = (lo ÷ step) * step
        stop      = ((hi - 1) ÷ step + 1) * step
        bin_edges = Float64.(start:step:stop)
        midpoints = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in 1:length(bin_edges)-1]
        labels    = ["$(Int(bin_edges[i]))–$(Int(bin_edges[i+1]))"
                     for i in 1:length(bin_edges)-1]
    end

    rot = (isnothing(bins) && hi - lo <= 15) ? 0.0 : π/4

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Edge Weight", ylabel="Count",
              title="Edge Weight Distribution",
              xticks=(midpoints, labels),
              xticklabelrotation=rot)
    hist!(ax, ws; bins=bin_edges, color=:steelblue)
    return fig
end

##
"""
    plot_std_graph(g::SimpleWeightedGraph; colormap=Reverse(:RdBu), node_size=30)

Plot the output of `weight_std_graph` as a network diagram using a circular layout.

Nodes are labeled 1–n. Edge color encodes the edge weight (std dev) on a
blue-to-red scale, with a colorbar on the right. Returns a `Makie.Figure`.
"""
function plot_std_graph(g::SimpleWeightedGraph; colormap=Reverse(:RdBu), node_size=30)
    n = nv(g)
    angles = range(0, 2π; length=n + 1)[1:n]
    pos = [Point2f(cos(a), sin(a)) for a in angles]

    ws = [e.weight for e in edges(g)]
    wmin, wmax = extrema(ws)

    segments = Point2f[]
    seg_weights = Float64[]
    for e in edges(g)
        push!(segments, pos[src(e)], pos[dst(e)])
        push!(seg_weights, e.weight, e.weight)
    end

    fig = Figure(size=(700, 600))
    ax = Axis(fig[1, 1]; aspect=DataAspect(), title="Edge Weight Std Dev Graph")
    hidedecorations!(ax)
    hidespines!(ax)

    ls = linesegments!(ax, segments; color=seg_weights, colormap=colormap,
                       colorrange=(wmin, wmax), linewidth=6)
    Label(fig[0, 2], "volatile"; fontsize=11, tellwidth=false)
    Colorbar(fig[1, 2], ls; label="Std Dev of Edge Weight", width=25, spinewidth=0)
    Label(fig[2, 2], "stable"; fontsize=11, tellwidth=false)

    scatter!(ax, pos; markersize=node_size, color=:white,
             strokecolor=:black, strokewidth=1.5)
    text!(ax, pos; text=string.(1:n), align=(:center, :center), fontsize=14)

    return fig
end
##
const _PARTITION_PALETTE = [:tomato, :steelblue, :seagreen, :gold,
                             :mediumpurple, :darkorange, :deepskyblue]

"""
    plot_partition(g::SimpleWeightedGraph, part::Vector{Int32}) -> Figure

Display a circular-layout plot of `g` with nodes colored by partition assignment.

# Arguments
- `g`: the graph that was partitioned
- `part`: partition assignment vector as returned by `METIS_partition`

# Returns
A `CairoMakie.Figure` with the partition plot.
"""
function plot_partition(g::SimpleWeightedGraph, part::Vector{Int32})::Figure
    num_partition = maximum(part)
    n = nv(g)
    θ = range(0, 2π, length=n+1)[1:n]
    xs, ys = cos.(θ), sin.(θ)
    node_colors = [_PARTITION_PALETTE[mod1(p, length(_PARTITION_PALETTE))] for p in part]

    fig = Figure(size=(500, 500))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="METIS partition (k=$num_partition)")
    hidedecorations!(ax)
    hidespines!(ax)

    for e in edges(g)
        u, v = src(e), dst(e)
        lines!(ax, [xs[u], xs[v]], [ys[u], ys[v]], color=:gray70, linewidth=1)
    end

    scatter!(ax, xs, ys, color=node_colors, markersize=35,
             strokewidth=1.5, strokecolor=:black)

    for i in 1:n
        text!(ax, xs[i], ys[i], text=string(i),
              align=(:center, :center), fontsize=13, color=:white)
    end

    return fig
end
##
# Fruchterman-Reingold spring layout; returns (n_nodes × 2) position matrix.
function _spring_layout(n_nodes::Int, edges::Vector{Tuple{Int,Int}}; iterations::Int = 200)::Matrix{Float64}
    pos = rand(n_nodes, 2) .* 2.0 .- 1.0
    k = sqrt(4.0 / max(n_nodes, 1))
    for iter in 1:iterations
        t = 0.1 * (1.0 - iter / iterations)
        disp = zeros(n_nodes, 2)
        for i in 1:n_nodes
            for j in 1:n_nodes
                i == j && continue
                delta = pos[i, :] .- pos[j, :]
                d = max(norm(delta), 1e-6)
                disp[i, :] .+= delta .* (k^2 / d^2)
            end
        end
        for (u, v) in edges
            delta = pos[u, :] .- pos[v, :]
            d = max(norm(delta), 1e-6)
            disp[u, :] .-= delta .* (d / k)
            disp[v, :] .+= delta .* (d / k)
        end
        for i in 1:n_nodes
            d = max(norm(disp[i, :]), 1e-6)
            pos[i, :] .+= disp[i, :] .* (min(d, t) / d)
        end
    end
    return pos
end

# Two-level partition-aware layout: partition centers on an outer circle, vertices
# evenly arranged on a sub-circle within each partition, fake vertices at the
# centroid of their connected real vertices.
# Sub-circle radius scales as sqrt(n_vertices) so blob area ∝ partition size.
function _partition_layout(n::Int, m::Int, edges::Vector{Tuple{Int,Int}}, parts::Vector{Int64})::Matrix{Float64}
    unique_parts = sort(unique(parts))
    k = length(unique_parts)
    pos = zeros(n + m, 2)
    # Per-partition sub-circle radius: area ∝ n_vertices
    base_r = 0.3
    radii = Dict(p => base_r * sqrt(count(==(p), parts)) for p in unique_parts)
    r_max = maximum(values(radii))
    # Outer radius: adjacent sub-circles stay separated by at least 0.3 units
    R = k == 1 ? 0.0 : (r_max + 0.15) / sin(π / k) * 1.2
    for (ci, p) in enumerate(unique_parts)
        θ_c = 2π * (ci - 1) / k
        cx, cy = R * cos(θ_c), R * sin(θ_c)
        r = radii[p]
        verts = findall(==(p), parts)
        nv = length(verts)
        for (j, v) in enumerate(verts)
            θ = nv == 1 ? 0.0 : 2π * (j - 1) / nv
            pos[v, 1] = cx + r * cos(θ)
            pos[v, 2] = cy + r * sin(θ)
        end
    end
    # fake vertices: centroid of their connected real vertices
    for i in 1:m
        fv = n + i
        connected = [v for (u, v) in edges if u == fv]
        if isempty(connected)
            pos[fv, :] .= 0.0
        else
            pos[fv, 1] = sum(pos[v, 1] for v in connected) / length(connected)
            pos[fv, 2] = sum(pos[v, 2] for v in connected) / length(connected)
        end
    end
    _repel_fake!(pos, n)
    return pos
end

# Iteratively push all nodes apart until no pair is closer than min_dist.
function _repel_all!(pos::Matrix{Float64}; min_dist::Float64 = 0.2, iterations::Int = 80)
    n = size(pos, 1)
    for _ in 1:iterations
        moved = false
        for i in 1:n, j in (i+1):n
            dx, dy = pos[i,1]-pos[j,1], pos[i,2]-pos[j,2]
            d = sqrt(dx^2 + dy^2)
            if d < min_dist && d > 1e-9
                f = (min_dist - d) * 0.5 / d
                pos[i,1] += dx*f;  pos[i,2] += dy*f
                pos[j,1] -= dx*f;  pos[j,2] -= dy*f
                moved = true
            end
        end
        moved || break
    end
end

# Like _repel_all! but only fake vertices (indices n_real+1:end) are allowed to move.
function _repel_fake!(pos::Matrix{Float64}, n_real::Int; min_dist::Float64 = 0.2, iterations::Int = 100)
    n = size(pos, 1)
    for _ in 1:iterations
        moved = false
        for i in (n_real+1):n
            for j in 1:n
                i == j && continue
                dx, dy = pos[i,1]-pos[j,1], pos[i,2]-pos[j,2]
                d = sqrt(dx^2 + dy^2)
                if d < min_dist && d > 1e-9
                    f = (min_dist - d) / d
                    pos[i,1] += dx*f;  pos[i,2] += dy*f
                    moved = true
                end
            end
        end
        moved || break
    end
end

# Unique curvature magnitudes per hyperedge (cycled); signs are optimised at draw time.
const _HE_MAGS = [0.20, 0.32, 0.12, 0.28, 0.16, 0.36, 0.24, 0.08, 0.40, 0.18]

# Standard open-segment intersection test.
function _seg_cross(ax, ay, bx, by, cx, cy, dx, dy)
    d1x, d1y = bx-ax, by-ay
    d2x, d2y = dx-cx, dy-cy
    den = d1x*d2y - d1y*d2x
    abs(den) < 1e-12 && return false
    t = ((cx-ax)*d2y - (cy-ay)*d2x) / den
    u = ((cx-ax)*d1y - (cy-ay)*d1x) / den
    return 0.0 < t < 1.0 && 0.0 < u < 1.0
end

# Count crossings between all inter-hyperedge Bézier pairs, each approximated as
# two line segments (start→ctrl, ctrl→end).
function _count_crossings(pos, he_edges, curves)
    n   = 0
    ids = sort(collect(keys(he_edges)))
    for a in eachindex(ids), b in (a+1):lastindex(ids)
        hi, hj = ids[a], ids[b]
        ki, kj = curves[hi], curves[hj]
        for (fvi, vi) in he_edges[hi], (fvj, vj) in he_edges[hj]
            p0x,p0y = pos[vi, 1],pos[vi, 2];   p2x,p2y = pos[fvi,1],pos[fvi,2]
            dx, dy  = p2x-p0x, p2y-p0y;        L  = sqrt(dx^2+dy^2);  L  < 1e-9 && continue
            m1x = (p0x+p2x)/2 + ki*(-dy/L);    m1y = (p0y+p2y)/2 + ki*(dx/L)
            q0x,q0y = pos[vj, 1],pos[vj, 2];   q2x,q2y = pos[fvj,1],pos[fvj,2]
            ex, ey  = q2x-q0x, q2y-q0y;        Lj = sqrt(ex^2+ey^2); Lj < 1e-9 && continue
            m2x = (q0x+q2x)/2 + kj*(-ey/Lj);  m2y = (q0y+q2y)/2 + kj*(ex/Lj)
            n += (_seg_cross(p0x,p0y,m1x,m1y, q0x,q0y,m2x,m2y) ? 1 : 0) +
                 (_seg_cross(p0x,p0y,m1x,m1y, m2x,m2y,q2x,q2y) ? 1 : 0) +
                 (_seg_cross(m1x,m1y,p2x,p2y, q0x,q0y,m2x,m2y) ? 1 : 0) +
                 (_seg_cross(m1x,m1y,p2x,p2y, m2x,m2y,q2x,q2y) ? 1 : 0)
        end
    end
    return n
end

# Greedy sign-flip optimiser: for each hyperedge flip its curvature sign if doing so
# reduces the total inter-hyperedge Bézier crossing count.
function _optimize_curvature_signs(pos::Matrix{Float64},
                                   edges::Vector{Tuple{Int,Int}}, n_real::Int)
    he_edges = Dict{Int, Vector{Tuple{Int,Int}}}()
    for (fv, v) in edges
        push!(get!(he_edges, fv - n_real, Tuple{Int,Int}[]), (fv, v))
    end
    isempty(he_edges) && return he_edges  # empty → return empty Dict
    ids    = sort(collect(keys(he_edges)))
    curves = Dict(he => _HE_MAGS[mod1(he, length(_HE_MAGS))] for he in ids)
    for _ in 1:20
        improved = false
        for he in ids
            curr = _count_crossings(pos, he_edges, curves)
            curves[he] *= -1
            _count_crossings(pos, he_edges, curves) >= curr ? (curves[he] *= -1) : (improved = true)
        end
        improved || break
    end
    return curves
end

# Draw each edge as a quadratic Bézier with a curvature whose sign is chosen to
# minimise crossings between edges from different hyperedges.
function _draw_edges!(ax::Axis, pos::Matrix{Float64},
                      edges::Vector{Tuple{Int,Int}}, n_real::Int;
                      n_pts::Int = 40, he_colors::Dict{Int,Any} = Dict{Int,Any}())
    curves = _optimize_curvature_signs(pos, edges, n_real)
    ts = range(0.0, 1.0, n_pts)
    for (fv, v) in edges
        he  = fv - n_real
        k   = get(curves, he, 0.20)
        col = get(he_colors, he, :gray70)
        p0x, p0y = pos[v,  1], pos[v,  2]
        p2x, p2y = pos[fv, 1], pos[fv, 2]
        dx, dy = p2x - p0x, p2y - p0y
        L = sqrt(dx^2 + dy^2);  L < 1e-9 && continue
        cx = (p0x + p2x) / 2 + k * (-dy / L)
        cy = (p0y + p2y) / 2 + k * ( dx / L)
        xs = @. (1-ts)^2 * p0x + 2ts*(1-ts)*cx + ts^2*p2x
        ys = @. (1-ts)^2 * p0y + 2ts*(1-ts)*cy + ts^2*p2y
        lines!(ax, xs, ys; color = col, linewidth = 1)
    end
end

# Periodic Catmull-Rom spline through n control points (n×2 matrix).
# Returns a closed Vector{Point2f}.
function _catmull_rom(P::Matrix{Float64}; steps::Int = 30)::Vector{Point2f}
    n = size(P, 1)
    pts = Point2f[]
    for i in 1:n
        p0 = P[mod1(i-1, n), :];  p1 = P[mod1(i,   n), :]
        p2 = P[mod1(i+1, n), :];  p3 = P[mod1(i+2, n), :]
        for s in 0:(steps-1)
            t = s / steps;  t2 = t^2;  t3 = t^3
            push!(pts, Point2f(
                0.5*((2p1[1]) + (-p0[1]+p2[1])*t + (2p0[1]-5p1[1]+4p2[1]-p3[1])*t2 + (-p0[1]+3p1[1]-3p2[1]+p3[1])*t3),
                0.5*((2p1[2]) + (-p0[2]+p2[2])*t + (2p0[2]-5p1[2]+4p2[2]-p3[2])*t2 + (-p0[2]+3p1[2]-3p2[2]+p3[2])*t3)
            ))
        end
    end
    push!(pts, pts[1])
    return pts
end

# Smooth closed blob enclosing vertex positions vpos, offset outward by pad.
# Falls back to a circle for fewer than 3 vertices.
function _blob_curve(vpos::Matrix{Float64}; pad::Float64 = 0.12)::Vector{Point2f}
    nv = size(vpos, 1)
    cx = sum(vpos[:, 1]) / nv
    cy = sum(vpos[:, 2]) / nv
    ctrl = Matrix{Float64}(undef, nv, 2)
    for i in 1:nv
        dx, dy = vpos[i,1] - cx, vpos[i,2] - cy
        d = sqrt(dx^2 + dy^2)
        scale = d < 1e-9 ? 1.0 : (d + pad) / d
        ctrl[i, 1] = cx + dx * scale
        ctrl[i, 2] = cy + dy * scale
    end
    perm = sortperm([atan(ctrl[i,2] - cy, ctrl[i,1] - cx) for i in 1:nv])
    ctrl = ctrl[perm, :]
    if nv < 3
        rad = maximum(sqrt((vpos[i,1]-cx)^2 + (vpos[i,2]-cy)^2) for i in 1:nv) + pad
        return [Point2f(cx + rad*cos(θ), cy + rad*sin(θ)) for θ in range(0, 2π, 121)]
    end
    return _catmull_rom(ctrl)
end

"""
    plot_hypergraph(h::KaHyPar.HyperGraph) -> Figure

Draw a hypergraph using the GraphBased representation.

# Arguments
- `h`: hypergraph to visualize

# Returns
A `CairoMakie.Figure` with real vertices drawn as labelled circles (steel blue,
label "qi" in white) and hyperedge fake vertices as diamonds (orange) with their
weight displayed beside them, connected by gray edges. A legend in the top-right
corner distinguishes vertex and hyperedge node markers.
"""
function plot_hypergraph(h::KaHyPar.HyperGraph)::Figure
    n = Int(h.n_vertices)
    m = length(h.edge_indices) - 1
    edges = Tuple{Int,Int}[]
    for i in 1:m
        fv = n + i
        start_idx = Int(h.edge_indices[i]) + 1
        end_idx   = Int(h.edge_indices[i + 1])
        for vid_idx in start_idx:end_idx
            v = Int(h.hyperedges[vid_idx]) + 1
            push!(edges, (fv, v))
        end
    end
    pos = _spring_layout(n + m, edges)
    _repel_all!(pos)

    universal_color = :darkorange
    variant_color   = :crimson
    univ_edge_col   = :gray60
    var_edge_col    = :salmon

    node_colors = [Int(h.e_weights[i]) >= 0 ? universal_color : variant_color for i in 1:m]
    he_colors   = Dict{Int,Any}(i => (Int(h.e_weights[i]) >= 0 ? univ_edge_col : var_edge_col) for i in 1:m)

    fig = Figure(size = (800, 600))
    ax  = Axis(fig[1, 1], aspect = DataAspect(), title = "Hypergraph")
    hidedecorations!(ax)
    hidespines!(ax)
    _draw_edges!(ax, pos, edges, n; he_colors)
    scatter!(ax, pos[1:n, 1], pos[1:n, 2];
             marker = :circle, markersize = 30,
             color = :steelblue, strokecolor = :white, strokewidth = 1)
    n > 0 && text!(ax, [Point2f(pos[i, 1], pos[i, 2]) for i in 1:n];
                   text = ["q$i" for i in 1:n], fontsize = 10,
                   color = :white, align = (:center, :center))
    m > 0 && scatter!(ax, pos[n+1:n+m, 1], pos[n+1:n+m, 2];
             marker = :diamond, markersize = 14,
             color = node_colors, strokecolor = :white, strokewidth = 1)
    m > 0 && text!(ax, [Point2f(pos[n+i, 1], pos[n+i, 2]) for i in 1:m];
                   text = ["$(Int(h.e_weights[i]))" for i in 1:m], fontsize = 11,
                   color = node_colors, align = (:left, :bottom))
    Legend(fig[1, 1],
           [MarkerElement(marker = :circle,  color = :steelblue,      strokecolor = :white, strokewidth = 1, markersize = 18),
            MarkerElement(marker = :diamond, color = universal_color,  strokecolor = :white, strokewidth = 1, markersize = 14),
            MarkerElement(marker = :diamond, color = variant_color,    strokecolor = :white, strokewidth = 1, markersize = 14)],
           ["Vertex", "Hyperedge (universal)", "Hyperedge (variant)"];
           tellheight = false, tellwidth = false, halign = :right, valign = :top, margin = (10, 10, 10, 10))
    return fig
end

const _PART_PALETTE = [:steelblue, :firebrick, :seagreen, :darkorange,
                       :mediumpurple, :darkcyan, :deeppink, :olive]

"""
    plot_hypergraph_partition(h::KaHyPar.HyperGraph, parts::Vector{Int64}) -> Figure

Draw a partitioned hypergraph using the GraphBased representation.

# Arguments
- `h`: hypergraph to visualize
- `parts`: 0-indexed partition IDs, one per real vertex (as returned by `KaHyPar.partition`)

# Returns
A `CairoMakie.Figure` with a dotted smooth blob drawn around each partition's
real vertices (no fill). The blob shape is a periodic Catmull-Rom spline through
the vertex positions offset outward. Vertices are labelled "qi" in white inside
each circle; hyperedge weights are displayed beside each diamond. A legend in the
top-right corner distinguishes vertex and hyperedge node markers.
"""
function plot_hypergraph_partition(h::KaHyPar.HyperGraph, parts::Vector{Int64})::Figure
    n = Int(h.n_vertices)
    m = length(h.edge_indices) - 1
    edges = Tuple{Int,Int}[]
    for i in 1:m
        fv = n + i
        start_idx = Int(h.edge_indices[i]) + 1
        end_idx   = Int(h.edge_indices[i + 1])
        for vid_idx in start_idx:end_idx
            v = Int(h.hyperedges[vid_idx]) + 1
            push!(edges, (fv, v))
        end
    end
    pos = _partition_layout(n, m, edges, parts)
    fig = Figure(size = (800, 600))
    ax  = Axis(fig[1, 1], aspect = DataAspect(), title = "Hypergraph Partition (number of partitions=$(maximum(parts)+1))")
    hidedecorations!(ax)
    hidespines!(ax)
    for (ci, p) in enumerate(sort(unique(parts)))
        color = _PART_PALETTE[mod1(ci, length(_PART_PALETTE))]
        verts = findall(==(p), parts)
        blob  = _blob_curve(pos[verts, :]; pad = 0.1)
        lines!(ax, blob; color = color, linewidth = 2, linestyle = :dot)
    end
    _draw_edges!(ax, pos, edges, n)
    scatter!(ax, pos[1:n, 1], pos[1:n, 2];
             marker = :circle, markersize = 30,
             color = :steelblue, strokecolor = :white, strokewidth = 1)
    n > 0 && text!(ax, [Point2f(pos[i, 1], pos[i, 2]) for i in 1:n];
                   text = ["q$i" for i in 1:n], fontsize = 10,
                   color = :white, align = (:center, :center))
    m > 0 && scatter!(ax, pos[n+1:n+m, 1], pos[n+1:n+m, 2];
             marker = :diamond, markersize = 14,
             color = :darkorange, strokecolor = :white, strokewidth = 1)
    m > 0 && text!(ax, [Point2f(pos[n+i, 1], pos[n+i, 2]) for i in 1:m];
                   text = ["$(Int(h.e_weights[i]))" for i in 1:m], fontsize = 11,
                   color = :darkorange, align = (:left, :bottom))
    Legend(fig[1, 1],
           [MarkerElement(marker = :circle,  color = :steelblue,  strokecolor = :white, strokewidth = 1, markersize = 18),
            MarkerElement(marker = :diamond, color = :darkorange, strokecolor = :white, strokewidth = 1, markersize = 14)],
           ["Vertex", "Hyperedge"];
           tellheight = false, tellwidth = false, halign = :right, valign = :top, margin = (10, 10, 10, 10))
    return fig
end

"""
    plot_hyperedge_frequency(freq::Dict{Vector{Int}, Float64}) -> Figure

Visualize hyperedge frequency as a two-panel figure: a bar chart of frequencies on top
and an incidence heatmap (hyperedge × vertex) below, sorted by descending frequency.

# Arguments
- `freq`: mapping from sorted 1-indexed vertex vectors to their fractional frequency,
  as returned by normalizing `hyperedge_frequency`

# Returns
A `CairoMakie.Figure` with the two linked panels.
"""
function plot_hyperedge_frequency(freq::Dict{Vector{Int}, Float64})::Figure
    sorted = sort(collect(freq), by = x -> x[2], rev = true)
    edges  = first.(sorted)
    freqs  = last.(sorted)
    n_e    = length(edges)
    n_v    = maximum(maximum.(edges))

    mat = zeros(n_e, n_v)
    for (j, edge) in enumerate(edges)
        for v in edge
            mat[j, v] = 1.0
        end
    end

    fig    = Figure(size = (900, 400))
    ax_bar = Axis(fig[1, 1], ylabel = "Frequency", xticklabelsvisible = false)
    ax_hm  = Axis(fig[2, 1], xlabel = "Hyperedge rank", ylabel = "Vertex",
                  yticks = (1:n_v, string.(1:n_v)))

    barplot!(ax_bar, 1:n_e, freqs)
    heatmap!(ax_hm, mat, colormap = :Blues)

    linkxaxes!(ax_bar, ax_hm)
    rowsize!(fig.layout, 1, Relative(0.25))
    return fig
end


end # module
