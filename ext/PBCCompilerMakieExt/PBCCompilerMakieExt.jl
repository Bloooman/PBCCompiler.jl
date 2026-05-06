module PBCCompilerMakieExt

using Makie
using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, _affectedqubits
using Graphs, SimpleWeightedGraphs
using Moshi.Match: @match

import PBCCompiler: circuitplot, circuitplot!, circuitplot_axis, plot_histogram, plot_interaction, plot_weight_histogram, plot_std_graph, plot_partition

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
    all_qubits = _affectedqubits(circuit)
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

        qubits = _affectedqubits(op)
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
        all_qubits = _affectedqubits(circuit)
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
    plot_interaction(weights::AbstractMatrix{<:Real}) -> Figure

Plot the interaction graph G(V,E) from an adjacency matrix.

Vertices V correspond to indices 1..n (row/column indices of `weights`).
Element `weights[i,j]` is the edge weight of the edge between vertex i and vertex j.
Only edges with nonzero weight are drawn; edge thickness and opacity scale with weight.
Edge weights are labelled at edge midpoints.

Returns a CairoMakie `Figure`.
"""
function plot_interaction(weights::AbstractMatrix{<:Real})
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

    max_w = maximum(weights)

    for i in 1:n, j in i+1:n
        w = weights[i, j]
        iszero(w) && continue
        frac = w / max_w
        lines!(ax, [xs[i], xs[j]], [ys[i], ys[j]];
               color=(:gray20, 0.3 + 0.7 * frac), linewidth=1 + 3 * frac)
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


end # module
