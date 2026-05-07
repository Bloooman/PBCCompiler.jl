"""
Plotting scaffolding for circuit visualization.

The actual implementation is provided by the PBCCompilerMakieExt extension
when Makie is loaded.
"""

"""
    circuitplot(circuit::Circuit)
    circuitplot!(ax, circuit::Circuit)

Create a visual representation of a quantum circuit.

Requires the Makie package to be loaded. Each operation is drawn as a colored
rectangle spanning the qubits it affects. Measurements show the classical bit
index, and conditional operations show their dependency.

# Plot Attributes
- `gatewidth`: Width of gate rectangles (default: 0.8)
- `qubitspacing`: Vertical spacing between qubit lines (default: 1.0)
- `wirecolor`: Color of qubit wire lines (default: :black)
- `wirelinewidth`: Width of qubit wire lines (default: 1.0)

# Gate Colors (by variant)
- `paulicolor`: Color for Pauli gates
- `measurementcolor`: Color for measurements
- `exphalfpicolor`: Color for exp(i*pi/2*P) gates
- `expquatpicolor`: Color for exp(i*pi/4*P) gates
- `expeighpicolor`: Color for exp(i*pi/8*P) gates
- `conditionalcolor`: Color for PauliConditional gates
- `bitconditionalcolor`: Color for BitConditional gates
"""
function circuitplot end

function circuitplot! end

"""
    circuitplot_axis(subfig, circuit; kwargs...)

Create a complete Makie figure panel with a circuit plot and appropriate axis settings.

Returns a tuple of (subfig, axis, plot).
"""
function circuitplot_axis end

##
"""
    plot_histogram(data)

Plot a histogram of integer-valued `data` using CairoMakie.

- Each bar is centered over its integer value, with the x-axis label appearing
  directly below the bar (not at bar edges).
- X-axis tick labels are shown as bitstrings.
- The frequency count is labeled on top of each bar.
"""
function plot_histogram end

"""
    plot_interaction(weights::AbstractMatrix{<:Real}) -> Figure

Plot the interaction graph G(V,E) from an adjacency matrix.

Vertices V correspond to indices 1..n (row/column indices of `weights`).
Element `weights[i,j]` is the edge weight of the edge between vertex i and vertex j.
Only edges with nonzero weight are drawn; edge thickness and opacity scale with weight.
Edge weights are labelled at edge midpoints.

Returns a CairoMakie `Figure`.
"""
function plot_interaction end

"""
    plot_weight_histogram(g::SimpleWeightedGraph; bins=nothing)

Plot a histogram of edge weights from `g`, assuming integer edge weights.

When `bins` is `nothing` (default), each distinct integer weight gets its own bar.
When `bins` is an integer, weights are grouped into approximately that many bins
with integer-aligned boundaries. Returns a `Makie.Figure`.
"""
function plot_weight_histogram end

"""
    plot_std_graph(g::SimpleWeightedGraph; colormap=Reverse(:RdBu), node_size=30)

Plot the output of `weight_std_graph` as a network diagram using a circular layout.

Nodes are labeled 1–n. Edge color encodes the edge weight (std dev) on a
blue-to-red scale, with a colorbar on the right. Returns a `Makie.Figure`.
"""
function plot_std_graph end
##
"""
    circuitplot(circuit::Circuit) -> PyObject

Render a Circuit as a Qiskit circuit diagram.

# Arguments
- `circuit`: sequence of `CircuitOp` operations to visualize

# Returns
A matplotlib `Figure` with the circuit diagram. Qubit count and
classical bit count are inferred from the maximum indices in the circuit.
"""
function circuitplot_qiskit end
##
"""
    plot_partition(g::SimpleWeightedGraph, part::Vector{Int32}) -> Figure

Display a circular-layout plot of `g` with nodes colored by partition assignment.

# Arguments
- `g`: the graph that was partitioned
- `part`: partition assignment vector as returned by `METIS_partition`

# Returns
A `CairoMakie.Figure` with the partition plot.
"""
function plot_partition end
##
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
function plot_hypergraph end

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
function plot_hypergraph_partition end
