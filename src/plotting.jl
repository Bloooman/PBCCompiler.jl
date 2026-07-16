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
    plot_histogram(data; nbits=nothing)

Plot a histogram of integer-valued `data` using CairoMakie.

- Each bar is centered over its integer value, with the x-axis label appearing
  directly below the bar (not at bar edges).
- X-axis tick labels are shown as bitstrings, zero-padded to `nbits`. If `nbits`
  is not given, it defaults to the number of bits needed to represent the
  largest observed value — note this under-pads whenever the true bit width
  isn't determinable from `data` alone (e.g. all outcomes happen to be 0),
  so pass `nbits` explicitly whenever the true bit width is known.
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
function plot_graph end

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
function plot_hyperedge_frequency end

##
"""
    plot_cooccurrence(h::KaHyPar.HyperGraph) -> Figure

Co-occurrence weighted-adjacency heatmap with hierarchical leaf ordering.

# Arguments
- `h`: KaHyPar hypergraph (`n_vertices`, `edge_indices`, `hyperedges`, `e_weights`)

# Returns
CairoMakie Figure. Rows/columns are nodes reordered by average-linkage clustering
on dissimilarity 1 - W/max(W). Tick labels show original 0-based node ids.
Diagonal is zeroed so the colormap reflects only inter-node co-occurrence.

W_ij = sum of `e_weights[e]` over all hyperedges `e` containing both node i and node j
(raw sum, not normalized by hyperedge size or count).

# Reading the plot
- **Bright cell at (i, j)**: nodes i and j co-appear in one or more heavy hyperedges.
- **Bright diagonal block**: a group of nodes that all co-appear frequently with each other —
  this is a community. The reordering groups such nodes together so communities emerge
  as contiguous blocks along the diagonal.
- **Off-diagonal brightness**: nodes i and j share hyperedges despite being in different
  communities (e.g. bridge nodes or cut hyperedges).
- Tick labels show the original 0-based node id at each reordered position.
"""
function plot_cooccurrence end

##
"""
    plot_incidence(h::KaHyPar.HyperGraph, parts::Vector{Int}) -> Figure

Incidence matrix heatmap with nodes and hyperedges reordered to expose partition structure.

# Arguments
- `h`: KaHyPar hypergraph
- `parts`: partition assignment vector of length `h.n_vertices` (0-based block ids, as returned by KaHyPar)

# Ordering
- **Nodes (x-axis)**: sorted by `parts` block id — all nodes in block 0 first, then
  block 1, etc. Red vertical lines mark the block boundaries.
- **Hyperedges (y-axis)**: non-cut hyperedges first, grouped by their partition block
  (well-defined since all their nodes are in the same block) and sorted by `e_weights`
  descending within each block. Cut hyperedges are placed at the top, also sorted by
  `e_weights` descending.

# Returns
CairoMakie Figure. Cell color encodes `e_weights[e]` for member nodes, white for
non-members. Cut hyperedges (members spanning >1 partition block) are highlighted
with a semi-transparent red overlay.

# Reading the plot
- **Colored cell at column i, row e**: node i is a member of hyperedge e; darkness
  reflects the hyperedge weight `e_weights[e]`.
- **White cell**: node i is not a member of hyperedge e.
- **Colored band confined within one vertical block**: a non-cut hyperedge — all its
  nodes are in the same partition.
- **Red-tinted row**: a cut hyperedge whose members span more than one partition block.
  Its colored cells will cross at least one red vertical line. These are the edges
  penalized by the min-cut objective.
"""
function plot_incidence end
