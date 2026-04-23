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
