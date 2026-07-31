# ext/PBCCompilerQuantikzExt.jl
module PBCCompilerQuantikzExt

using PBCCompiler
using Quantikz
using PBCCompiler: Circuit, CircuitOp
using QuantumClifford: PauliOperator, @P_str
using Moshi.Match: @match

import PBCCompiler: circuitplot_quantikz, circuitstring_quantikz

const P = typeof(P"XYZ")

# Default `fill=` colors, one per operation kind. Chosen to match the palette the
# Qiskit renderer uses so the two backends stay visually comparable.
const DEFAULT_COLORS = Dict(
    :measurement => "cyan!25",
    :halfpi      => "black!15",
    :quatpi      => "orange!45",
    :eighpi      => "yellow!45",
    :control     => "red!30",
    :target      => "orange!45",
)

# The quantikz environment options. Deliberately *not* Quantikz.jl's default,
# whose `transparent` key suppresses every `fill=` we emit.
const DEFAULT_OPTIONS = "row sep={0.8cm,between origins}"

##
# One drawn cell: a stack of (Pauli letter, subscript) parts plus a fill color.
# The stack is a list because a single qubit can carry more than one role, e.g.
# both the control and the target leg of a `PauliConditional`.
struct PauliBox
    parts::Vector{Tuple{Char,String}}
    fill::String
end

function body(b::PauliBox)
    join((isempty(sub) ? string(l) : "$(l)_{$(sub)}" for (l, sub) in b.parts), "\\,")
end

cell(b::PauliBox) = isempty(b.fill) ? "\\gate{$(body(b))}" :
                    "\\gate[style={fill=$(b.fill)}]{$(body(b))}"

"""
A single PBC operation drawn as one labelled box per *supported* qubit, joined by
a vertical rail. Qubits the operation does not act on stay clear, unlike the
single spanning rectangle that Quantikz's own multi-qubit ops (and Qiskit) draw.
"""
struct PauliBoxOp <: Quantikz.QuantikzOp
    boxes::Vector{PauliBox}
    targets::Vector{Int}
    """Classical bit to draw a `\\cwbend` to, or `nothing`."""
    bit::Union{Int,Nothing}
end

# Quantikz reserves the min:max span of these when packing columns, which is what
# keeps another operation out of the wires this one's rail passes through.
Quantikz.affectedqubits(op::PauliBoxOp) = op.targets
Quantikz.affectedbits(op::PauliBoxOp) = isnothing(op.bit) ? Int[] : [op.bit]

function Quantikz.update_table!(qtable, step, op::PauliBoxOp)
    table = Quantikz.qubitsview(qtable)
    perm = sortperm(op.targets)
    targets = op.targets[perm]
    boxes = op.boxes[perm]

    table[targets[1], step] = cell(boxes[1])
    # `\vqw` offsets are relative to the previously linked cell, not to the top one.
    prev = targets[1]
    for (t, b) in zip(@view(targets[2:end]), @view(boxes[2:end]))
        table[t, step] = cell(b) * "\\vqw{$(prev - t)}"
        prev = t
    end

    if !isnothing(op.bit)
        layoutbit = if Quantikz.classicalbitslayout[] == :expanded
            op.bit
        elseif Quantikz.classicalbitslayout[] == :compressed
            1
        else
            Quantikz.classicalbitslayout_error()
        end
        anchor = targets[end]
        Quantikz.bitsview(qtable)[layoutbit, step] =
            "\\cwbend{$(-(qtable.qubits - anchor) - qtable.ancillaries - layoutbit)}"
    end
    qtable
end

##
# QuantumClifford writes identity as '_' in `string(p)`; normalize to 'I'.
_pauli_chars(p::P) = [c == '_' ? 'I' : c for c in string(p) if c in "XYZ_"]

"""
Drop identity factors, returning the surviving Pauli letters paired with the
qubits they act on. An all-identity Pauli has nothing to drop, so it is returned
whole rather than as an empty (undrawable) operation.
"""
function _active_pauli_qubits(p::P, qubits::Vector{Int})
    chars = _pauli_chars(p)
    mask = findall(!=('I'), chars)
    isempty(mask) && return (chars, qubits)
    return (chars[mask], qubits[mask])
end

"""Merge boxes that landed on the same qubit, concatenating their parts."""
function _merge_boxes(boxes::Vector{PauliBox}, targets::Vector{Int})
    allunique(targets) && return (boxes, targets)
    order = Int[]
    byqubit = Dict{Int,PauliBox}()
    for (b, t) in zip(boxes, targets)
        if haskey(byqubit, t)
            prev = byqubit[t]
            byqubit[t] = PauliBox(vcat(prev.parts, b.parts), prev.fill)
        else
            byqubit[t] = b
            push!(order, t)
        end
    end
    return ([byqubit[t] for t in order], order)
end

"""
Build the box stack for a Pauli acting on `qubits`.

`sub` is placed nowhere (`showangles=:none`), on the first box (`:first`), or on
every box (`:all`). The default is `:none` because the fill color already
encodes the rotation angle, and the label only widens the boxes.
"""
function _boxes(pauli::P, qubits::Vector{Int}, sub::String, fill::String, showangles::Symbol)
    (chars, active) = _active_pauli_qubits(pauli, qubits)
    boxes = map(enumerate(chars)) do (i, c)
        s = showangles === :all ? sub :
            showangles === :first && i == 1 ? sub : ""
        PauliBox([(c, s)], fill)
    end
    return (boxes, copy(active))
end

"""Append `sub` to the subscript of the box at position `i`."""
function _annotate(boxes::Vector{PauliBox}, i::Int, sub::String)
    b = boxes[i]
    (l, s) = b.parts[end]
    parts = copy(b.parts)
    parts[end] = (l, isempty(s) ? sub : "$(s),$(sub)")
    boxes[i] = PauliBox(parts, b.fill)
    return boxes
end

"""Translate one `CircuitOp` into a `PauliBoxOp`."""
function _to_op(op, colors, showangles)
    @match op begin
        CircuitOp.ExpHalfPiPauli(; pauli, qubits) => begin
            (boxes, targets) = _boxes(pauli, qubits, "\\pi/2", colors[:halfpi], showangles)
            PauliBoxOp(boxes, targets, nothing)
        end
        CircuitOp.ExpQuatPiPauli(; pauli, qubits) => begin
            (boxes, targets) = _boxes(pauli, qubits, "\\pi/4", colors[:quatpi], showangles)
            PauliBoxOp(boxes, targets, nothing)
        end
        CircuitOp.ExpEighPiPauli(; pauli, qubits) => begin
            (boxes, targets) = _boxes(pauli, qubits, "\\pi/8", colors[:eighpi], showangles)
            PauliBoxOp(boxes, targets, nothing)
        end
        CircuitOp.Measurement(; pauli, bit, qubits) => begin
            (boxes, targets) = _boxes(pauli, qubits, "", colors[:measurement], showangles)
            # The bend leaves from the bottom-most box; label the outcome bit there.
            _annotate(boxes, argmax(targets), "\\mathrm{c$(bit)}")
            PauliBoxOp(boxes, targets, bit)
        end
        CircuitOp.PauliConditional(; control_pauli, control_qubits, target_pauli, target_qubits) => begin
            (cb, cq) = _boxes(control_pauli, control_qubits, "\\mathrm{c}", colors[:control], :first)
            (tb, tq) = _boxes(target_pauli, target_qubits, "\\mathrm{t}", colors[:target], :first)
            (boxes, targets) = _merge_boxes(vcat(cb, tb), vcat(cq, tq))
            PauliBoxOp(boxes, targets, nothing)
        end
        CircuitOp.BitConditional(; op=inner, bit) => begin
            inner_op = _to_op(inner, colors, showangles)
            boxes = copy(inner_op.boxes)
            # "?c9" reads as "if c9", distinguishing a conditioning bit from a
            # measurement's output bit. Upright so it does not read as math.
            _annotate(boxes, argmin(inner_op.targets), "\\mathrm{?c$(bit)}")
            # An inner measurement keeps its outcome-bit label (added above) while
            # the bend is redirected to the conditioning bit.
            PauliBoxOp(boxes, inner_op.targets, bit)
        end
    end
end

_colors(overrides) = isempty(overrides) ? DEFAULT_COLORS : merge(DEFAULT_COLORS, overrides)

"""
Write `q1..qn` / `c1..cm` labels into the leading padding column of a table.

The classical section is a single wire under Quantikz's `:compressed` bit layout,
in which case it is labelled `c` rather than given a bit index.
"""
function _label_wires!(qtable)
    table = qtable.table
    for q in 1:qtable.qubits
        table[q, 1] = "\\lstick{\$q_{$(q)}\$}"
    end
    offset = qtable.qubits + qtable.ancillaries
    nbits = size(table, 1) - offset
    for b in 1:nbits
        table[offset + b, 1] = nbits == 1 ? "\\lstick{\$c\$}" : "\\lstick{\$c_{$(b)}\$}"
    end
    return qtable
end

##
"""
    circuitstring_quantikz(circuit::Circuit; kwargs...) -> String

Render a `Circuit` as a LaTeX `quantikz` environment.

Each operation is drawn as one labelled box per qubit it actually acts on, joined
by a vertical rail; qubits in between that the operation does not touch stay
clear. Box fill encodes the operation kind and the subscript the rotation angle,
classical bit, or conditioning bit.

# Keyword arguments
- `colors`: overrides for the per-kind fill colors, e.g. `Dict(:quatpi => "red!40")`.
  Keys are `:measurement`, `:halfpi`, `:quatpi`, `:eighpi`, `:control`, `:target`.
- `showangles`: `:none` (default) omits the rotation angle, since the box fill
  already encodes it; `:first` labels the topmost box of each operation, `:all`
  labels every box.
- `showlabels`: prefix each wire with its index, `q1..qn` and `c1..cm` (default
  `true`).
- `options`: the option list for the `quantikz` environment. Note that including
  `transparent` (Quantikz.jl's own default) suppresses all box fills.

# Returns
The LaTeX source as a `String`, ready to paste into a document.
"""
function circuitstring_quantikz(circuit::Circuit;
        colors::AbstractDict = Dict{Symbol,String}(),
        showangles::Symbol = :none,
        showlabels::Bool = true,
        options::AbstractString = DEFAULT_OPTIONS)
    showangles in (:first, :all, :none) ||
        throw(ArgumentError("showangles must be :first, :all or :none, got :$(showangles)"))
    # Quantikz cannot lay out a table with no rows; an empty circuit draws nothing.
    isempty(circuit) && return "\\begin{quantikz}[$(options)]\n\\end{quantikz}"
    c = _colors(colors)
    ops = [_to_op(op, c, showangles) for op in circuit]
    # Size the diagram from the *declared* qubits, so trailing wires that only
    # ever carry identity factors are still drawn.
    nqubits = maximum(PBCCompiler.max_affected_qubit, circuit)
    qtable = Quantikz.circuit2table(ops, nqubits)
    showlabels && _label_wires!(qtable)
    return Quantikz.table2string(qtable; quantikzoptions=String(options))
end

"""
    circuitplot_quantikz(circuit::Circuit; kwargs...)
    circuitplot_quantikz(circuit::Circuit, filename::AbstractString; kwargs...)

Render a `Circuit` as a quantikz circuit diagram.

Each operation is drawn as one labelled box per qubit it actually acts on, joined
by a vertical rail, so the support of an operation is readable off the diagram
directly. See [`circuitstring_quantikz`](@ref) for the keyword arguments
controlling colors and labels.

The one-argument form returns an image. The two-argument form writes to
`filename`, dispatching on its extension: `.tex` writes the LaTeX source, `.pdf`
the compiled document, anything else an image file.

# Keyword arguments
- `scale`: magnification of the rendered diagram (default `5`).
- plus every keyword of [`circuitstring_quantikz`](@ref).

# Examples
```julia
using Quantikz, PBCCompiler
circuitplot_quantikz(circuit)
circuitplot_quantikz(circuit, "circuit.pdf")
```
"""
function circuitplot_quantikz(circuit::Circuit; scale = 5, kwargs...)
    return Quantikz.string2image(circuitstring_quantikz(circuit; kwargs...); scale=scale)
end

function circuitplot_quantikz(circuit::Circuit, filename::AbstractString; scale = 5, kwargs...)
    str = circuitstring_quantikz(circuit; kwargs...)
    ext = lowercase(splitext(filename)[2])
    if ext == ".tex"
        write(filename, str)
    elseif ext == ".pdf"
        Quantikz.string2image(str; scale=scale, _workaround_savefile=filename)
    else
        Quantikz.save(filename, Quantikz.string2image(str; scale=scale))
    end
    return filename
end

end
