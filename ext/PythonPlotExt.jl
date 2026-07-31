# ext/PythonPlotExt.jl
module PythonPlotExt

using PBCCompiler
using PythonCall
using PBCCompiler: Circuit, CircuitOp
using QuantumClifford: PauliOperator, @P_str
using Moshi.Data: isa_variant
using Moshi.Match: @match

import PBCCompiler: circuitplot_qiskit

##
const P = typeof(P"XYZ")

"""
    circuitplot(circuit::Circuit) -> PyObject

Render a Circuit as a Qiskit circuit diagram.

# Arguments
- `circuit`: sequence of `CircuitOp` operations to visualize

# Returns
A matplotlib `Figure` with the circuit diagram. Qubit count and
classical bit count are inferred from the maximum indices in the circuit.
"""
function circuitplot_qiskit(circuit::Circuit)
    n_q = Ref(0)
    n_c = Ref(0)

    function scan(op)
        @match op begin
            CircuitOp.Measurement(; bit, qubits) => begin
                n_q[] = max(n_q[], maximum(qubits))
                n_c[] = max(n_c[], bit)
            end
            CircuitOp.ExpHalfPiPauli(; qubits) => begin
                n_q[] = max(n_q[], maximum(qubits))
            end
            CircuitOp.ExpQuatPiPauli(; qubits) => begin
                n_q[] = max(n_q[], maximum(qubits))
            end
            CircuitOp.ExpEighPiPauli(; qubits) => begin
                n_q[] = max(n_q[], maximum(qubits))
            end
            CircuitOp.PauliConditional(; control_qubits, target_qubits) => begin
                n_q[] = max(n_q[], maximum(control_qubits), maximum(target_qubits))
            end
            CircuitOp.BitConditional(; op=inner, bit) => begin
                n_c[] = max(n_c[], bit)
                scan(inner)
            end
        end
    end
    for op in circuit
        scan(op)
    end

    qc_mod = pyimport("qiskit.circuit")
    Gate = qc_mod.Gate
    IfElseOp = qc_mod.IfElseOp
    qc = qc_mod.QuantumCircuit(n_q[], n_c[])

    # QuantumClifford uses '_' for identity in string(); normalize to 'I'
    pauli_chars(p::P) = [c == '_' ? 'I' : c for c in string(p) if c in "XYZ_"]
    py0(q::Vector{Int}) = pylist(q .- 1)  # Julia 1-based → Python 0-based

    # Return only the non-identity Pauli chars and their corresponding qubits
    function active_pauli_qubits(p::P, qubits::Vector{Int})
        chars = pauli_chars(p)
        mask = findall(!=('I'), chars)
        isempty(mask) && return (String(chars), qubits)
        return (String(chars[mask]), qubits[mask])
    end

    # Return the set of active (non-identity) qubit indices for an op
    function op_qubits(op)
        @match op begin
            CircuitOp.Measurement(; pauli, qubits) =>
                Set((active_pauli_qubits(pauli, qubits))[2])
            CircuitOp.ExpHalfPiPauli(; pauli, qubits) =>
                Set((active_pauli_qubits(pauli, qubits))[2])
            CircuitOp.ExpQuatPiPauli(; pauli, qubits) =>
                Set((active_pauli_qubits(pauli, qubits))[2])
            CircuitOp.ExpEighPiPauli(; pauli, qubits) =>
                Set((active_pauli_qubits(pauli, qubits))[2])
            CircuitOp.PauliConditional(; control_pauli, control_qubits, target_pauli, target_qubits) => begin
                (_, aq_ctrl) = active_pauli_qubits(control_pauli, control_qubits)
                (_, aq_tgt)  = active_pauli_qubits(target_pauli, target_qubits)
                Set(vcat(aq_ctrl, aq_tgt))
            end
            CircuitOp.BitConditional(; op=inner) => op_qubits(inner)
        end
    end

    half_pi_names    = Set{String}()
    quat_pi_names    = Set{String}()
    eigh_pi_names    = Set{String}()
    pauli_cond_names = Set{String}()
    meas_names       = Set{String}()

    function make_gate(op)
        @match op begin
            CircuitOp.ExpHalfPiPauli(; pauli, qubits) => begin
                (label, active_q) = active_pauli_qubits(pauli, qubits)
                name = "\$R_{$(label)}\$\n\$\\pi/2\$"
                push!(half_pi_names, name)
                (Gate(name, length(active_q), pylist([])), active_q)
            end
            CircuitOp.ExpQuatPiPauli(; pauli, qubits) => begin
                (label, active_q) = active_pauli_qubits(pauli, qubits)
                name = "\$R_{$(label)}\$\n\$\\pi/4\$"
                push!(quat_pi_names, name)
                (Gate(name, length(active_q), pylist([])), active_q)
            end
            CircuitOp.ExpEighPiPauli(; pauli, qubits) => begin
                (label, active_q) = active_pauli_qubits(pauli, qubits)
                name = "\$R_{$(label)}\$\n\$\\pi/8\$"
                push!(eigh_pi_names, name)
                (Gate(name, length(active_q), pylist([])), active_q)
            end
            CircuitOp.Measurement(; pauli, bit, qubits) => begin
                # Label from the *filtered* Pauli so its letters line up 1:1 with
                # the gate's qubit arguments (and hence the wire index labels).
                (label, active_q) = active_pauli_qubits(pauli, qubits)
                name = "M[$(label)]→c$(bit)"
                push!(meas_names, name)
                (Gate(name, length(active_q), pylist([])), active_q)
            end
            CircuitOp.PauliConditional(; control_pauli, control_qubits, target_pauli, target_qubits) => begin
                (ctrl_label, active_ctrl) = active_pauli_qubits(control_pauli, control_qubits)
                (tgt_label, active_tgt)   = active_pauli_qubits(target_pauli, target_qubits)
                all_q = vcat(active_ctrl, active_tgt)
                name = "C($(ctrl_label))\n↓\nT($(tgt_label))"
                push!(pauli_cond_names, name)
                (Gate(name, length(all_q), pylist([])), all_q)
            end
            CircuitOp.BitConditional(; op=inner) => begin
                make_gate(inner)
            end
        end
    end

    # Greedy layer grouping: consecutive ops with no active-qubit overlap share a layer
    layers = Vector{Vector{eltype(circuit)}}()
    current_layer = eltype(circuit)[]
    current_qubits = Set{Int}()
    for op in circuit
        qs = op_qubits(op)
        if isempty(intersect(qs, current_qubits))
            push!(current_layer, op)
            union!(current_qubits, qs)
        else
            push!(layers, current_layer)
            current_layer = [op]
            current_qubits = qs
        end
    end
    isempty(current_layer) || push!(layers, current_layer)

    for (i, layer) in enumerate(layers)
        i > 1 && qc.barrier()
        for op in layer
            if isa_variant(op, CircuitOp.BitConditional)
                (gate, q) = make_gate(op.op)
                body = qc_mod.QuantumCircuit(length(q))
                body.append(gate, pylist(collect(0:length(q)-1)))
                if_op = IfElseOp((qc.clbits[op.bit - 1], pybool(true)), body, pybuiltins.None)
                qc.append(if_op, py0(q))
            else
                (gate, q) = make_gate(op)
                qc.append(gate, py0(q))
            end
        end
    end

    display_colors = pydict()
    for name in half_pi_names
        display_colors[pystr(name)] = pytuple([pystr("#808080"), pystr("#000000")])
    end
    for name in quat_pi_names
        display_colors[pystr(name)] = pytuple([pystr("#FF8C00"), pystr("#000000")])
    end
    for name in eigh_pi_names
        display_colors[pystr(name)] = pytuple([pystr("#FFD700"), pystr("#000000")])
    end
    for name in pauli_cond_names
        display_colors[pystr(name)] = pytuple([pystr("#FF8C00"), pystr("#000000")])
    end
    for name in meas_names
        display_colors[pystr(name)] = pytuple([pystr("#2196F3"), pystr("#000000")])
    end
    style = pydict(Dict("displaycolor" => display_colors))

    return qc.draw(output=pystr("mpl"), style=style, plot_barriers=pybool(false))
end

end
