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
    pauli_label(p::P) = String(pauli_chars(p))
    py0(q::Vector{Int}) = pylist(q .- 1)  # Julia 1-based → Python 0-based

    # Strip leading/trailing identity from both Pauli chars and qubit list
    function trim_pauli(p::P, qubits::Vector{Int})
        chars = pauli_chars(p)
        lo = findfirst(!=('I'), chars)
        hi = findlast(!=('I'), chars)
        isnothing(lo) && return (String(chars), qubits)
        return (String(chars[lo:hi]), qubits[lo:hi])
    end

    half_pi_names    = Set{String}()
    quat_pi_names    = Set{String}()
    eigh_pi_names    = Set{String}()
    pauli_cond_names = Set{String}()
    meas_names       = Set{String}()

    function make_gate(op)
        @match op begin
            CircuitOp.ExpHalfPiPauli(; pauli, qubits) => begin
                (label, active_q) = trim_pauli(pauli, qubits)
                name = "\$R_{$(label)}\$\n\$\\pi/2\$"
                push!(half_pi_names, name)
                (Gate(name, length(active_q), pylist([])), active_q)
            end
            CircuitOp.ExpQuatPiPauli(; pauli, qubits) => begin
                (label, active_q) = trim_pauli(pauli, qubits)
                name = "\$R_{$(label)}\$\n\$\\pi/4\$"
                push!(quat_pi_names, name)
                (Gate(name, length(active_q), pylist([])), active_q)
            end
            CircuitOp.ExpEighPiPauli(; pauli, qubits) => begin
                (label, active_q) = trim_pauli(pauli, qubits)
                name = "\$R_{$(label)}\$\n\$\\pi/8\$"
                push!(eigh_pi_names, name)
                (Gate(name, length(active_q), pylist([])), active_q)
            end
            CircuitOp.Measurement(; pauli, bit, qubits) => begin
                name = "M[$(pauli_label(pauli))]→c$(bit)"
                push!(meas_names, name)
                (Gate(name, length(qubits), pylist([])), qubits)
            end
            CircuitOp.PauliConditional(; control_pauli, control_qubits, target_pauli, target_qubits) => begin
                all_q = vcat(control_qubits, target_qubits)
                name = "C($(pauli_label(control_pauli)))\n↓\nT($(pauli_label(target_pauli)))"
                push!(pauli_cond_names, name)
                (Gate(name, length(all_q), pylist([])), all_q)
            end
            CircuitOp.BitConditional(; op=inner) => begin
                make_gate(inner)
            end
        end
    end

    for op in circuit
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

    return qc.draw(output=pystr("mpl"), style=style)
end

end
