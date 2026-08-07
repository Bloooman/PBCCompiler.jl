
"""
Functions for querying which qubits are affected by circuit operations.
"""
##
using Moshi.Match: @match
"""
    affectedqubits(op::CircuitOp.Type) -> Vector{Int}

Return the sorted list of qubit indices affected by a circuit operation.

# Examples
```julia
op = Pauli(P"XY", [1, 2])
affectedqubits(op)  # returns [1, 2]

op = PauliConditional(P"X", [1], P"Z", [3])
affectedqubits(op)  # returns [1, 3]
```
"""
function affectedqubits(op::CircuitOp.Type)
    qubits = @match op begin
        CircuitOp.Measurement(pauli, bit, qubits) => qubits
        CircuitOp.Pauli(pauli, qubits) => qubits
        CircuitOp.ExpHalfPiPauli(pauli, qubits) => qubits
        CircuitOp.ExpQuatPiPauli(pauli, qubits) => qubits
        CircuitOp.ExpEighPiPauli(pauli, qubits) => qubits
        CircuitOp.PrepMagic(qubit, qubits) => vcat([qubit], qubits)
        CircuitOp.PauliConditional(cp, cq, tp, tq) => vcat(cq, tq)
        CircuitOp.BitConditional(inner_op, bit) => affectedqubits(inner_op)
    end
    return sort(unique(qubits))
end

"""
    max_affected_qubit(op::CircuitOp.Type) -> Int

Return the largest qubit index affected by a circuit operation, without
allocating the sorted index list that [`affectedqubits`](@ref) builds.
"""
function max_affected_qubit(op::CircuitOp.Type)
    @match op begin
        CircuitOp.PrepMagic(qubit, qubits) => max(qubit, maximum(qubits))
        CircuitOp.PauliConditional(cp, cq, tp, tq) => max(maximum(cq), maximum(tq))
        CircuitOp.BitConditional(inner_op, bit) => max_affected_qubit(inner_op)
        _ => maximum(op.qubits)
    end
end

"""
    affectedqubits(circuit::Circuit) -> Vector{Int}

Return the sorted list of all qubit indices affected by any operation in the circuit.
"""
function affectedqubits(circuit::Circuit)
    isempty(circuit) && return Int[]
    # Accumulate into one vector rather than `reduce(vcat, ...)`, which
    # reallocates and recopies the whole result once per operation
    all_qubits = Int[]
    for op in circuit
        append!(all_qubits, affectedqubits(op))
    end
    return sort!(unique!(all_qubits))
end
