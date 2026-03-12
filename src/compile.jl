"""
Functions for compiling quantum circuits by moving measurement operations to the beginning of the circuit.
"""
struct PauliQubitMismatchError <: Exception
    msg::String
end

function validate_CircuitOp(op::CircuitOp.Type)
    p=affectedpaulis(op)
    q=affectedqubits(op)
    name=variant_name(op)
    if length(p) != length(q)
        throw(PauliQubitMismatchError("$name($p, $q): The length of the Pauli string is not the same as the number of affected qubits. Please check the input operation."))
    end
end

function validate_circuit(circuit::Circuit)
    for op in circuit
        validate_CircuitOp(op)
    end
end

function find_measurement_indices(circuit::Circuit)
    measurement_indices = []
    for (index, op) in enumerate(circuit)
        if isa_variant(op,CircuitOp.Measurement)
            push!(measurement_indices, index)
        end
    end
    return measurement_indices
end

function find_nonclifford_indices(circuit::Circuit)
    nonclifford_indices = []
    for (index, op) in enumerate(circuit)
        if isa_variant(op, CircuitOp.ExpEighPiPauli)
            push!(nonclifford_indices, index)
        end
    end
    return nonclifford_indices
end

function compilation(circuit::Circuit)
    validate_circuit(circuit)
    if find_nonclifford_indices(circuit) != []
        for index in find_nonclifford_indices(circuit)
            circuit=traversal(circuit, conjugate, :left, 1, index-1)
        end
    end
    for index in find_measurement_indices(circuit)
        circuit=traversal(circuit, conjugate, :left, 1, index-1)
    end
    return circuit
end
