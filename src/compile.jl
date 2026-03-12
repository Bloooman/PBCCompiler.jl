"""
Functions for compiling quantum circuits by moving measurement operations to the beginning of the circuit.
"""

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
    for index in find_nonclifford_indices(circuit)
        circuit=traversal(circuit, conjugate, :left, 1, index-1)
    end
    for index in find_measurement_indices(circuit)
        circuit=traversal(circuit, conjugate, :left, 1, index-1)
    end
    return circuit
end
