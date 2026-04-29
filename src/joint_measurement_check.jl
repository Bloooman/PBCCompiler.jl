"""
Helper functions to Perform Joint measurement on Pauli Product Measurements
"""

function _validate_input(circuit::Circuit, input::Stabilizer)
    if _get_circuit_width(circuit)<length(input[1])
        throw(ArgumentError("Input state has more qubits than circuit input"))
    else
        nothing
    end
end

function _find_BitConditional_indices(circuit::Circuit)
    BitConditional_indices = []
    for (index, op) in enumerate(circuit)
        if isa_variant(op, CircuitOp.BitConditional)
            push!(BitConditional_indices, index)
        end
    end
    return BitConditional_indices
end

"""Perform Commutativity/Dependencies Check on Pauli Product Measurement with Stabilizer list"""
function _check_PPM(s::Stabilizer,op::CircuitOp.Type, num_qubits::Int)
    if !isa_variant(op,CircuitOp.Measurement)
        return nothing
    else
        Paulilen = num_qubits
        Pauli=embed(Paulilen, op.qubits, op.pauli)
        return project!(copy(s),Pauli)
    end
end

"""
0x00 denotes +1 eigenvalue, 0x02 denotes -1 eigenvalue
0 denotes +1 eigenvalue, 1 denotes -1 eigenvalue
false denotes +1 eigenvalue, true denotes -1 eigenvalue
"""

"""
Perform Joint Measurement on CircuitOp if it's a CircuitOp.Measurement
Store results as corresponding measurement type: classical_random_result, classical_deterministic_result, quantum_result
"""
function _get_measurement_result(s::Stabilizer, op::CircuitOp.Type, num_qubits::Int)
    len=length(s)
    projection = _check_PPM(s, op, num_qubits)
    if projection === nothing
        return nothing
    else
        if projection[3] === nothing
            if projection[2]<=len
                result = rand(Bool[0,1])
                return (classical_random_result(result),projection[2])
            else
                result = Bool(quantum_measurement(op)>>1)
                return (quantum_result(result),projection[2])
            end
        else
            result = Bool(projection[3]>>1)
            return (classical_deterministic_result(result),projection[2])
        end
    end
end

function quantum_measurement(op::CircuitOp.Type)
    print("Enter quantum measurement result: ")
    measurement_result = parse(Int,readline())
    if abs(measurement_result) != 1
        throw(ArgumentError("Measurement Result can only be +1 or -1!!!"))
    else
        if measurement_result == 1
            return 0x00
        elseif measurement_result == -1
            return 0x02
        else
            return nothing
        end
    end
end

"""Resolve conditional circuit operations defined by CircuitOp.BitConditional"""
function resolve_conditionals(compstate::ComputerState)
    CS=compstate
    circuit=CS.circuit
    MS=CS.memory_state
    creg=MS.classical_register
    index=_find_BitConditional_indices(circuit)
    for i in index
        @debug("Start resoving BitConditional at $i")
        operation=circuit[i]
        control_bit=creg[operation.bit]
        if control_bit !== nothing
            if control_bit
                @debug("$i has a controlled bit")
                splice!(circuit, i, [operation.op])
                @debug("$i Resolved")
                preprocess_circuit(circuit)
                break
            else
                deleteat!(circuit, i)
                @debug("No correction needed")
                preprocess_circuit(circuit)
                break
            end
        else
            @debug("Control Bit undetermined")
            nothing
        end
    end
end
