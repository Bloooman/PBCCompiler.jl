"""
Helper functions to Perform Joint measurement on Pauli Product Measurements
"""

function validate_input(circuit::Circuit, input::Stabilizer)
    if get_circuit_width(circuit)<length(input[1])
        throw(ArgumentError("Input state has more qubits than circuit input"))
    else
        nothing
    end
end

"""Perform Commutativity/Dependencies Check on Pauli Product Measurement with Stabilizer list"""
function check_PPM(s::Stabilizer,op::CircuitOp.Type, num_qubits::Int)
    if !isa_variant(op,CircuitOp.Measurement)
        return nothing
    else
        paulilen = num_qubits
        pauli=embed(paulilen, op.qubits, op.pauli)
        return project!(copy(s),pauli)
    end
end

"""
Perform Joint Measurement on CircuitOp if it's a CircuitOp.Measurement
Store results as corresponding measurement type: classical_random_result, classical_deterministic_result, quantum_result
"""
function get_measurement_result(state::S, op::CircuitOp.Type) where S <: AbstractSimState
    @debug "Measuring" op _group=:api
    ms=state.memory_state
    s=ms.stabilizer_group
    num_qubits = get_circuit_width(state.circuit)
    len=length(s)
    projection = check_PPM(s, op, num_qubits)
    if projection === nothing
        return nothing
    else
        if projection[3] === nothing
            if projection[2]<=len
                result = rand(Bool[0,1])
                return (classical_random_result(op.pauli, result),projection[2])
            else
                (quantum_state, result) = quantum_measurement(state, op, num_qubits)
                return (quantum_result(op.pauli, result),projection[2], quantum_state)
            end
        else
            result = Bool(projection[3]>>1)
            return (classical_deterministic_result(op.pauli, result),projection[2])
        end
    end
end

"""Perform quantum measurement simulation on given state using QuantumClifford.jl backend"""
function quantum_measurement(state::ComputerState, op::CircuitOp.Type, num_qubits::Int)
    magicqubits = collect(num_qubits-state.num_gadgets+1:num_qubits)
    quantum_state = state.memory_state.quantum_memory
    if quantum_state === nothing
        throw(ArgumentError("Magic State not initiated"))
    else
        real_p=op.pauli[magicqubits]
        bit_result = projectrand!(quantum_state, real_p)[2]
        result=Bool(bit_result>>1)
        return (quantum_state, result)
    end
end

"""Perform quantum measurement simulation on given state using classical sampling according to weight determined by user named outcome_probs"""
function quantum_measurement(state::DummyState, op::CircuitOp.Type, num_qubits::Int)
    quantum_state = state.memory_state.quantum_memory
    result = wsample([false,true],state.outcome_probs)
    return (quantum_state,result)
end

"""Resolve conditional circuit operations defined by CircuitOp.BitConditional"""
function resolve_conditionals(state::S) where S <: AbstractSimState
    circuit=state.circuit
    ms=state.memory_state
    creg=ms.classical_register
    index=find_BitConditional_indices(circuit)
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
