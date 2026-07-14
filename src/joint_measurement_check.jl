"""
Helper functions to check the first PPM in circuit, determine MeasurementResultType: ClassicalDetermRes, ClassicalRandomRes, QuantumRes
"""
##
using QuantumClifford: project!, Stabilizer, one, GeneralizedStabilizer, tensor_pow, apply!, pcT, projectrand!
using Moshi.Data: variant_name, isa_variant
using StatsBase: wsample
using Accessors: @reset
##

function validate_input(circuit::Circuit, input::Stabilizer)
    if get_circuit_width(circuit)<length(input[1])
        throw(ArgumentError("Input state has more qubits than circuit input"))
    else
        nothing
    end
end

function create_hadamard_basis_state(num_qubit::Int)
    n = num_qubit

    generators = one(Stabilizer, n; basis=:X)

    return Stabilizer(generators)
end

function create_magic_state(num_magic::Int)
    n=num_magic

    generators = GeneralizedStabilizer(create_hadamard_basis_state(n))

    T = tensor_pow(pcT,n)

    apply!(generators,T)

    return generators

end

"""
    check_PPM(s::Stabilizer,op::CircuitOp.Type, num_qubits::Int) -> Tuple{MixedStabilizer, Int64, Any}

Function that performs commutation and dependency checks on Pauli Product Measurement according to stabilizer list

It returns
- a stabilizer that might not be in canonical form
- the index of the row where the non-commuting operator was (that row is now equal to pauli; its phase is not updated and for a faithful measurement simulation it needs to be randomized by the user)
- and the result of the projection if there was no non-commuting operator (nothing otherwise)

In the case that PPM is commuting and independent, second field will return index = length(Stabilizer)+1, third field nothing
"""
function check_PPM(s::Stabilizer,op::CircuitOp.Type, num_qubits::Int)
    if !isa_variant(op,CircuitOp.Measurement)
        return nothing
    else
        Paulilen = num_qubits
        Pauli = embed(Paulilen, op.qubits, op.pauli)
        return project!(copy(s),Pauli)
    end
end

"""
    get_measurement_result(state::S, op::CircuitOp.Type) where S <: AbstractSimState -> Tuple{MeasurementResult, Int64, GeneralizedStabilizer}

Perform Joint Measurement on a CircuitOp.Measurement

It returns
- a MeasurementResult with corresponding MeasurementResultType (true denotes 0, false denotes 1)
- the index of the row where the non-commuting operator was (that row is now equal to pauli; its phase is not updated and for a faithful measurement simulation it needs to be randomized by the user)
- a GeneralizedStabilizer represents the quantum state after measurement
"""
function get_measurement_result(state::CompilerState, op::CircuitOp.Type)
    @debug "Measuring" op _group=:api
    rt = state.runtime
    check_list=state.stabilizer_group
    num_qubits = get_circuit_width(state.circuit)
    len=length(check_list)
    projection = check_PPM(check_list, op, num_qubits)
    if projection === nothing
        return nothing
    else
        if projection[3] === nothing
            if projection[2]<=len
                result = rand(Bool[0,1])
                @debug "This measurement outputs Classical Random Result" _group=:api
                q_1=[1:get_circuit_width(state.circuit);]
                Q_1=ExpQuatPiPauli(check_list[projection[2]],q_1)
                p_2=(-1)^result*op.pauli
                Q_2=ExpQuatPiPauli(p_2,op.qubits)
                pushfirst!(state.circuit,Q_1,Q_2,Q_1)
                preprocess_circuit(state.circuit)
                # Store the full-width Pauli so all result types share absolute qubit positions
                return (ClassicalRandomRes(embed(size(state.stabilizer_group)[2], op.qubits, op.pauli), result), state)
            else
                (rt, result) = quantum_measurement(rt, op, num_qubits)
                # The joint observable factors as (data part) ⊗ (magic part).
                # quantum_measurement projects only the magic part, so the
                # data part's eigenvalue under the current stabilizer group
                # (-1 for e.g. a -Z input row) must multiply the outcome.
                result ⊻= data_part_eigenvalue(state, op, num_qubits)
                paulistring=embed(size(state.stabilizer_group)[2], op.qubits, op.pauli)
                # The post-measurement state is stabilized by the SIGNED observable:
                # outcome -1 (result=true) stabilizes -P, not +P. Recording +P
                # unconditionally flips later dependent outcomes on half the shots.
                a_stabilizer= Stabilizer([(-1)^result * paulistring])
                check_list=vcat(check_list,a_stabilizer)
                @reset state.stabilizer_group = check_list
                @reset state.runtime = rt
                # Store the full-width Pauli so downstream absolute-index slicing
                # (to_result) sees the correct qubit positions
                return (QuantumRes(paulistring, result), state)
            end
        else
            result = Bool(projection[3]>>1)
            @debug "This measurement outputs Classical Deterministic Result" _group=:api
            # Store the full-width Pauli so all result types share absolute qubit positions
            return (ClassicalDetermRes(embed(size(state.stabilizer_group)[2], op.qubits, op.pauli), result), state)
        end
    end
end

function get_measurement_result(state::CompilerState{TraversalRuntime}, op::CircuitOp.Type)
    @debug "Measuring" op _group=:api
    rt = state.runtime
    outcome_probs = [1-rt.p1_outcome_probs, rt.p1_outcome_probs]
    result = wsample([false,true],outcome_probs)
    check_list=state.stabilizer_group
    num_qubits = get_circuit_width(state.circuit)
    len=length(check_list)
    projection = check_PPM(check_list, op, num_qubits)
    if projection === nothing
        return nothing
    else
        if projection[3] === nothing
            if projection[2]<=len
                @debug "This measurement outputs Classical Random Result" _group=:api
                q_1=[1:get_circuit_width(state.circuit);]
                Q_1=ExpQuatPiPauli(check_list[projection[2]],q_1)
                p_2=(-1)^result*op.pauli
                Q_2=ExpQuatPiPauli(p_2,op.qubits)
                pushfirst!(state.circuit,Q_1,Q_2,Q_1)
                preprocess_circuit(state.circuit)
                # Store the full-width Pauli so all result types share absolute qubit positions
                return (ClassicalRandomRes(embed(size(state.stabilizer_group)[2], op.qubits, op.pauli), result), state)
            else
                paulistring=embed(size(state.stabilizer_group)[2], op.qubits, op.pauli)
                # See the SimRuntime method: the recorded stabilizer row must carry
                # the measured sign
                a_stabilizer= Stabilizer([(-1)^result * paulistring])
                check_list=vcat(check_list,a_stabilizer)
                @reset state.stabilizer_group = check_list
                # Store the full-width Pauli so downstream absolute-index slicing
                # (to_result) sees the correct qubit positions
                return (QuantumRes(paulistring, result), state)
            end
        else
            result = Bool(projection[3]>>1)
            @debug "This measurement outputs Classical Deterministic Result" _group=:api
            # Store the full-width Pauli so all result types share absolute qubit positions
            return (ClassicalDetermRes(embed(size(state.stabilizer_group)[2], op.qubits, op.pauli), result), state)
        end
    end
end
##
"""
    quantum_measurement(state::SimRuntime, op::CircuitOp.Type, num_qubits::Int) -> Tuple{SimRuntime, Bool}
Perform quantum measurement simulation on given state using QuantumClifford.jl backend
"""
function quantum_measurement(rt::SimRuntime, op::CircuitOp.Type, num_qubits::Int)
    quantum_state = rt.quantum_memory
    if quantum_state === nothing
        throw(ArgumentError("Magic State not initiated"))
    end
    magicqubits = collect(num_qubits-length(quantum_state.stab)+1:num_qubits)
    # op.pauli is indexed by position within op.qubits; embed it to the full
    # register before slicing by absolute qubit indices (out-of-bounds
    # PauliOperator indexing silently yields identity instead of throwing)
    real_p=embed(num_qubits, op.qubits, op.pauli)[magicqubits]
    bit_result = projectrand!(quantum_state, real_p)[2]
    result=Bool(bit_result>>1)
    rt = @reset rt.quantum_memory = quantum_state
    return (rt, result)
end

"""
    data_part_eigenvalue(state::CompilerState, op::CircuitOp.Type, num_qubits::Int) -> Bool

Return the eigenvalue bit (`true` denotes -1) of the data-register part of the
joint gadget observable `op` under the current stabilizer group.

A gadget measurement whose data part anticommutes with the group takes the
random branch before reaching the quantum branch, and the group is full rank
over the data register, so here the data part always has a definite sign —
`-1` exactly when the input state carries it (e.g. a `-Z` input row). The
observable's own phase is excluded: it is already carried by the magic-part
slice inside `quantum_measurement`. Runtimes without a magic-state memory
(e.g. `DummyRuntime`) return `false`.
"""
function data_part_eigenvalue(state::CompilerState{SimRuntime}, op::CircuitOp.Type, num_qubits::Int)
    memory = state.runtime.quantum_memory
    memory === nothing && return false
    num_data = num_qubits - length(memory.stab)
    data_p = embed(num_qubits, op.qubits, op.pauli)[1:num_data]
    data_p.phase[] = 0x00
    any(i -> let (x, z) = data_p[i]; x || z end, 1:num_data) || return false
    data_full = embed(num_qubits, collect(1:num_data), data_p)
    projection = project!(copy(state.stabilizer_group), data_full)
    projection[3] === nothing &&
        error("Data part $data_full of gadget measurement $op has no definite eigenvalue under the stabilizer group")
    return Bool(projection[3] >> 1)
end

data_part_eigenvalue(state::CompilerState, op::CircuitOp.Type, num_qubits::Int) = false

"""
    quantum_measurement(state::DummyRuntime, op::CircuitOp.Type, num_qubits::Int) -> Tuple{DummyRuntime, Bool}
Perform quantum measurement simulation on given state using classical sampling according to weight determined by user named outcome_probs
"""
function quantum_measurement(rt::DummyRuntime, op::CircuitOp.Type, num_qubits::Int)
    outcome_probs = [1-rt.p1_outcome_probs, rt.p1_outcome_probs]
    result = wsample([false,true],outcome_probs)
    return (rt,result)
end

"""Resolve conditional circuit operations defined by CircuitOp.BitConditional within circuit defined in state.circuit"""
function resolve_conditionals(state::CompilerState)
    circuit=state.circuit
    creg=state.classical_register
    # Keep resolving until a full scan finds no determined BitConditional.
    # Indices go stale after each splice!/deleteat! (and preprocess_circuit can
    # reorder the circuit), so re-scan from scratch after every resolution.
    resolved = true
    while resolved
        resolved = false
        for i in find_variant_indices(circuit, BitConditional)
            @debug("Start resolving BitConditional at $i")
            operation=circuit[i]
            control_bit=creg[operation.bit]
            if control_bit !== nothing
                if control_bit
                    @debug("$i has a controlled bit")
                    splice!(circuit, i, [operation.op])
                    @debug("$i Resolved")
                else
                    deleteat!(circuit, i)
                    @debug("No correction needed")
                end
                preprocess_circuit(circuit)
                resolved = true
                break
            else
                @debug("Control Bit undetermined")
            end
        end
    end
end
