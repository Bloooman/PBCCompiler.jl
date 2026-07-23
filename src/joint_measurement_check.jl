"""
Helper functions to check the first PPM in circuit, determine MeasurementResultType: ClassicalDetermRes, ClassicalRandomRes, QuantumRes
"""
##
using QuantumClifford: project!, Stabilizer, one, GeneralizedStabilizer, tensor_pow, apply!, pcT, projectrand!, nqubits, comm, stabilizerview, phases
using Moshi.Data: variant_name, isa_variant
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
    get_measurement_result(state::CompilerState, op::CircuitOp.Type) -> Union{Tuple{MeasurementResult.Type, CompilerState}, Nothing}

Perform Joint Measurement on a CircuitOp.Measurement

Returns `nothing` if `op` is not a `Measurement`; otherwise returns
- a MeasurementResult with corresponding MeasurementResultType (true denotes 0, false denotes 1)
- the updated `CompilerState` (its tableau reflects the post-measurement state;
  for a coin-flip outcome the update happens through compensating rotations
  spliced into the circuit instead)
"""
function get_measurement_result(state::CompilerState, op::CircuitOp.Type)
    @debug "Measuring" op _group=:api
    isa_variant(op, CircuitOp.Measurement) || return nothing
    rt = state.runtime
    md = state.stabilizer_group
    # The tableau spans the full register, so its width is the circuit width
    # and is available in O(1)
    num_qubits = nqubits(md)
    sv = stabilizerview(md)
    # Full-width Pauli: all result types share absolute qubit positions
    pauli = embed(num_qubits, op.qubits, op.pauli)
    # An anticommuting stabilizer row means a coin-flip outcome. Detect it by a
    # row scan so the tableau stays untouched: this branch updates the state by
    # splicing compensating rotations into the circuit, not by projecting.
    anticom = findfirst(i -> comm(pauli, sv, i) != 0x0, 1:length(sv))
    if anticom !== nothing
        result = rand(Bool)
        @debug "This measurement outputs Classical Random Result" _group=:api
        q_1=[1:num_qubits;]
        # copy: the row is a view into the tableau, which later projections
        # mutate in place
        Q_1=ExpQuatPiPauli(copy(sv[anticom]),q_1)
        p_2=(-1)^result*op.pauli
        Q_2=ExpQuatPiPauli(p_2,op.qubits)
        pushfirst!(state.circuit,Q_1,Q_2,Q_1)
        # The inserted rotations are Clifford, so commuting them out is all
        # the preprocessing pipeline would do here
        absorb_cliffords!(state.circuit)
        return (ClassicalRandomRes(pauli, result), state)
    end
    projection = project!(md, pauli)
    if projection[3] === nothing
        # Commuting and independent of the group: project! grew the rank in
        # place, adding `pauli` as the row at index projection[2]
        (rt, result) = quantum_measurement(rt, op, num_qubits)
        # The joint observable factors as (data part) ⊗ (magic part).
        # quantum_measurement projects only the magic part, so the
        # data part's eigenvalue under the current stabilizer group
        # (-1 for e.g. a -Z input row) must multiply the outcome.
        result ⊻= data_part_eigenvalue(state, op, num_qubits)
        # The post-measurement state is stabilized by the SIGNED observable:
        # outcome -1 (result=true) stabilizes -P, not +P. Recording +P
        # unconditionally flips later dependent outcomes on half the shots.
        phs = phases(stabilizerview(md))
        phs[projection[2]] = (phs[projection[2]] + (result ? 0x2 : 0x0)) & 0x3
        @reset state.runtime = rt
        return (QuantumRes(pauli, result), state)
    else
        result = Bool(projection[3]>>1)
        @debug "This measurement outputs Classical Deterministic Result" _group=:api
        return (ClassicalDetermRes(pauli, result), state)
    end
end

function get_measurement_result(state::CompilerState{TraversalRuntime}, op::CircuitOp.Type)
    @debug "Measuring" op _group=:api
    isa_variant(op, CircuitOp.Measurement) || return nothing
    result = rand() < state.runtime.p1_outcome_probs
    md = state.stabilizer_group
    # The tableau spans the full register, so its width is the circuit width
    # and is available in O(1)
    num_qubits = nqubits(md)
    sv = stabilizerview(md)
    # Full-width Pauli: all result types share absolute qubit positions
    pauli = embed(num_qubits, op.qubits, op.pauli)
    # See the SimRuntime method for the branch structure
    anticom = findfirst(i -> comm(pauli, sv, i) != 0x0, 1:length(sv))
    if anticom !== nothing
        @debug "This measurement outputs Classical Random Result" _group=:api
        q_1=[1:num_qubits;]
        # copy: the row is a view into the tableau, which later projections
        # mutate in place
        Q_1=ExpQuatPiPauli(copy(sv[anticom]),q_1)
        p_2=(-1)^result*op.pauli
        Q_2=ExpQuatPiPauli(p_2,op.qubits)
        pushfirst!(state.circuit,Q_1,Q_2,Q_1)
        # The inserted rotations are Clifford, so commuting them out is all
        # the preprocessing pipeline would do here
        absorb_cliffords!(state.circuit)
        return (ClassicalRandomRes(pauli, result), state)
    end
    projection = project!(md, pauli)
    if projection[3] === nothing
        # The recorded stabilizer row must carry the measured sign
        phs = phases(stabilizerview(md))
        phs[projection[2]] = (phs[projection[2]] + (result ? 0x2 : 0x0)) & 0x3
        return (QuantumRes(pauli, result), state)
    else
        result = Bool(projection[3]>>1)
        @debug "This measurement outputs Classical Deterministic Result" _group=:api
        return (ClassicalDetermRes(pauli, result), state)
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
    magicqubits = num_qubits-length(quantum_state.stab)+1:num_qubits
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
    result = rand() < rt.p1_outcome_probs
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
                inner = operation.op
                if control_bit
                    @debug("$i has a controlled bit")
                    splice!(circuit, i, [inner])
                    @debug("$i Resolved")
                else
                    deleteat!(circuit, i)
                    @debug("No correction needed")
                end
                # A resolved Clifford rotation (the only kind gadgetize emits)
                # just needs commuting past the remaining measurements; other
                # inner ops fall back to the full pipeline
                if !control_bit || isa_variant(inner, ExpHalfPiPauli) || isa_variant(inner, ExpQuatPiPauli)
                    absorb_cliffords!(circuit)
                else
                    preprocess_circuit(circuit)
                end
                resolved = true
                break
            else
                @debug("Control Bit undetermined")
            end
        end
    end
end
