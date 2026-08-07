"""
Helper functions to check the first PPM in circuit, determine MeasurementResultType: ClassicalDetermRes, ClassicalRandomRes, QuantumRes
"""
##
using QuantumClifford: project!, Stabilizer, one, GeneralizedStabilizer, apply!, pcT, projectrand!, nqubits, comm, stabilizerview, phases, UnitaryPauliChannel, invsparsity, MixedDestabilizer, rank
using Moshi.Data: variant_name, isa_variant
using Accessors: @reset
##

"""
    validate_input(circuit::Circuit, input::Stabilizer)

Check that a user-supplied input state can be used with `circuit`: it must fit
within the circuit's qubit register, and it must be fully stabilized — one
independent generator per qubit.

Throws an `ArgumentError` describing the violation; returns `nothing` otherwise.
"""
function validate_input(circuit::Circuit, input::Stabilizer)
    num_qubits = nqubits(input)
    if get_circuit_width(circuit) < num_qubits
        throw(ArgumentError("Input state has more qubits than circuit input"))
    end
    # Gadget measurements factor into a data part and a magic part, and
    # `data_part_eigenvalue` requires the data part to have a definite eigenvalue
    # under the stabilizer group — which holds only at full rank. Reject a
    # mixed/underdetermined state here, where the message can name the cause,
    # rather than deep inside a gadget measurement. `rank` (not the row count)
    # also catches generators that are present but linearly dependent.
    state_rank = rank(MixedDestabilizer(input))
    if state_rank != num_qubits
        throw(ArgumentError(
            "Input state must be fully stabilized: it has rank $state_rank over " *
            "$num_qubits qubits ($(length(input)) generators given). Mixed or " *
            "underdetermined input states are not supported by the compilation pipeline."))
    end
    return nothing
end

function create_hadamard_basis_state(num_qubit::Int)
    n = num_qubit

    generators = one(Stabilizer, n; basis=:X)

    return Stabilizer(generators)
end

function create_magic_state(num_magic::Int)
    # The T gates are NOT applied here: they are deferred to
    # `quantum_measurement`, which activates each magic qubit right before the
    # first measurement touching it. This keeps the chi-expansion of the
    # GeneralizedStabilizer at 4^(live qubits) instead of 4^(total), which is
    # exact because T_i commutes with every op not touching qubit i.
    return GeneralizedStabilizer(create_hadamard_basis_state(num_magic))
end

const EmbeddedPcT = typeof(UnitaryPauliChannel(map(p -> embed(1, 1, p), pcT.paulis), pcT.weights))

# One entry per (register width, magic qubit). The channel depends only on those
# two numbers, so rebuilding it on every T activation -- and on every shot of a
# sampling run -- was pure repeat work. Guarded by a lock: the cache outlives any
# single `run`, and nothing else stops two threads from sampling concurrently.
const EMBEDDED_PCT_CACHE = Dict{Tuple{Int,Int}, EmbeddedPcT}()
const EMBEDDED_PCT_LOCK = ReentrantLock()

"""Embed the single-qubit `pcT` channel on qubit `i` of an `n`-qubit register."""
function embedded_pcT(n::Int, i::Int)
    @lock EMBEDDED_PCT_LOCK get!(EMBEDDED_PCT_CACHE, (n, i)) do
        UnitaryPauliChannel(map(p -> embed(n, i, p), pcT.paulis), pcT.weights)
    end
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
    isa_variant(op, CircuitOp.Measurement) || return nothing
    rt = state.runtime
    md = state.stabilizer_group
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
        return (ClassicalDetermRes(pauli, result), state)
    end
end

function get_measurement_result(state::CompilerState{<:StabilizerRuntime}, op::CircuitOp.Type)
    isa_variant(op, CircuitOp.Measurement) || return nothing
    rt = state.runtime
    md = state.stabilizer_group
    num_qubits = nqubits(md)
    pauli = embed(num_qubits, op.qubits, op.pauli)
    projection = project!(md, pauli)
    if projection[3] === nothing
        (rt, result) = quantum_measurement(rt, op, num_qubits)
        phs = phases(stabilizerview(md))
        phs[projection[2]] = (phs[projection[2]] + (result ? 0x2 : 0x0)) & 0x3
        @reset state.runtime = rt
        return (QuantumRes(pauli, result), state)
    else
        # Determined outcome: `project!` reports anticom index 0 here
        result = Bool(projection[3]>>1)
        return (ClassicalDetermRes(pauli, result), state)
    end
end

function get_measurement_result(state::CompilerState{TraversalRuntime}, op::CircuitOp.Type)
    isa_variant(op, CircuitOp.Measurement) || return nothing
    result = rand() < state.runtime.p1_outcome_probs
    md = state.stabilizer_group
    num_qubits = nqubits(md)
    sv = stabilizerview(md)
    pauli = embed(num_qubits, op.qubits, op.pauli)
    # See the SimRuntime method for the branch structure
    anticom = findfirst(i -> comm(pauli, sv, i) != 0x0, 1:length(sv))
    if anticom !== nothing
        q_1=[1:num_qubits;]
        # copy: the row is a view into the tableau, which later projections
        # mutate in place
        Q_1=ExpQuatPiPauli(copy(sv[anticom]),q_1)
        p_2=(-1)^result*op.pauli
        Q_2=ExpQuatPiPauli(p_2,op.qubits)
        pushfirst!(state.circuit,Q_1,Q_2,Q_1)
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
        return (ClassicalDetermRes(pauli, result), state)
    end
end
##
"""
    quantum_measurement(state::SimRuntime, op::CircuitOp.Type, num_qubits::Int) -> Tuple{SimRuntime, Bool}
Perform quantum measurement simulation on given state using QuantumClifford.jl backend
"""
function quantum_measurement(rt::S, op::CircuitOp.Type, num_qubits::Int) where S<:AbstractRuntime
    quantum_state = rt.quantum_memory
    if quantum_state === nothing
        throw(ArgumentError("Magic State not initiated"))
    end
    magicqubits = num_qubits-length(quantum_state.stab)+1:num_qubits
    # op.pauli is indexed by position within op.qubits; embed it to the full
    # register before slicing by absolute qubit indices (out-of-bounds
    # PauliOperator indexing silently yields identity instead of throwing)
    real_p=embed(num_qubits, op.qubits, op.pauli)[magicqubits]
    # Deferred T gates: activate every magic qubit this measurement touches
    # whose T has not been applied yet (see `create_magic_state`). The
    # activation bit stays set after the qubit collapses back to a stabilizer
    # state — reapplying T there would be wrong
    num_magic = length(real_p)
    for k in 1:num_magic
        (x, z) = real_p[k]
        if (x || z) && !rt.activated[k]
            apply!(quantum_state, embedded_pcT(num_magic, k))
            rt.activated[k] = true
        end
    end
    bit_result = projectrand!(quantum_state, real_p)[2]
    result=Bool(bit_result>>1)
    append!(rt.invsparsity_history,invsparsity(quantum_state))
    rt = @reset rt.quantum_memory = quantum_state
    return (rt, result)
end

function quantum_measurement(rt::StabilizerRuntime, op::CircuitOp.Type, num_qubits::Int)
    num_input_qubits = num_qubits - length(rt.activated)
    quantum_state = rt.quantum_memory
    if quantum_state === nothing
        throw(ArgumentError("Magic State not initiated"))
    end
    real_p=embed(num_qubits, op.qubits, op.pauli)
    # Deferred T gates: activate every magic qubit this measurement touches whose
    # T has not been applied yet (see `create_magic_state`). Select those qubits
    # by value: `op.qubits` is the operation's support, so slicing it at the
    # register-partition offset only lands on the magic qubits when the support
    # happens to be the contiguous run 1:num_qubits — which is the case after a
    # rotation has been conjugated to full width, but not for one that reached
    # `gadgetize` unconjugated. There the slice comes out empty and the T is
    # silently skipped.
    for k in op.qubits
        k > num_input_qubits || continue
        (x, z) = real_p[k]
        if (x || z) && !rt.activated[k-num_input_qubits]
            apply!(quantum_state, embedded_pcT(num_qubits, k))
            rt.activated[k-num_input_qubits] = true
        end
    end
    bit_result = projectrand!(quantum_state, real_p)[2]
    result=Bool(bit_result>>1)
    append!(rt.invsparsity_history,invsparsity(quantum_state))
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
function data_part_eigenvalue(state::CompilerState, op::CircuitOp.Type, num_qubits::Int)
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

data_part_eigenvalue(state::CompilerState{DummyRuntime}, op::CircuitOp.Type, num_qubits::Int) = false

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
    while true
        # Find the first BitConditional whose control bit is known. The previous
        # form built the full index vector of every BitConditional and then broke
        # out at the first resolvable one, so all but one entry was discarded
        i = findfirst(op -> isa_variant(op, BitConditional) && creg[op.bit] !== nothing, circuit)
        i === nothing && return
        operation = circuit[i]
        control_bit = creg[operation.bit]
        inner = operation.op
        @debug "Resolving BitConditional at $i" control_bit _group=:api
        if control_bit
            circuit[i] = inner
        else
            deleteat!(circuit, i)
        end
        # A resolved Clifford rotation (the only kind gadgetize emits)
        # just needs commuting past the remaining measurements; other
        # inner ops fall back to the full pipeline
        if !control_bit || isa_variant(inner, ExpHalfPiPauli) || isa_variant(inner, ExpQuatPiPauli)
            absorb_cliffords!(circuit)
        else
            preprocess_circuit(circuit)
        end
    end
end
