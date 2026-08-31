using Accessors: @reset
using QuantumClifford: MixedDestabilizer, Tableau, nqubits, stabilizerview
"""
    preprocess_circuit(circuit::Circuit) -> Circuit

Process a Clifford + T circuit into a Pauli Product circuit by appropriately commuting
all Clifford gates past the nonClifford-gates and absorbing them in the Pauli Product Measurements.
Then replace all nonClifford Pauli Product Rotations with gadgets.
"""
function preprocess_circuit(circuit::Circuit)
    if isempty(circuit)
        return circuit
    end
    remove_pauliconditional(circuit)
    group_nonclifford(circuit)
    merge_ops(circuit)
    remove_clifford(circuit)
    remove_nonclifford(circuit)
    remove_post_measurement(circuit)
end

"""
    build_compilerstate(input_circuit::Circuit, rt::S, input_state::Union{Stabilizer, Nothing}=nothing) where S <: AbstractRuntime

Build the initial `CompilerState` from an input circuit and input state.

# Arguments
- `input_circuit`: input circuit for compilation
- `rt`: runtime for compilation
- `input_state`: initial qubit state; defaults to the all-|0⟩ state
"""
function build_compilerstate(input_circuit::Circuit, rt::S, input_state::Union{Stabilizer, Nothing}=nothing) where S <: AbstractRuntime
    isempty(input_circuit) && return _empty_state(input_circuit, rt)
    validate_circuit(input_circuit)
    if isnothing(input_state)
        input_state=Stabilizer(one(Stabilizer, get_circuit_width(input_circuit); basis=:Z))
    else
        validate_input(input_circuit,input_state)
    end
    num_input_qubits = get_circuit_width(input_circuit)
    circuit = copy(input_circuit)
    # The one and only preprocessing run: both the shared execution state and
    # the runtime's magic-state memory are derived from this same circuit, so
    # they cannot desync
    preprocess_circuit(circuit)
    # A circuit can cancel to nothing (e.g. a rotation followed by its inverse)
    isempty(circuit) && return _empty_state(circuit, rt)
    shared = build_shared(circuit, input_state)
    rt_data = build_rt_data(circuit, input_state, num_input_qubits, rt)
    CompilerState(; shared..., runtime=rt_data)
end

"""Compiler state for a circuit with nothing to execute."""
function _empty_state(circuit::Circuit, rt::AbstractRuntime)
    # MixedDestabilizer(::Stabilizer) rejects zero-qubit input, so build the
    # empty tableau explicitly
    empty_group = MixedDestabilizer(Tableau(zeros(UInt8, 0), 0, zeros(UInt64, 0, 0)), 0)
    CompilerState(MeasurementResult.Type[], empty_group, Union{Nothing,Bool}[], circuit, 1, rt)
end
##
"""Build the runtime-independent parts of the `CompilerState` from an already-preprocessed circuit."""
function build_shared(preprocessed::Circuit, input_state::Stabilizer)
    num_bits=get_bit_number(preprocessed)
    measres=Vector{MeasurementResult.Type}(undef, num_bits)
    creg=Array{Union{Nothing, Bool}}(nothing, num_bits)
    # The destabilizer half lets each measurement project in O(n^2) instead of
    # recanonicalizing the whole group
    stabgroup=MixedDestabilizer(make_stabilizer_list(input_state, preprocessed))
    return (measurement_results = measres, classical_register = creg, stabilizer_group = stabgroup, circuit = preprocessed, instruction_pointer = 1)
end
##
"""
Number of gadget (magic-state) qubits a preprocessed circuit needs: gadgetization
appends its magic qubits above the input width, so the gadget count is the
width difference.
"""
_gadget_count(preprocessed::Circuit, num_input_qubits::Int) = get_circuit_width(preprocessed) - num_input_qubits

"""
Build the runtime data from an already-preprocessed circuit. For `SimRuntime`,
allocate one magic-state qubit per gadget; gadgetization appends its magic
qubits above the input width, so the gadget count is the width difference.
"""
function build_rt_data(preprocessed::Circuit, input_state::Stabilizer, num_input_qubits::Int, rt::Union{SimRuntime, StabilizerRuntime, HybridRuntime})
    num_qubits = get_circuit_width(preprocessed)
    num_gadgets = _gadget_count(preprocessed, num_input_qubits)
    input = make_stabilizer_list(input_state, preprocessed)
    generators = one(Stabilizer, num_qubits; basis=:X)
    state = Stabilizer(generators)
    state[1:num_input_qubits] = input
    quantum_memory = GeneralizedStabilizer(state)
    @debug "Number of gadgets inserted" num_gadgets _group=:api
    @reset rt.quantum_memory=quantum_memory
    # Unlike SimRuntime -- which allocates no magic register at all when there
    # are no gadgets, so `nothing` is the honest value there -- this runtime
    # always holds the full register. An empty BitVector keeps
    # `num_input_qubits = num_qubits - length(activated)` correct on gadget-free
    # circuits instead of throwing `MethodError: length(::Nothing)`
    @reset rt.activated = falses(num_gadgets)
    # See the SimRuntime method: fresh vector so shots do not concatenate
    @reset rt.invsparsity_history = Int[]
    return rt
end

function build_rt_data(preprocessed::Circuit, input_state::Stabilizer, num_input_qubits::Int, rt::collapseRuntime)
    num_qubits = get_circuit_width(preprocessed)
    num_gadgets = _gadget_count(preprocessed, num_input_qubits)
    input = make_stabilizer_list(input_state, preprocessed)
    generators = one(Stabilizer, num_qubits; basis=:X)
    state = Stabilizer(generators)
    state[1:num_input_qubits] = input
    quantum_memory = GeneralizedStabilizer(state)
    @debug "Number of gadgets inserted" num_gadgets _group=:api
    @reset rt.quantum_memory=quantum_memory
    # Unlike SimRuntime -- which allocates no magic register at all when there
    # are no gadgets, so `nothing` is the honest value there -- this runtime
    # always holds the full register. An empty BitVector keeps
    # `num_input_qubits = num_qubits - length(activated)` correct on gadget-free
    # circuits instead of throwing `MethodError: length(::Nothing)`
    @reset rt.activated = falses(num_gadgets)
    @reset rt.collapsed = falses(num_gadgets)
    # See the SimRuntime method: fresh vector so shots do not concatenate
    @reset rt.invsparsity_history = Int[]
    return rt
end

function build_rt_data(preprocessed::Circuit, input_state::Stabilizer, num_input_qubits::Int, rt::DummyRuntime)
    num_gadgets = _gadget_count(preprocessed, num_input_qubits)
    @reset rt.activated = num_gadgets==0 ? nothing : falses(num_gadgets)
    return rt
end

function build_rt_data(preprocessed::Circuit, input_state::Stabilizer, num_input_qubits::Int, rt::DummyStabilizerRuntime)
    num_gadgets = _gadget_count(preprocessed, num_input_qubits)
    @reset rt.activated = falses(num_gadgets)
    return rt
end

function build_rt_data(preprocessed::Circuit, input_state::Stabilizer, num_input_qubits::Int, rt::R) where R<:AbstractRuntime
    return rt
end

function transition(state::CompilerState{<:HybridRuntime})
     rt = state.runtime
     max = rt.maximum_measurement_support
     (isnothing(max) || isnothing(rt.activated)) && return state
     num_input_qubits = nqubits(state.stabilizer_group) - num_gadget_qubits(rt)
     count(rt.activated) < max - num_input_qubits && return state
     n_measurements = count(i -> isassigned(state.measurement_results, i), 1:length(state.measurement_results))
     @reset state.runtime = HybridStabilizerRuntime(rt.quantum_memory, rt.activated, rt.invsparsity_history, copy(rt.activated), n_measurements)
     return state
end

function transition(state::CompilerState{<:AbstractStabilizerRuntime})
    return state
end

"""
    do_quantum_step(state::CompilerState) -> CompilerState

Perform the next joint measurement and update the `CompilerState` accordingly.

`meas_list` are the positions of the `Measurement` ops in `state.circuit`;
callers that already computed it (the `run` loop) pass it in to avoid a second
scan of the circuit.
"""
function do_quantum_step(state::CompilerState, meas_list::Vector{Int}=find_variant_indices(state.circuit, Measurement))
    circuit = state.circuit
    i=state.instruction_pointer
    op=circuit[meas_list[i]]
    bit_index=op.bit
    (meas_result, state)=get_measurement_result(state, op)
    @debug "Measurement $i" meas_result.pauli meas_result.result _group=:api
    state.measurement_results[i]=meas_result
    state.classical_register[bit_index]=meas_result.result
    @reset state.instruction_pointer = i+1
end

"""
Restrict each `QuantumRes` in `quantum` to the qubits in `magicqubits`, folding
away any magic qubit whose Pauli support is fully accounted for by an earlier
entry in the list (a joint measurement whose only remaining non-identity
support is one already-resolved magic qubit needs no further QPU operation).

When `embed_width` is given, each kept entry's Pauli is placed back at that
full register width (identity elsewhere) instead of returned restricted to
`magicqubits` alone -- used by `to_result(::CompilerState{<:HybridStabilizerRuntime})`
so every `QPU_workload` entry has the same width as the post-transition ones
it's concatenated with.
"""
function _magic_block_qpu_load(quantum, magicqubits, embed_width::Union{Int,Nothing}=nothing)
    excluded_local = Set{Int}()
    qpu_load = Vector{MeasurementResult.Type}()

    for mr in quantum
        magic_p = mr.pauli[magicqubits]
        for i in excluded_local
            magic_p[i] = (false, false)
        end
        non_id = findall(i -> let (x, z) = magic_p[i]; x || z end, 1:length(magic_p))
        if length(non_id) == 1
            push!(excluded_local, non_id[1])
        else
            p = embed_width === nothing ? magic_p : embed(embed_width, magicqubits, magic_p)
            push!(qpu_load, QuantumRes(p, mr.result))
        end
    end
    return qpu_load
end

function to_result(state::CompilerState)
    num_qubits = nqubits(state.stabilizer_group)
    meas_result = state.measurement_results
    assigned = [meas_result[i] for i in 1:length(meas_result) if isassigned(meas_result, i)]
    quantum = filter(mr -> isa_variant(mr, QuantumRes), assigned)
    magicqubits = num_qubits - length(quantum) + 1 : num_qubits
    qpu_load = _magic_block_qpu_load(quantum, magicqubits)

    # CompilationResult keeps the plain Stabilizer representation (stable
    # serialization format); extract it from the working tableau
    CompilationResult(state.measurement_results, qpu_load, copy(stabilizerview(state.stabilizer_group)), length(qpu_load))
end

"""
Number of magic-state (gadget) qubits a `StabilizerRuntime`/`DummyStabilizerRuntime`/`HybridStabilizerRuntime`
holds -- the trailing block of the register, above the data qubits.

`activated` is `nothing` on a state from [`_empty_state`](@ref), which never runs
`build_rt_data`; a circuit with nothing to execute has no gadgets, hence 0.

`SimRuntime`/`DummyRuntime` do not need this: the generic `to_result` derives
its magic block from the `QuantumRes` count rather than from `activated`.
"""
num_gadget_qubits(rt::AbstractRuntime) =
    rt.activated === nothing ? 0 : length(rt.activated)

"""
`to_result` for `collapseRuntime`. The generic `to_result(state::CompilerState)`
sizes the magic-qubit window from the `QuantumRes` count, which assumes every
gadget measurement lands as `QuantumRes` -- `collapseRuntime` breaks that
assumption on purpose, reclassifying a magic qubit's measurement as
`ClassicalBiasedRes` once its support has collapsed. Sized from
`num_gadget_qubits(state.runtime)` instead (like `AbstractStabilizerRuntime`'s
method) so `magicqubits`/`QPU_workload` stay correctly sized regardless of how
many measurements collapsed.

`ClassicalBiasedRes` entries are excluded from `QPU_workload`, same as the
generic method excludes `ClassicalDetermRes`/`ClassicalRandomRes` -- whether
that's the right semantics needs further discussion.
"""
function to_result(state::CompilerState{<:collapseRuntime})
    num_qubits = nqubits(state.stabilizer_group)
    meas_result = state.measurement_results
    assigned = [meas_result[i] for i in 1:length(meas_result) if isassigned(meas_result, i)]
    quantum = filter(mr -> isa_variant(mr, QuantumRes), assigned)
    num_gadgets = num_gadget_qubits(state.runtime)
    magicqubits = num_qubits - num_gadgets + 1 : num_qubits
    qpu_load = _magic_block_qpu_load(quantum, magicqubits)

    # CompilationResult keeps the plain Stabilizer representation (stable
    # serialization format); extract it from the working tableau
    CompilationResult(state.measurement_results, qpu_load, copy(stabilizerview(state.stabilizer_group)), length(qpu_load))
end

function to_result(state::CompilerState{<:AbstractStabilizerRuntime})
    num_qubits = nqubits(state.stabilizer_group)
    meas_result = state.measurement_results
    assigned = [meas_result[i] for i in 1:length(meas_result) if isassigned(meas_result, i)]
    quantum = filter(mr -> isa_variant(mr, QuantumRes), assigned)
    num_input_qubits = num_qubits - num_gadget_qubits(state.runtime)

    keep_qubits = collect(1:num_input_qubits)
    s_sub = stabilizerview(state.stabilizer_group)[:, keep_qubits]
    non_trivial_rows = [i for i in 1:length(s_sub) if !iszero(s_sub[i].xz)]
    s_clean = s_sub[non_trivial_rows]

    # CompilationResult keeps the plain Stabilizer representation (stable
    # serialization format); extract it from the working tableau
    CompilationResult(state.measurement_results, quantum, s_clean, length(quantum))
end

"""
`to_result` for a `HybridRuntime` that converted mid-run. Measurements taken
before the transition were `SimRuntime`-style: only their magic-qubit support
is real QPU work (the data part was already resolved classically), so those
`QuantumRes` entries go through the same magic-block restriction the generic
`to_result` applies. Measurements taken after the transition were genuine
whole-register `StabilizerRuntime` projections, so they're kept as-is, same
as the plain `AbstractStabilizerRuntime` method. The result tableau keeps the
data qubits plus whichever magic qubits were already live in `quantum_memory`
at the transition point ([`HybridStabilizerRuntime.activated_at_transition`](@ref)) --
magic qubits only touched after the transition are pure `StabilizerRuntime`
territory and don't carry the same "live quantum resource" meaning.

Pre-transition entries are embedded back to full register width (identity on
every qubit their restricted-and-locally-absorbed Pauli doesn't cover)
instead of kept at the narrow magic-block width, so every entry in the
returned `QPU_workload` has the same length as the post-transition ones it's
concatenated with.
"""
function to_result(state::CompilerState{<:HybridStabilizerRuntime})
    rt = state.runtime
    num_qubits = nqubits(state.stabilizer_group)
    meas_result = state.measurement_results
    assigned_idx = [i for i in 1:length(meas_result) if isassigned(meas_result, i)]
    num_input_qubits = num_qubits - num_gadget_qubits(rt)

    pre_idx = filter(i -> i <= rt.n_measurements_at_transition, assigned_idx)
    post_idx = filter(i -> i > rt.n_measurements_at_transition, assigned_idx)
    pre_quantum = filter(mr -> isa_variant(mr, QuantumRes), meas_result[pre_idx])
    post_quantum = filter(mr -> isa_variant(mr, QuantumRes), meas_result[post_idx])

    magicqubits = num_input_qubits+1:num_qubits
    qpu_load = vcat(_magic_block_qpu_load(pre_quantum, magicqubits, num_qubits), post_quantum)

    keep_qubits = vcat(1:num_input_qubits, num_input_qubits .+ findall(rt.activated_at_transition))
    s_sub = stabilizerview(state.stabilizer_group)[:, keep_qubits]
    non_trivial_rows = [i for i in 1:length(s_sub) if !iszero(s_sub[i].xz)]
    s_clean = s_sub[non_trivial_rows]

    # CompilationResult keeps the plain Stabilizer representation (stable
    # serialization format); extract it from the working tableau
    CompilationResult(state.measurement_results, qpu_load, s_clean, length(qpu_load))
end

"""
Whether every measurement in `state.circuit` has already been performed --
the shared stopping condition for [`run`](@ref)/[`execute!`](@ref) and for any
caller driving `execute!` to completion by hand. Bounding by both the
measurement count and the result-vector length prevents out-of-bounds access
when bit indices and measurement counts disagree.
"""
_execution_complete(state::CompilerState) =
    state.instruction_pointer > min(length(find_variant_indices(state.circuit, Measurement)), length(state.measurement_results))

"""
    run(input_circuit::Circuit, rt::AbstractRuntime, input_state::Union{Stabilizer, Nothing}=nothing) -> CompilerState

Compute/compile the provided circuit against an input state (described by a stabilizer
group) and return the final `CompilerState`.

This function is not exported because it shadows `Base.run`; call it as
`PBCCompiler.run`.

# Arguments
- `input_circuit`: input circuit for compilation
- `rt`: runtime that supplies measurement outcomes (`SimRuntime`, `StabilizerRuntime`, `DummyRuntime`, `DummyStabilizerRuntime`, `HybridRuntime`, `TraversalRuntime`)
- `input_state`: initial qubit state; defaults to the all-|0⟩ state
"""
function run(input_circuit::Circuit, rt::S, input_state::Union{Stabilizer, Nothing}=nothing) where S <: AbstractRuntime
    state = build_compilerstate(input_circuit, rt, input_state)
    while !_execution_complete(state)
        state = execute!(state)
    end
    return state
end

function run(input_circuit::Circuit, rt::HybridRuntime, input_state::Union{Stabilizer, Nothing}=nothing)
    state = build_compilerstate(input_circuit, rt, input_state)
    while !_execution_complete(state)
        state = execute!(state)
        state = transition(state)
    end
    return state
end

"""
    execute!(state::CompilerState) -> CompilerState

Perform the next remaining measurement step of an already-compiled
`CompilerState` and return the updated state. A no-op (returns `state`
unchanged) once every measurement has been performed.

`state` is consumed: its tableau, circuit and registers are updated in place.
To run several shots off one compilation, pass a `copy` of the compiled state
on each shot, looping `execute!` until the shot is complete.

# Examples
```julia
compiled = build_compilerstate(circuit, SimRuntime())
shots = map(1:100) do _
    s = copy(compiled)
    while !PBCCompiler._execution_complete(s)
        s = execute!(s)
    end
    s
end
```
"""
function execute!(state::CompilerState)
    resolve_conditionals(state)
    _execution_complete(state) && return state
    meas_list = find_variant_indices(state.circuit, Measurement)
    state=do_quantum_step(state, meas_list)
    return state
end
