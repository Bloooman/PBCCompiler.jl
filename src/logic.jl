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
Build the runtime data from an already-preprocessed circuit. For `SimRuntime`,
allocate one magic-state qubit per gadget; gadgetization appends its magic
qubits above the input width, so the gadget count is the width difference.
"""
function build_rt_data(preprocessed::Circuit, input_state::Stabilizer, num_input_qubits::Int, rt::SimRuntime)
    num_gadgets = get_circuit_width(preprocessed) - num_input_qubits
    quantum_memory = num_gadgets==0 ? nothing : create_magic_state(num_gadgets)
    @debug "Number of gadgets inserted" num_gadgets _group=:api
    @reset rt.quantum_memory=quantum_memory
    @reset rt.activated = num_gadgets==0 ? nothing : falses(num_gadgets)
    # Fresh vector, not `empty!`: callers reuse one runtime across shots, and
    # appending into the caller's own vector would concatenate every shot's
    # telemetry into one series
    @reset rt.invsparsity_history = Int[]
    return rt
end

function build_rt_data(preprocessed::Circuit, input_state::Stabilizer, num_input_qubits::Int, rt::StabilizerRuntime)
    num_qubits = get_circuit_width(preprocessed)
    num_gadgets = num_qubits - num_input_qubits
    input = make_stabilizer_list(input_state, preprocessed)
    generators = one(Stabilizer, num_qubits; basis=:X)
    state = Stabilizer(generators)
    state[1:num_input_qubits] = input
    quantum_memory = GeneralizedStabilizer(state)
    @debug "Number of gadgets inserted" num_gadgets _group=:api
    @reset rt.quantum_memory=quantum_memory
    @reset rt.activated = num_gadgets==0 ? nothing : falses(num_gadgets)
    # See the SimRuntime method: fresh vector so shots do not concatenate
    @reset rt.invsparsity_history = Int[]
    return rt
end

function build_rt_data(preprocessed::Circuit, input_state::Stabilizer, num_input_qubits::Int, rt::R) where R<:AbstractRuntime
    return rt
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

function to_result(state::CompilerState)
    num_qubits = nqubits(state.stabilizer_group)
    quantum = filter(mr -> isa_variant(mr, QuantumRes), state.measurement_results)
    magicqubits = num_qubits - length(quantum) + 1 : num_qubits

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
            push!(qpu_load, QuantumRes(magic_p, mr.result))
        end
    end

    # CompilationResult keeps the plain Stabilizer representation (stable
    # serialization format); extract it from the working tableau
    CompilationResult(state.measurement_results, qpu_load, copy(stabilizerview(state.stabilizer_group)), length(qpu_load))
end

function to_result(state::CompilerState{<:StabilizerRuntime})
    num_qubits = nqubits(state.stabilizer_group)
    quantum = filter(mr -> isa_variant(mr, QuantumRes), state.measurement_results)
    num_input_qubits = num_qubits - length(state.runtime.activated)
    magicqubits = num_input_qubits + 1 : num_qubits

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
            push!(qpu_load, QuantumRes(magic_p, mr.result))
        end
    end

    keep_qubits = collect(1:num_input_qubits)
    s_sub = stabilizerview(state.stabilizer_group)[:, keep_qubits]
    non_trivial_rows = [i for i in 1:length(s_sub) if !iszero(s_sub[i].xz)]
    s_clean = s_sub[non_trivial_rows]

    # CompilationResult keeps the plain Stabilizer representation (stable
    # serialization format); extract it from the working tableau
    CompilationResult(state.measurement_results, qpu_load, s_clean, length(qpu_load))
end

"""
    run(input_circuit::Circuit, rt::AbstractRuntime, input_state::Union{Stabilizer, Nothing}=nothing) -> CompilerState

Compute/compile the provided circuit against an input state (described by a stabilizer
group) and return the final `CompilerState`.

This function is not exported because it shadows `Base.run`; call it as
`PBCCompiler.run`.

# Arguments
- `input_circuit`: input circuit for compilation
- `rt`: runtime that supplies measurement outcomes (`SimRuntime`, `DummyRuntime`, `TraversalRuntime`)
- `input_state`: initial qubit state; defaults to the all-|0⟩ state
"""
function run(input_circuit::Circuit, rt::S, input_state::Union{Stabilizer, Nothing}=nothing) where S <: AbstractRuntime
    execute!(build_compilerstate(input_circuit, rt, input_state))
end

"""
    execute!(state::CompilerState) -> CompilerState

Perform every remaining measurement of an already-compiled `CompilerState` and
return the final state.

`state` is consumed: its tableau, circuit and registers are updated in place.
To run several shots off one compilation, pass a `copy` of the compiled state
on each shot.

# Examples
```julia
compiled = build_compilerstate(circuit, SimRuntime())
shots = [execute!(copy(compiled)) for _ in 1:100]
```
"""
function execute!(state::CompilerState)
    while !isempty(state.circuit)
        resolve_conditionals(state)
        # Stop once every measurement has been performed; bounding by both the
        # measurement count and the result-vector length prevents out-of-bounds
        # access when bit indices and measurement counts disagree
        meas_list = find_variant_indices(state.circuit, Measurement)
        if state.instruction_pointer > min(length(meas_list), length(state.measurement_results))
            break
        end
        state=do_quantum_step(state, meas_list)
    end
    return state
end
