#This file contains the main logic for compute and compile input circuit into PBC Circuit
using Accessors: @reset
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
    if !isempty(input_circuit)
        validate_circuit(input_circuit)
        if isnothing(input_state)
            input_state=Stabilizer(one(Stabilizer, get_circuit_width(input_circuit); basis=:Z))
        else
            validate_input(input_circuit,input_state)
            input_state=input_state
        end
        shared = build_shared(input_circuit, input_state)
        rt_data = build_rt_data(input_circuit, input_state, rt)
        CompilerState(; shared..., runtime=rt_data)
    else
        CompilerState(MeasurementResult.Type[], S"", Union{Nothing,Bool}[], input_circuit, 1, rt)
    end
end
##
function build_shared(input_circuit::Circuit, input_state::Stabilizer)
    circuit = copy(input_circuit)
    preprocess_circuit(circuit)
    num_bits=get_bit_number(circuit)
    measres=Vector{MeasurementResult.Type}(undef, num_bits)
    creg=Array{Union{Nothing, Bool}}(nothing, num_bits)
    stabgroup=make_stabilizer_list(input_state, circuit)
    return (measurement_results = measres, classical_register = creg, stabilizer_group = stabgroup, circuit = circuit, instruction_pointer = 1)
end
##
function build_rt_data(input_circuit::Circuit, input_state::Stabilizer, rt::SimRuntime)
    num_pauli_qubits=get_circuit_width(input_circuit)
    circuit = copy(input_circuit)
    preprocess_circuit(circuit)
    stabilizier_group=make_stabilizer_list(input_state, circuit)
    num_gadgets=size(stabilizier_group)[2]-num_pauli_qubits
    quantum_memory = num_gadgets==0 ? nothing : create_magic_state(num_gadgets)
    @debug "Number of gadgets inserted" num_gadgets _group=:api
    @reset rt.quantum_memory=quantum_memory
    return rt
end

function build_rt_data(input_circuit::Circuit, input_state::Stabilizer, rt::R) where R<:AbstractRuntime
    return rt
end

"""
    do_quantum_step(state::CompilerState) -> CompilerState

Perform the next joint measurement and update the `CompilerState` accordingly.
"""
function do_quantum_step(state::CompilerState)
    circuit = state.circuit
    i=state.instruction_pointer
    @debug "Now working with $i th measurement" _group=:api
    meas_list = find_variant_indices(circuit, Measurement)
    op=circuit[meas_list[i]]
    bit_index=op.bit
    (meas_result, state)=get_measurement_result(state, op)
    @debug "Measurement result is $(meas_result.result)" _group=:api
    state.measurement_results[i]=meas_result
    state.classical_register[bit_index]=meas_result.result
    @reset state.instruction_pointer = i+1
end

function to_result(state::CompilerState)
    num_qubits = get_circuit_width(state.circuit)
    quantum = filter(mr -> isa_variant(mr, QuantumRes), state.measurement_results)
    magicqubits = collect(num_qubits - length(quantum) + 1 : num_qubits)

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

    CompilationResult(state.measurement_results, qpu_load, state.stabilizer_group, length(qpu_load))
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
    state = build_compilerstate(input_circuit, rt, input_state)
    while !isempty(state.circuit)
        @debug "Working on $(state.instruction_pointer) th PPM" _group=:api
        resolve_conditionals(state)
        # Stop once every measurement has been performed; bounding by both the
        # measurement count and the result-vector length prevents out-of-bounds
        # access when bit indices and measurement counts disagree
        num_meas = length(find_variant_indices(state.circuit, Measurement))
        if state.instruction_pointer > min(num_meas, length(state.measurement_results))
            break
        end
        state=do_quantum_step(state)
        @debug "Performed $(state.instruction_pointer) th PPM" _group=:api
        @debug "Current classical register: $(state.classical_register)" _group=:api
    end
    @debug "Compute/Compile Complete" _group=:api
    return state
end
