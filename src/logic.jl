#This file contains the main logic for compute and compile input circuit into PBC Circuit
using Accessors: @reset
"""
    preprocess_circuit(circuit::Circuit) -> Circuit

Process a Clifford + T circuit into a Pauli Product circuit by appropriately commuting
all Clifford gates past the nonClifford-gates and absorbing them in the Pauli Product Measurements.
Then replace all nonClifford Pauli Product Rotations with gadgets.
"""
function preprocess_circuit(circuit::Circuit)
    if isempty(circuit) || length(circuit) < 2
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
    get_state(circuit::Circuit, input_state::Union{Stabilizer, Nothing}=nothing, dummy::Bool=false) -> S where S <: AbstractRuntime

Get initial Execution State using input circuit and input state.

# Fields
- `circuit`: input circuit for compilation
- `input_state`: initial qubit state
- `dummy`: flag for running dummy simulation
- `outcome_probs`: 2-element distribution vector [p_p1, p_m1]. p_p1 is the probability measuring +1; p_m1 is the probability measuring -1
"""
function get_state(circuit::Circuit, input_state::Union{Stabilizer, Nothing}=nothing; dummy::Bool=false, outcome_probs::Vector{Int}=[1,1])
    validate_circuit(circuit)
    num_pauli_qubits=get_circuit_width(circuit)
    @debug "Initial number of qubits $num_pauli_qubits" _group=:api
    if isnothing(input_state)
        state=Stabilizer(one(Stabilizer, num_pauli_qubits; basis=:Z))
    else
        validate_input(circuit,input_state)
        state=input_state
    end
    preprocess_circuit(circuit)
    stabilzier_group=make_stabilizer_list(state, circuit)
    num_gadgets=size(stabilzier_group)[2]-num_pauli_qubits
    @debug "Number of gadgets inserted" num_gadgets _group=:api
    magicstate = num_gadgets==0 || dummy ? nothing : create_magic_state(num_gadgets)
    num_bits=get_bit_number(circuit)
    MeasRes=Vector{MeasurementResult.Type}(undef, num_bits)
    creg=Array{Union{Nothing, Bool}}(nothing, num_bits)
    cs=CompilerState(MeasRes, stabilzier_group, creg, circuit, num_gadgets, 1)
    if dummy
        return DummyRuntime(cs, outcome_probs)
    else
        return SimRuntime(cs, magicstate)
    end
end

"""
    do_quantum_step(state::S) -> S where S <: AbstractSimState

Perform the next joint measurement and update ComputerState accordingly
"""
function do_quantum_step(state::S) where S <: AbstractRuntime
    cs = state.compiler_state
    circuit = cs.circuit
    i=cs.instruction_pointer
    @debug "Now working with $i th measurement" _group=:api
    meas_list = find_variant_indices(circuit,Measurement)
    op=circuit[meas_list[i]]
    bit_index=op.bit
    (meas_result,state)=get_measurement_result(state, op)
    @debug "Measurement result is $res" _group=:api
    state.compiler_state.measurement_results[i]=meas_result
    state.compiler_state.classical_register[bit_index]=meas_result.result
    @reset state.compiler_state.instruction_pointer = i+1
end

"""
    run(input_circuit::Circuit, input_state::Union{Stabilizer, Nothing}=nothing; dummy::Bool=false, outcome_probs::Vector{Int}=[1,1]) -> S where S <: AbstractRuntime

Run compute/compile with provided circuit and input state(described by stabilizer group)
# Fields
- `circuit`: input circuit for compilation
- `input_state`: initial qubit state
- `dummy`: flag for running dummy simulation
- `outcome_probs`: 2-element distribution vector [p_p1, p_m1]. p_p1 is the probability measuring +1; p_m1 is the probability measuring -1
"""
function run(input_circuit::Circuit, input_state::Union{Stabilizer, Nothing}=nothing; dummy::Bool=false, outcome_probs::Vector{Int}=[1,1])
    state = get_state(input_circuit, input_state; dummy=dummy, outcome_probs=outcome_probs)
    len=length(state.compiler_state.classical_register)
    while true && !isempty(state.compiler_state.circuit)
        @debug "Working on $(state.compiler_state.instruction_pointer) th PPM" _group=:api
        resolve_conditionals(state)
        state=do_quantum_step(state)
        @debug "Performed $(state.compiler_state.instruction_pointer) th PPM" _group=:api
        @debug "Current classical register: $(state.compiler_state.classical_register)" _group=:api
        if state.compiler_state.instruction_pointer>len
            break
        end
    end
    @debug "Compute/Compile Complete" _group=:api
    pbc_circuit=[]
    for i in 1:length(state.compiler_state.circuit)
        op=state.compiler_state.circuit[i]
        pauli=state.compiler_state.measurement_results[i].pauli
        new_op=CircuitOp.Measurement(pauli,op.bit,op.qubits)
        push!(pbc_circuit,new_op)
    end
    return @reset state.compiler_state.circuit = pbc_circuit
    @debug "Result returned" _group=:api
end
