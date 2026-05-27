module PBCCompiler

using Moshi.Data: @data, variant_name, isa_variant
using Moshi.Match: @match
using Moshi.Derive: @derive
using QuantumClifford: PauliOperator, @P_str, comm, embed, ⊗, random_pauli, tensor, @S_str, Stabilizer, project!, GeneralizedStabilizer, projectrand!
using Random: randstring
using StatsBase: sample
using Accessors: @reset
##

"""TODO docstring"""
const P = typeof(P"XYZ")

"""TODO docstring"""
@data CircuitOp begin
    """Measurement of pauli string P (ie., + XY) on qubits in vector at field "qubits" (ie.,[1,3]), measurement result is stored in classical bit denoted in "bit" """
    struct Measurement
        pauli::P
        bit::Int
        qubits::Vector{Int}
    end
    """TODO docstring"""
    struct Pauli
        pauli::P
        qubits::Vector{Int}
    end
    """Perform Pauli Product Rotation(PPR) in the form of Pφ = exp(−iP φ), where P is pauli string, φ is an angle Perform pi/2 PPR on qubits denoted in Vector qubits"""
    struct ExpHalfPiPauli
        pauli::P
        qubits::Vector{Int}
    end
    """Perform pi/4 PPR on qubits denoted in Vector qubits"""
    struct ExpQuatPiPauli
        pauli::P
        qubits::Vector{Int}
    end
    """Perform pi/8 PPR on qubits denoted in Vector qubits"""
    struct ExpEighPiPauli
        pauli::P
        qubits::Vector{Int}
    end
    """TODO docstring"""
    struct PrepMagic
        qubit::Int
        qubits::Vector{Int}
    end
    """Perform a (pi/2) Pauli rotation (defined by target_pauli) on the target qubits, conditional on the control qubits falling into the -1 eigenspace of control_pauli"""
    struct PauliConditional
        control_pauli::P
        control_qubits::Vector{Int}
        target_pauli::P
        target_qubits::Vector{Int}
    end
    """TODO docstring"""
    struct BitConditional
        op::CircuitOp
        bit::Int
    end
end

@derive CircuitOp[Hash, Eq, Show]

"""TODO docstring"""
const Circuit = Vector{CircuitOp.Type}

using .CircuitOp: Measurement, Pauli, ExpHalfPiPauli, ExpQuatPiPauli, ExpEighPiPauli, PrepMagic, PauliConditional, BitConditional

include("traversal.jl")
include("affectedqubits.jl")
include("plotting.jl")
include("pair_transformation.jl")
include("preprocess.jl")
include("Random_Circuit.jl")
##

"""TODO docstring"""
function make_counter()
    var = Ref{Int}()
    var[] = 0
    return function counter()
        var[] += 1
        return var[]
    end
end

##

"""
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

"""Reweite P1-controlled-P2 gates as C(P1, P2) = (P1 ⊗ P2)π/4 · (1 ⊗ P2)−π/4 · (P1 ⊗ 1)−π/4."""
function remove_pauliconditional(circuit::Circuit)
    indices=find_variant_indices(circuit,PauliConditional)
    for i in reverse(indices)
        op=circuit[i]
        @match op begin
            CircuitOp.PauliConditional(cp, cq, tp, tq) => begin
                op_1=CircuitOp.ExpQuatPiPauli(-cp, cq)
                op_2=CircuitOp.ExpQuatPiPauli(-tp, tq)
                op_3=CircuitOp.ExpQuatPiPauli(cp⊗tp, sort(union(cq, tq)))
                splice!(circuit, i, (op_3, op_2, op_1))
            end
            _ => nothing
        end
    end
end

"""TODO docstring"""
function group_nonclifford(circuit::Circuit)
    if find_variant_indices(circuit,ExpEighPiPauli) != []
        for index in find_variant_indices(circuit,ExpEighPiPauli)
            circuit=traversal(circuit, conjugate_noncliff, :left, 1, index-1)
        end
    end
end

"""Identifies and combines identical Pauli rotations:
    For example, two PPR (π/8) on the same Pauli operator P are merged into a single Clifford-level PPR (π/4).
    A rotation and its inverse, PPR (π/8) and PPR (−π/8), cancel each other out completely and are removed."""
function merge_ops(circuit::Circuit)
    traversal(circuit,merge_rotations, :left, 1, :end)
end

"""TODO docstring"""
function remove_clifford(circuit::Circuit)
    validate_circuit(circuit)
    for index in find_variant_indices(circuit,Measurement)
        circuit=traversal(circuit, conjugate_measurement, :left, 1, index-1)
    end
    return circuit
end

"""TODO docstring"""
function remove_nonclifford(circuit::Circuit)
    indices=find_variant_indices(circuit,ExpEighPiPauli)
    num_input_qubit=get_circuit_width(circuit)
    num_magic_state=0
    for i in reverse(indices)
        num_magic_state+=1
        op=circuit[i]
        gadget = gadgetize(op, num_input_qubit, num_magic_state)
        splice!(circuit, i, gadget)
    end
end

"""TODO docstring"""
function remove_post_measurement(circuit::Circuit)
    # remove all gates after the last measurement
    index=maximum(find_variant_indices(circuit,Measurement))
    resize!(circuit, index)
end


##

"""ADT representing different types of measurement result"""
@data MeasurementResultType begin
    """Denoting measurement results that classically determined by a coin flip"""
    ClassicalDetermRes
    """Denoting measurement results that are classically determined by stored eigenvalues of stabilizers"""
    ClassicalRandomRes
    """Denoting measurement results that require performing actual quantum measurement"""
    QuantumRes
end

using .MeasurementResultType: ClassicalDetermRes, ClassicalRandomRes, QuantumRes

"""Struct holding measurement result value and its type"""
struct MeasurementResult
    """Single bit measurement result in boolean"""
    result::Union{Bool,Nothing}
    """Measurement result type of this result (ClassicalDetermRes, ClassicalRandomRes, QuantumRes)"""
    result_type::MeasurementResultType.Type
end

@derive MeasurementResultType[Hash, Eq, Show]

"""TODO docstring"""
classical_deterministic_result(m::Union{Bool,Nothing}) = MeasurementResult(m, ClassicalDetermRes())
"""TODO docstring"""
classical_random_result(m::Union{Bool,Nothing}) = MeasurementResult(m, ClassicalRandomRes())
"""TODO docstring"""
quantum_result(m::Union{Bool,Nothing}) = MeasurementResult(m, QuantumRes())

"""Struct that contains information describing current quantum state"""
struct MemoryState
    """Vector that holds all MeasurementResult"""
    measurement_results::Vector{MeasurementResult}
    """Stabilizer object that describes current quantum state"""
    stabilizer_group::Stabilizer
    """GeneralizedStabilizer object holding current quantum state within quantum computer"""
    quantum_memory::Union{GeneralizedStabilizer, Nothing}
    """Vector that holds all classical bits storing corresponding measurement results"""
    classical_register::Vector{Union{Nothing,Bool}}
end

"""Abstract base type for simulator execution states."""
abstract type AbstractSimState end

"""Execution state for real backend simulation."""
struct ComputerState <: AbstractSimState
    """Contain current circuit object"""
    circuit::Circuit
    """Number of gadgets inserted to replace nonclifford circuit operations"""
    num_gadgets::Int
    """Denote the Pauli Product Measurement that is being processed"""
    instruction_pointer::Int
    """Contain current quantum state"""
    memory_state::MemoryState
end

"""Execution state for dummy simulation. No real backend is contacted."""
struct DummyState <: AbstractSimState
    """Contain current circuit object"""
    circuit::Circuit
    """Number of gadgets inserted to replace nonclifford circuit operations"""
    num_gadgets::Int
    """Denote the Pauli Product Measurement that is being processed"""
    instruction_pointer::Int
    """Contain current quantum state"""
    memory_state::MemoryState
    """Weight vector describes sampling probability between +1 and -1 measurement results"""
    outcome_probs::Vector{Int}
end

include("joint_measurement_check.jl")
##

"""
    get_CompState(circuit::Circuit, input_state::Union{Stabilizer, Nothing}=nothing, dummy::Bool=false) -> AbstractSimState

Get initial Execution State using input circuit and input state.

# Fields
- `circuit`: input circuit for compilation
- `input_state`: initial qubit state
- `dummy`: flag for running dummy simulation
- `outcome_probs`: 2-element distribution vector [p_p1, p_m1]. p_p1 is the probability measuring +1; p_m1 is the probability measuring -1
"""
function get_CompState(circuit::Circuit, input_state::Union{Stabilizer, Nothing}=nothing; dummy::Bool=false, outcome_probs::Vector{Int}=[1,1])
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
    MeasRes=Vector{MeasurementResult}(undef, num_bits)
    creg=Array{Union{Nothing, Bool}}(nothing, num_bits)
    ms=MemoryState(MeasRes, stabilzier_group, magicstate, creg)
    if dummy
        return DummyState(circuit, num_gadgets, 1, ms, outcome_probs)
    else
        return ComputerState(circuit,num_gadgets, 1, ms)
    end
end

"""TODO docstring"""
abstract type QuantumRuntime end

"""TODO docstring -- all measurements return `nothing` and classically-trackable states are set as if result was `false`."""
struct MockRuntime <: QuantumRuntime end

"""Perform the next joint measurement and update ComputerState accordingly"""
function do_quantum_step(state::S, runtime::Type{<:QuantumRuntime}=MockRuntime) where S <: AbstractSimState
    # run the quantum measurement, appropriately updating MemoryState
    circuit = state.circuit
    i=state.instruction_pointer
    ms=state.memory_state
    @debug "Now working with $i th measurement" _group=:api
    meas_list = find_variant_indices(circuit,Measurement)
    meas_i=circuit[meas_list[i]]
    bit_index=meas_i.bit
    checklist=ms.stabilizer_group
    res=get_measurement_result(state, meas_i)
    @debug "Measurement result is $res" _group=:api
    MR=res[1]
    j=res[2]
    ms.measurement_results[i]=MR
    ms.classical_register[bit_index]=MR.result
    @match MR.result_type begin
        ClassicalDetermRes() => nothing
        QuantumRes() => begin
            @debug "This measurement outputs Quantum Result" _group=:api
            quantum_state=res[3]
            paulistring=embed(size(ms.stabilizer_group)[2], meas_i.qubits, meas_i.pauli)
            a_stabilizer= Stabilizer([paulistring])
            stabilizer_group=vcat(ms.stabilizer_group,a_stabilizer)
            @reset ms.stabilizer_group = stabilizer_group
            @reset ms.quantum_memory = quantum_state
        end
        ClassicalRandomRes() => begin
            @debug "This measurement outputs Classical Random Result" _group=:api
            q_1=[1:get_circuit_width(circuit);]
            Q_1=ExpQuatPiPauli(checklist[j],q_1)
            p_2=(-1)^MR.result*meas_i.pauli
            Q_2=ExpQuatPiPauli(p_2,meas_i.qubits)
            pushfirst!(circuit,Q_1,Q_2,Q_1)
            preprocess_circuit(circuit)
        end
    end
    @reset state.instruction_pointer = i+1
    @reset state.memory_state = ms
end

"""Run compute/compile with provided circuit and input state(described by stabilizer group)"""
function run(input_circuit::Circuit, input_state::Union{Stabilizer, Nothing}=nothing; dummy::Bool=false, outcome_probs::Vector{Int}=[1,1])
    state = get_CompState(input_circuit, input_state; dummy=dummy, outcome_probs=outcome_probs)
    len=length(state.memory_state.classical_register)
    while true && !isempty(state.circuit)
        @debug "Working on $(state.instruction_pointer) th PPM" _group=:api
        resolve_conditionals(state)
        state=do_quantum_step(state)
        @debug "Performed $(state.instruction_pointer) th PPM" _group=:api
        @debug "Current classical register: $(state.memory_state.classical_register)" _group=:api
        if state.instruction_pointer>len
            break
        end
    end
    @debug "Compute/Compile Complete" _group=:api
    pbc_circuit=[]
    for i in 1:length(state.circuit)
        op=state.circuit[i]
        pauli=state.memory_state.measurement_results[i].pauli
        new_op=CircuitOp.Measurement(pauli,op.bit,op.qubits)
        push!(pbc_circuit,new_op)
    end
    return @reset state.circuit = pbc_circuit
    @debug "Result returned" _group=:api
end

end # module PBCCompiler
