module PBCCompiler

using Moshi.Data: @data, variant_name, isa_variant
using Moshi.Match: @match
using Moshi.Derive: @derive
using QuantumClifford: PauliOperator, @P_str, comm, embed, ⊗, random_pauli, tensor, @S_str, Stabilizer, project!
using Random: randstring
using StatsBase: sample

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

"""TODO docstring"""
struct MemoryState
    """TODO docstring"""
    measurement_results::Dict{Int,MeasurementResult}
    """TODO docstring"""
    pauli_qubits::Vector{Int}
    """TODO docstring"""
    pauli_state::P
    """TODO docstring"""
    magic_qubits::Vector{Int}
    """TODO docstring"""
    magic_state::Any
end

"""Struct that contains information describing current quantum state"""
struct test_MemoryState
    """Vector that contains index of all data qubits that hold circuit input"""
    pauli_qubits::Vector{Int}
    """Vector that contains index of all qubits that hold magic states"""
    magic_qubits::Vector{Int}
    """Vector that holds all MeasurementResult"""
    measurement_results::Vector{MeasurementResult}
    """Stabilizer object that describes current quantum state"""
    StabilizerGroup::Stabilizer
    """Vector that holds all classical bits storing corresponding measurement results"""
    classical_register::Vector{Union{Nothing,Bool}}
end


"""Struct that contains current state of compiler"""
struct ComputerState
    """Contain current circuit object"""
    circuit::Circuit
    """Denote the Pauli Product Measurement that is being processed"""
    instruction_pointer::Int
    """Contain current quantum state"""
    memory_state::test_MemoryState
end

include("joint_measurement_check.jl")
##

"""Get initial ComputerState using input circuit and input state"""
function get_CompState(circuit::Circuit, input_state::Stabilizer)
    num_pauli_qubits=get_circuit_width(circuit)
    pauliqubits=Int[1:num_pauli_qubits;]
    preprocess_circuit(circuit)
    magicqubits=Int[num_pauli_qubits+1: get_circuit_width(circuit);]
    num_bits=get_bit_number(circuit)
    MeasRes=Vector{MeasurementResult}(undef, num_bits)
    creg=Array{Union{Nothing, Bool}}(nothing, num_bits)
    stabilzier_group=make_stabilizer_list(input_state, circuit)
    MS=test_MemoryState(pauliqubits, magicqubits, MeasRes, stabilzier_group, creg)
    CS=ComputerState(circuit, 1, MS)
    return CS
end

"""TODO docstring"""
function next_quantum_step(compstate::ComputerState)
    while true
        # resolve conditionals
        # find next measurement -- if there is none, return nothing
        # commute measurement through preceding gates
        # check, given knowledge of the memory, whether the measurement is known or 50/50 random
        #   - if yes, store measurement result and update classically-trackable computer state
        #   - if not, break and return the measurement to perform on the quantum computer
    end
end

"""TODO docstring"""
abstract type QuantumRuntime end

"""TODO docstring -- all measurements return `nothing` and classically-trackable states are set as if result was `false`."""
struct MockRuntime <: QuantumRuntime end

"""Perform the next joint measurement and update ComputerState accordingly"""
function do_quantum_step(compstate::ComputerState, runtime::Type{<:QuantumRuntime}=MockRuntime)
    # run the quantum measurement, appropriately updating MemoryState
    circuit=compstate.circuit
    i=compstate.instruction_pointer
    MS=compstate.memory_state
    @debug("Now working with $i th measurement")
    meas_list = find_variant_indices(circuit,Measurement)
    meas_i=circuit[meas_list[i]]
    bit_index=meas_i.bit
    checklist=MS.StabilizerGroup
    (MR,j)=get_measurement_result(checklist, meas_i, get_circuit_width(circuit))
    MS.measurement_results[bit_index]=MR
    MS.classical_register[bit_index]=MR.result
    @match MR.result_type begin
        ClassicalDetermRes() => nothing
        QuantumRes() => begin
            @debug("This measurement outputs Quantum Result")
            paulistring=embed(size(MS.StabilizerGroup)[2], meas_i.qubits, meas_i.pauli)
            a_stabilizer= Stabilizer([paulistring])
            StabilizerGroup=vcat(MS.StabilizerGroup,a_stabilizer)
            MS=test_MemoryState(MS.pauli_qubits, MS.magic_qubits, MS.measurement_results, StabilizerGroup, MS.classical_register)
        end
        ClassicalRandomRes() => begin
            @debug("This measurement outputs Classical Random Result")
            q_1=[1:get_circuit_width(circuit);]
            Q_1=ExpQuatPiPauli(checklist[j],q_1)
            p_2=(-1)^MR.result*meas_i.pauli
            Q_2=ExpQuatPiPauli(p_2,meas_i.qubits)
            pushfirst!(circuit,Q_1,Q_2,Q_1)
            preprocess_circuit(circuit)
        end
    end
    i=i+1
    return ComputerState(circuit, i, MS)
end

"""Run compute/compile with provided circuit and input state(described by stabilizer group)"""
function run(input_circuit::Circuit, input_state::Stabilizer)
    # run preprocessing
    # prepare ComputerState
    validate_circuit(input_circuit)
    validate_input(input_circuit,input_state)
    CS = get_CompState(input_circuit, input_state)
    len=length(CS.memory_state.classical_register)
    while true && !isempty(CS.circuit)
        @debug("Working on $(CS.instruction_pointer) th PPM")
        # run next_quantum_step and do_quantum_step until there is no next step
        resolve_conditionals(CS)
        @debug("After BitConditional resolved, the circuit becomes: \n$(join(CS.circuit, "\n"))")
        CS=do_quantum_step(CS)
        @debug("Performed $pointer th PPM")
        @debug("After PPM resolved, the circuit becomes: \n$(join(CS.circuit, "\n"))")
        if CS.instruction_pointer>len
            break
        end
        @debug("Current classical register: $(CS.memory_state.classical_register)")
    end
    return CS
end

end # module PBCCompiler
