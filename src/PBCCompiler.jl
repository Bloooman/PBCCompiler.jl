module PBCCompiler

using Moshi.Data: @data, variant_name, isa_variant
using Moshi.Match: @match
using Moshi.Derive: @derive
using QuantumClifford: PauliOperator, @P_str, comm, embed, ⊗, random_pauli, tensor, Stabilizer
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
function _group_nonclifford(circuit::Circuit)
    if _find_variant_indices(circuit,ExpEighPiPauli) != []
        for index in _find_variant_indices(circuit,ExpEighPiPauli)
            circuit=traversal(circuit, conjugate, :left, 1, index-1)
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
function _remove_clifford(circuit::Circuit)
    _validate_circuit(circuit)
    for index in _find_variant_indices(circuit,Measurement)
        circuit=traversal(circuit, conjugate, :left, 1, index-1)
    end
    return circuit
end

"""TODO docstring"""
function _remove_nonclifford(circuit::Circuit)
    indices=_find_variant_indices(circuit,ExpEighPiPauli)
    num_input_qubit=_get_circuit_width(circuit)
    num_magic_state=0
    for i in reverse(indices)
        num_magic_state+=1
        _gadgetize(circuit, i, num_input_qubit, num_magic_state)
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

"""TODO docstring"""
classical_deterministic_result(m::Union{Bool,Nothing}) = MeasurementResult(m, ClassicalDetermRes)
"""TODO docstring"""
classical_random_result(m::Union{Bool,Nothing}) = MeasurementResult(m, ClassicalRandomRes)
"""TODO docstring"""
quantum_result(m::Union{Bool,Nothing}) = MeasurementResult(m, QuantumRes)

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

"""TODO docstring"""
struct ComputerState
    """Contain current circuit object"""
    circuit::Circuit
    """Denote the Pauli Product Measurement that is being processed"""
    instruction_pointer::Int
    """Contain current quantum state"""
    memory_state::MemoryState
end

##

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

"""TODO docstring"""
function do_quantum_step(compstate::ComputerState, measurement, runtime::Type{<:QuantumRuntime}=MockRuntime)
    # run the quantum measurement, appropriately updating MemoryState
end

"""TODO docstring"""
function run(circuit::Circuit)
    # run preprocessing
    # prepare ComputerState
    while true
        # run next_quantum_step and do_quantum_step until there is no next step
    end
    return # measurement results
end

end # module PBCCompiler
