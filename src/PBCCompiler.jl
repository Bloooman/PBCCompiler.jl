module PBCCompiler

using Moshi.Data: @data, variant_name
using Moshi.Match: @match
using QuantumClifford: PauliOperator, @P_str, embed, comm, ⊗

##

"""TODO docstring"""
const P = typeof(P"XYZ")

"""TODO docstring"""
@data CircuitOp begin
    """TODO docstring"""
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
    """TODO docstring"""
    struct ExpHalfPiPauli
        pauli::P
        qubits::Vector{Int}
    end
    """TODO docstring"""
    struct ExpQuatPiPauli
        pauli::P
        qubits::Vector{Int}
    end
    """TODO docstring"""
    struct ExpEighPiPauli
        pauli::P
        qubits::Vector{Int}
    end
    """TODO docstring"""
    struct PrepMagic
        qubit::Int
        qubits::Vector{Int}
    end
    """TODO docstring"""
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

"""TODO docstring"""
const Circuit = Vector{CircuitOp.Type}

using .CircuitOp: Measurement, Pauli, ExpHalfPiPauli, ExpQuatPiPauli, ExpEighPiPauli, PrepMagic, PauliConditional, BitConditional

include("traversal.jl")
include("affectedqubits.jl")
include("plotting.jl")
include("pair_transformation.jl")
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
    remove_pauliconditional(circuit)
    commute_nonclifford_to_front(circuit)
    group_nonclifford(circuit)
    commute_measurements_to_end(circuit)
    remove_nonclifford(circuit)
    remove_post_measurement(circuit)
end

"""TODO docstring"""
function remove_nonclifford(circuit::Circuit)
    # introduce magic states and related measurements and conditional gates
end

"""TODO docstring"""
function remove_pauliconditional(circuit::Circuit)
    # turn pauli conditionals into expquatpipauli
end

"""TODO docstring"""
function remove_clifford(circuit::Circuit)
    # turn
end

"""TODO docstring"""
function remove_post_measurement(circuit::Circuit)
    # remove all gates after the last measurement
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
    memory_state::test_MemoryState
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
