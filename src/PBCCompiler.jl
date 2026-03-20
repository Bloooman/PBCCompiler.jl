module PBCCompiler

using Moshi.Data: @data, variant_name, isa_variant
using Moshi.Match: @match
using QuantumClifford: PauliOperator, @P_str, comm, embed, ⊗, random_pauli, tensor, Stabilizer
using Random: randstring
using StatsBase: sample

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

"""TODO docstring"""
function preprocess_circuit(circuit::Circuit)
    remove_pauliconditional(circuit)
    group_nonclifford(circuit)
    merge_ops(circuit)
    remove_clifford(circuit)
    remove_nonclifford(circuit)
    remove_post_measurement(circuit)
end

"""TODO docstring"""
function remove_nonclifford(circuit::Circuit)
    num_non_clifford=length(find_nonclifford_indices(circuit))
    for i in 1:num_non_clifford
        gadgetize(circuit)
    end
end

"""TODO docstring"""
function merge_ops(circuit::Circuit)
    traversal(circuit,merge_rotations, :left, 1, :end)
end

"""TODO docstring"""
function remove_pauliconditional(circuit::Circuit)
    # turn pauli conditionals into expquatpipauli
end

"""TODO docstring"""
function remove_clifford(circuit::Circuit)
    validate_circuit(circuit)
    for index in find_measurement_indices(circuit)
        circuit=traversal(circuit, conjugate, :left, 1, index-1)
    end
    return circuit
end

"""TODO docstring"""
function remove_post_measurement(circuit::Circuit)
    # remove all gates after the last measurement
    index=maximum(find_measurement_indices(circuit))
    resize!(circuit, index)
end

##

"""TODO docstring"""
@data MeasurementResultType begin
    """TODO docstring"""
    ClassicalDetermRes
    """TODO docstring"""
    ClassicalRandomRes
    """TODO docstring"""
    QuantumRes
end

using .MeasurementResultType: ClassicalDetermRes, ClassicalRandomRes, QuantumRes

"""TODO docstring"""
struct MeasurementResult
    """TODO docstring"""
    result::Union{Bool,Nothing}
    """TODO docstring"""
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
    """TODO docstring"""
    circuit::Circuit
    """TODO docstring"""
    instruction_pointer::Int
    """TODO docstring"""
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
