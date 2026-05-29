using QuantumClifford: PauliOperator, @P_str, Stabilizer, GeneralizedStabilizer
using Moshi.Data: @data, variant_name
using Moshi.Derive: @derive
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
##

"""ADT representing different types of measurement result"""
@data MeasurementResult begin
    """Denoting measurement results that classically determined by a coin flip"""
    struct ClassicalDetermRes
        """Corresponding Pauli String of Measurement"""
        pauli::PauliOperator
        """Single bit measurement result in boolean"""
        result::Union{Bool,Nothing}
    end
    """Denoting measurement results that are classically determined by stored eigenvalues of stabilizers"""
    struct ClassicalRandomRes
        """Corresponding Pauli String of Measurement"""
        pauli::PauliOperator
        """Single bit measurement result in boolean"""
        result::Union{Bool,Nothing}
    end
    """Denoting measurement results that require performing actual quantum measurement"""
    struct QuantumRes
        """Corresponding Pauli String of Measurement"""
        pauli::PauliOperator
        """Single bit measurement result in boolean"""
        result::Union{Bool,Nothing}
    end
end

using .MeasurementResult: ClassicalDetermRes, ClassicalRandomRes, QuantumRes

@derive MeasurementResult[Hash, Eq, Show]

"""Struct that contains information describing current compiler state"""
struct CompilerState
    """Vector that holds all MeasurementResult"""
    measurement_results::Vector{MeasurementResult.Type}
    """Stabilizer object that describes current quantum state"""
    stabilizer_group::Stabilizer
    """Vector that holds all classical bits storing corresponding measurement results"""
    classical_register::Vector{Union{Nothing,Bool}}
    """Contain current circuit object"""
    circuit::Circuit
    """Number of gadgets inserted to replace nonclifford circuit operations"""
    num_gadgets::Int
    """Denote the Pauli Product Measurement that is being processed"""
    instruction_pointer::Int
end
##
abstract type AbstractRuntime end

struct SimRuntime <: AbstractRuntime
    """Contain current compiler state"""
    compiler_state::CompilerState
    """GeneralizedStabilizer object holding current quantum state within quantum computer"""
    quantum_memory::Union{GeneralizedStabilizer, Nothing}
end

struct DummyRuntime <: AbstractRuntime
    """Contain current compiler state"""
    compiler_state::CompilerState
    """Weight vector describes sampling probability between +1 and -1 measurement results"""
    outcome_probs::Vector{Int}
end
##
function _result_type_str(t)
    t == :ClassicalDetermRes && return "ClassicalDeterministic"
    t == :ClassicalRandomRes && return "ClassicalRandom"
    t == :QuantumRes         && return "Quantum"
    return string(t)
end

_bool_str(::Nothing) = "nothing"
_bool_str(b::Bool) = string(b)

function _magic_pauli_str(p::PauliOperator, magic_qubits::Vector{Int})
    phase_char = p.phase[] in (0x00, 0x01) ? '+' : '-'
    chars = map(magic_qubits) do i
        x, z = p[i]
        x && z ? 'Y' : x ? 'X' : z ? 'Z' : '_'
    end
    return string(phase_char, String(chars))
end

"""
    show(io, result::ComputerState)

Pretty-print the first four fields of `result.memory_state`:
`pauli_qubits`, `magic_qubits`, `measurement_results`, and `StabilizerGroup`.
"""
function Base.show(io::IO, result::S) where S <: AbstractRuntime
    num_qubits = get_circuit_width(result.compiler_state.circuit)
    magicqubits = collect(num_qubits-result.compiler_state.num_gadgets+1:num_qubits)
    cs = result.compiler_state
    println(io, "Compilation Result")
    if !isdefined(cs, :measurement_results)
        println(io, "  Measurements: undefined")
        println(io, "  Quantum Measurement Results: undefined")
    else
        n = length(cs.measurement_results)
        println(io, "  Measurements ($n):")
        for i in 1:n
            if !isassigned(cs.measurement_results, i)
                println(io, "    [$i] undefined")
                continue
            end
            m = cs.measurement_results[i]
            println(io, "    [$i] ", m.pauli,
                        "  →  ", _bool_str(m.result),
                        "  (", _result_type_str(variant_name(m)), ")")
        end
        quantum = filter(i -> isassigned(cs.measurement_results, i) &&
                              isa_variant(cs.measurement_results[i], QuantumRes),
                         1:n)
        println(io, "  Quantum Measurement Results ($(length(quantum))):")
        for (j, i) in enumerate(quantum)
            m = cs.measurement_results[i]
            println(io, "    [$j] ", _magic_pauli_str(m.pauli, magicqubits),
                        "  →  ", _bool_str(m.result))
        end
    end
    print(io, "  Stabilizer Group:\n")
    show(io, cs.stabilizer_group)
end
