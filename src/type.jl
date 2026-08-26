using QuantumClifford: PauliOperator, @P_str, Stabilizer, GeneralizedStabilizer, MixedDestabilizer, stabilizerview
using Moshi.Data: @data, variant_name
using Moshi.Derive: @derive
##
"""Type of a QuantumClifford Pauli operator as produced by the `P"..."` string macro."""
const P = typeof(P"XYZ")

"""
Algebraic data type of the circuit operations understood by the compiler.

Variants: `Measurement`, `Pauli`, `ExpHalfPiPauli`, `ExpQuatPiPauli`, `ExpEighPiPauli`,
`PrepMagic`, `PauliConditional`, and `BitConditional`. Use
`Moshi.Data.isa_variant(op, CircuitOp.Pauli)` (not `op isa CircuitOp.Pauli`) to check
which variant an operation is.
"""
@data CircuitOp begin
    """Measurement of pauli string P (ie., + XY) on qubits in vector at field "qubits" (ie.,[1,3]), measurement result is stored in classical bit denoted in "bit" """
    struct Measurement
        pauli::P
        bit::Int
        qubits::Vector{Int}
    end
    """Apply the Pauli gate given by string `pauli` to the qubits in `qubits`."""
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
    """Prepare a magic (|T⟩) ancilla state on `qubit`, entangling with the register in `qubits`."""
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
    """Apply the wrapped operation `op` only if classical register bit `bit` reads 1."""
    struct BitConditional
        op::CircuitOp
        bit::Int
    end
end

@derive CircuitOp[Hash, Eq, Show]

"""Sequence of abstract circuit operations."""
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
        result::Bool
    end
    """Denoting measurement results that are classically determined by stored eigenvalues of stabilizers"""
    struct ClassicalRandomRes
        """Corresponding Pauli String of Measurement"""
        pauli::PauliOperator
        """Single bit measurement result in boolean"""
        result::Bool
    end
    """Denoting measurement results that require performing actual quantum measurement"""
    struct QuantumRes
        """Corresponding Pauli String of Measurement"""
        pauli::PauliOperator
        """Single bit measurement result in boolean"""
        result::Bool
    end
end

using .MeasurementResult: ClassicalDetermRes, ClassicalRandomRes, QuantumRes

@derive MeasurementResult[Hash, Eq, Show]
##
"""Supertype of the measurement backends a circuit can be compiled/computed against."""
abstract type AbstractRuntime end

"""Runtime that simulates magic-state measurements with QuantumClifford's `GeneralizedStabilizer`."""
struct SimRuntime{Q} <: AbstractRuntime
    # Q is the concrete GeneralizedStabilizer instantiation (or Nothing before a
    # magic register exists). Declaring the field as the unparameterized
    # `GeneralizedStabilizer` instead would make it non-concrete, which is
    # enough on its own to infer `run` and `do_quantum_step` as `Any`.
    """GeneralizedStabilizer object holding current quantum state within quantum computer"""
    quantum_memory::Q
    """
    Magic qubits whose deferred T gate has already been applied. The magic
    register starts as the stabilizer state |+>^n and each T gate is applied
    lazily, right before the first measurement touching its qubit, keeping the
    chi-expansion of `quantum_memory` at 4^(live qubits) instead of 4^(total)
    """
    activated::Union{BitVector, Nothing}
    """Number of non-zero elements in the density matrix at each simulation"""
    invsparsity_history::Vector{Int}
end

SimRuntime() = SimRuntime(nothing, nothing, Int[])

"""Supertype of runtimes that classify every non-deterministic outcome by
projecting the full register (no anticommuting coin-flip branch)."""
abstract type AbstractStabilizerRuntime <: AbstractRuntime end

"""
Runtime that simulates the full register — data qubits and magic qubits together —
with QuantumClifford's `GeneralizedStabilizer`.

Unlike [`SimRuntime`](@ref) — whose `quantum_memory` is the same full-width
register, but whose data-qubit part sits inert and unmeasured (the sign of the
data part is instead read off `state.stabilizer_group`) — this runtime
performs an actual quantum measurement for every non-deterministic outcome,
and resolves an anticommuting stabilizer row by projecting rather than by
splicing compensating rotations into the circuit. It therefore never produces
a `ClassicalRandomRes` — every non-deterministic outcome is recorded as a
`QuantumRes`.

That is the intended design, not a defect. The consequence is that
`QuantumRes` here means "not determined by the stabilizer group", which is a
broader class than `SimRuntime`'s "needs the magic register": a measurement that
is random purely because of the data qubits is a `QuantumRes` under this runtime
and a `ClassicalRandomRes` under `SimRuntime`. Expect a substantially larger
`QPU_workload` as a result (roughly 5x on random small circuits). Measurement
*outcomes* are unaffected and agree with `SimRuntime`; only the classification
and the metrics derived from it differ, so do not compare QPU-workload figures
across the two runtimes.
"""
struct StabilizerRuntime{Q} <: AbstractStabilizerRuntime
    # See `SimRuntime` on why the memory type is a parameter
    """GeneralizedStabilizer object holding current quantum state within quantum computer"""
    quantum_memory::Q
    """
    Magic qubits whose deferred T gate has already been applied. The magic
    register starts as the stabilizer state |+>^n and each T gate is applied
    lazily, right before the first measurement touching its qubit, keeping the
    chi-expansion of `quantum_memory` at 4^(live qubits) instead of 4^(total)
    """
    activated::Union{BitVector, Nothing}
    """Number of non-zero elements in the density matrix at each simulation"""
    invsparsity_history::Vector{Int}
end

StabilizerRuntime() = StabilizerRuntime(nothing, nothing, Int[])

"""Runtime that replaces quantum measurements with classical coin flips of a
fixed bias, while still tracking which magic qubits a measurement touched
(mirrors `SimRuntime.activated`) for parity/diagnostics."""
struct DummyRuntime <: AbstractRuntime
    """Probability of sampling the +1 measurement outcome (the -1 outcome has probability `1 - p1_outcome_probs`)"""
    p1_outcome_probs::Float64
    """Magic qubits touched by a measurement so far, mirroring
    `SimRuntime.activated`; `nothing` before `build_rt_data` runs (or when the
    circuit has no gadgets)."""
    activated::Union{BitVector, Nothing}
end

DummyRuntime() = DummyRuntime(0.5, nothing)
DummyRuntime(p::Float64) = DummyRuntime(p, nothing)

"""Cheap stand-in for `StabilizerRuntime`: same control flow (no anticommuting
coin-flip branch), same `activated` bookkeeping of which magic qubits a
measurement touched, but replaces the actual quantum measurement with a
classical coin flip of a fixed bias instead of simulating the register.
`to_result`/`QPU_workload` extraction work the same way they do for
`StabilizerRuntime`; only the outcome bias is not physically faithful."""
struct DummyStabilizerRuntime <: AbstractStabilizerRuntime
    """Probability of sampling the +1 measurement outcome (the -1 outcome has probability `1 - p1_outcome_probs`)"""
    p1_outcome_probs::Float64
    """Magic qubits touched by a measurement so far, mirroring
    `StabilizerRuntime.activated`; `nothing` before `build_rt_data` runs."""
    activated::Union{BitVector, Nothing}
end

DummyStabilizerRuntime() = DummyStabilizerRuntime(0.5, nothing)
DummyStabilizerRuntime(p::Float64) = DummyStabilizerRuntime(p, nothing)

"""
Runtime that starts out simulating the magic register like `SimRuntime`, then
converts in place into a `StabilizerRuntime` once `maximum_measurement_support`
activated magic qubits have been reached.

`PBCCompiler.run` drives this conversion by calling `transition` after every
measurement step. With `maximum_measurement_support = nothing` (the default),
it never converts and behaves exactly like `SimRuntime` for the whole run.
"""
struct HybridRuntime{Q} <: AbstractRuntime
    """GeneralizedStabilizer object holding current quantum state within quantum computer"""
    quantum_memory::Q
    """Magic qubits whose deferred T gate has already been applied, mirroring `SimRuntime.activated`"""
    activated::Union{BitVector, Nothing}
    """Number of non-zero elements in the density matrix at each simulation"""
    invsparsity_history::Vector{Int}
    """Number of activated magic qubits at which this runtime converts into a `StabilizerRuntime`; `nothing` means never convert"""
    maximum_measurement_support::Union{Int, Nothing}
end

HybridRuntime() = HybridRuntime(nothing, nothing, Int[], nothing)
HybridRuntime(m::Int64) = HybridRuntime(nothing, nothing, Int[], m)

"""Runtime that samples every measurement outcome classically with a fixed bias, used to traverse compilation branches deterministically."""
struct TraversalRuntime <: AbstractRuntime
    """Probability of sampling the +1 measurement outcome (the -1 outcome has probability `1 - p1_outcome_probs`)"""
    p1_outcome_probs::Float64
end

TraversalRuntime() = TraversalRuntime(1)
##
"""Struct that contains information describing current compiler state"""
Base.@kwdef struct CompilerState{R<:AbstractRuntime, T<:MixedDestabilizer}
    # R comes first so existing `CompilerState{<:SomeRuntime}` dispatch keeps
    # working. T is a parameter for the same reason as SimRuntime's Q: an
    # unparameterized `MixedDestabilizer` field is not concrete, which makes
    # `stabilizerview(state.stabilizer_group)` non-inferrable and turns the
    # per-row commutation scan in `get_measurement_result` into one dynamic
    # dispatch per stabilizer row per measurement.
    """Vector that holds all MeasurementResult in temporal order -- first measurement result is the first CircuitOp.Measurement being measured"""
    measurement_results::Vector{MeasurementResult.Type}
    """
    MixedDestabilizer tableau that describes the current quantum state
    It spans n qubits where n is the number of total qubits (magic and stabilizer state qubits)
    Its rank is below n before compilation is finished; the destabilizer half
    makes measurement projections cheap (no re-canonicalization per measurement)
    """
    stabilizer_group::T
    """Result of Measurement(..., bit, ...) is stored in classical_register[bit]"""
    classical_register::Vector{Union{Nothing,Bool}}
    """Contain current circuit object"""
    circuit::Circuit
    """Denote the Pauli Product Measurement that is being processed"""
    instruction_pointer::Int
    """Runtime backend used to obtain measurement results"""
    runtime::R
end

"""
    copy(rt::AbstractRuntime) -> AbstractRuntime

Copy a runtime deeply enough that measurements on the copy cannot be observed
through the original.

`TraversalRuntime` holds only a probability, so it is returned as-is.
`SimRuntime`, `StabilizerRuntime`, `DummyRuntime`, `DummyStabilizerRuntime`, and
`HybridRuntime` are immutable structs but their *contents* are not: `activated`
is set element-wise (and, for the runtimes with a `quantum_memory`,
`quantum_memory` is projected in place and `invsparsity_history` is appended
to). Sharing any of those between two states makes the second one skip its
deferred T gate (or, for the dummies, misreport which magic qubits were
touched) — a silently wrong outcome with no error.
"""
Base.copy(rt::AbstractRuntime) = rt

# `SimRuntime`/`StabilizerRuntime`/`DummyRuntime`/`DummyStabilizerRuntime`/
# `HybridRuntime` all copy the same way -- copy every field, nil-guarded since
# `activated`/`quantum_memory` can legitimately be `nothing` -- and differ only
# in their field lists, so one reflection-driven method covers all five
# instead of one hand-written copy per type.
const _RuntimeWithMutableFields = Union{SimRuntime,StabilizerRuntime,DummyRuntime,DummyStabilizerRuntime,HybridRuntime}

function Base.copy(rt::R) where R<:_RuntimeWithMutableFields
    fields = ntuple(fieldcount(R)) do i
        f = getfield(rt, i)
        f === nothing ? nothing : copy(f)
    end
    return R(fields...)
end

function Base.copy(s::CompilerState)
    CompilerState(
        copy(s.measurement_results),
        copy(s.stabilizer_group),
        copy(s.classical_register),
        copy(s.circuit),
        s.instruction_pointer,
        copy(s.runtime)
    )
end
##
"""Final outcome of compiling/running a circuit, produced by [`to_result`](@ref)."""
struct CompilationResult
    """Vector that holds all MeasurementResult in temporal order -- first measurement result is the first CircuitOp.Measurement being measured"""
    measurement_results::Vector{MeasurementResult.Type}
    """Joint Pauli measurements that must be executed on the QPU, restricted to the magic-state qubits"""
    QPU_workload::Vector{MeasurementResult.Type}
    """
    Stabilizer object that describes current quantum state
    It has n columns where n is the number of total qubits (magic and stabilizer state qubits)
    The tableau is not square (full-rank) before compilation is finished
    """
    stabilizer_group::Stabilizer
    """Number of sequential joint measurements the QPU must perform"""
    QPUDuration::Int
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

"""
Print one indented line per measurement result in `measurement_results`
(`"undefined"` for unassigned slots), shared by the `CompilerState` and
`CompilationResult` `show` methods.
"""
function _print_measurements(io::IO, measurement_results)
    n = length(measurement_results)
    println(io, "  Measurements ($n):")
    for i in 1:n
        if !isassigned(measurement_results, i)
            println(io, "    [$i] undefined")
            continue
        end
        m = measurement_results[i]
        println(io, "    [$i] ", m.pauli,
                    "  →  ", _bool_str(m.result),
                    "  (", _result_type_str(variant_name(m)), ")")
    end
end

"""
    show(io::IO, result::CompilerState)

Pretty-print all debug-relevant fields of a `CompilerState`.

# Arguments
- `io`: output stream
- `result`: compiler state to display

# Returns
Nothing; writes to `io`.
"""
function Base.show(io::IO, result::CompilerState)
    println(io, "CompilerState [debug]")
    println(io, "  Instruction Pointer: ", result.instruction_pointer)
    _print_measurements(io, result.measurement_results)
    reg = join(map(_bool_str, result.classical_register), ", ")
    println(io, "  Classical Register: [", reg, "]")
    print(io, "  Stabilizer Group:\n")
    show(io, stabilizerview(result.stabilizer_group))
end

"""
    show(io::IO, r::CompilationResult)

Pretty-print a `CompilationResult`.

# Arguments
- `io`: output stream
- `r`: compilation result to display

# Returns
Nothing; writes to `io`.
"""
function Base.show(io::IO, r::CompilationResult)
    println(io, "Compilation Result")
    _print_measurements(io, r.measurement_results)
    nq = length(r.QPU_workload)
    println(io, "  QPU Workload ($nq):")
    for (j, m) in enumerate(r.QPU_workload)
        println(io, "    [$j] ", m.pauli,
                    "  →  ", _bool_str(m.result))
    end
    println(io, "  QPU Duration: ", r.QPUDuration)
    print(io, "  Stabilizer Group:\n")
    show(io, r.stabilizer_group)
end
