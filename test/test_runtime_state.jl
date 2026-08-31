@testitem "CompilerState copy isolates runtime state" tags=[:runtime] begin
##
using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, Measurement, ExpEighPiPauli, SimRuntime,
    DummyRuntime, DummyStabilizerRuntime, HybridRuntime, StabilizerRuntime,
    DummyHybridRuntime, DummyHybridStabilizerRuntime, collapseRuntime,
    build_compilerstate, do_quantum_step
using QuantumClifford: @P_str

# One pi/8 rotation -> one magic-state gadget, so the runtime carries a
# GeneralizedStabilizer, an `activated` bit, and telemetry to be isolated.
gadget_circuit() = Circuit(CircuitOp.Type[
    CircuitOp.ExpEighPiPauli(P"Z", [1]),
    CircuitOp.Measurement(P"Z", 1, [1]),
])

@testset "copy does not alias magic-state memory" begin
    state = build_compilerstate(gadget_circuit(), SimRuntime(), nothing)
    @test state.runtime.quantum_memory !== nothing
    @test !any(state.runtime.activated)

    original_chi = length(state.runtime.quantum_memory.destabweights)

    other = copy(state)
    @test other.runtime !== state.runtime
    @test other.runtime.quantum_memory !== state.runtime.quantum_memory
    @test other.runtime.activated !== state.runtime.activated
    @test other.runtime.invsparsity_history !== state.runtime.invsparsity_history

    # Drive a full measurement on the copy. It activates a magic qubit (applying
    # the deferred T), projects the magic register, and records telemetry.
    other = do_quantum_step(other)
    @test any(other.runtime.activated)

    # None of that may be visible through the original, which has not measured
    # anything yet. Before the deep copy this failed: the original saw the
    # activation bit set and would skip its own deferred T gate.
    @test !any(state.runtime.activated)
    @test isempty(state.runtime.invsparsity_history)
    @test length(state.runtime.quantum_memory.destabweights) == original_chi
end

@testset "copy does not alias activated for DummyRuntime" begin
    state = build_compilerstate(gadget_circuit(), DummyRuntime(), nothing)
    @test !any(state.runtime.activated)

    other = copy(state)
    @test other.runtime !== state.runtime
    @test other.runtime.activated !== state.runtime.activated

    other = do_quantum_step(other)
    @test any(other.runtime.activated)
    # Must not be visible through the original -- same aliasing bug class as
    # SimRuntime's magic register, now for the dummy's activation bits.
    @test !any(state.runtime.activated)
end

@testset "copy does not alias activated for DummyStabilizerRuntime" begin
    state = build_compilerstate(gadget_circuit(), DummyStabilizerRuntime(), nothing)
    @test !any(state.runtime.activated)

    other = copy(state)
    @test other.runtime !== state.runtime
    @test other.runtime.activated !== state.runtime.activated

    other = do_quantum_step(other)
    @test any(other.runtime.activated)
    @test !any(state.runtime.activated)
end

@testset "copy does not alias magic-state memory for HybridRuntime" begin
    state = build_compilerstate(gadget_circuit(), HybridRuntime(), nothing)
    @test state.runtime.quantum_memory !== nothing
    @test !any(state.runtime.activated)

    original_chi = length(state.runtime.quantum_memory.destabweights)

    other = copy(state)
    @test other.runtime !== state.runtime
    @test other.runtime.quantum_memory !== state.runtime.quantum_memory
    @test other.runtime.activated !== state.runtime.activated
    @test other.runtime.invsparsity_history !== state.runtime.invsparsity_history

    other = do_quantum_step(other)
    @test any(other.runtime.activated)

    # Before HybridRuntime was added to `_RuntimeWithMutableFields`, `copy`
    # fell through to the identity fallback and this state was shared.
    @test !any(state.runtime.activated)
    @test isempty(state.runtime.invsparsity_history)
    @test length(state.runtime.quantum_memory.destabweights) == original_chi
end

@testset "copy does not alias activated for DummyHybridRuntime" begin
    state = build_compilerstate(gadget_circuit(), DummyHybridRuntime(), nothing)
    @test !any(state.runtime.activated)

    other = copy(state)
    @test other.runtime !== state.runtime
    @test other.runtime.activated !== state.runtime.activated

    other = do_quantum_step(other)
    @test any(other.runtime.activated)
    @test !any(state.runtime.activated)
end

@testset "copy does not alias activated for DummyHybridStabilizerRuntime" begin
    state = build_compilerstate(gadget_circuit(), DummyHybridRuntime(0), nothing)
    other = PBCCompiler.transition(state)
    @test other.runtime isa DummyHybridStabilizerRuntime
    @test !any(other.runtime.activated)

    copied = copy(other)
    @test copied.runtime !== other.runtime
    @test copied.runtime.activated !== other.runtime.activated

    copied = do_quantum_step(copied)
    @test any(copied.runtime.activated)
    @test !any(other.runtime.activated)
end

@testset "copy does not alias magic-state memory for collapseRuntime" begin
    state = build_compilerstate(gadget_circuit(), collapseRuntime(), nothing)
    @test state.runtime.quantum_memory !== nothing
    @test !any(state.runtime.activated)
    @test !any(state.runtime.collapsed)

    original_chi = length(state.runtime.quantum_memory.destabweights)

    other = copy(state)
    @test other.runtime !== state.runtime
    @test other.runtime.quantum_memory !== state.runtime.quantum_memory
    @test other.runtime.activated !== state.runtime.activated
    @test other.runtime.collapsed !== state.runtime.collapsed
    @test other.runtime.invsparsity_history !== state.runtime.invsparsity_history

    other = do_quantum_step(other)
    @test any(other.runtime.activated)

    # Before collapseRuntime was added to `_RuntimeWithMutableFields`, `copy`
    # fell through to the identity fallback and this state was shared -- a
    # second shot off `state` would have continued mutating the first shot's
    # already-projected magic register instead of starting fresh.
    @test !any(state.runtime.activated)
    @test isempty(state.runtime.invsparsity_history)
    @test length(state.runtime.quantum_memory.destabweights) == original_chi
end

@testset "copy leaves the tableau and circuit independent" begin
    state = build_compilerstate(gadget_circuit(), SimRuntime(), nothing)
    other = copy(state)
    @test other.stabilizer_group !== state.stabilizer_group
    @test other.circuit !== state.circuit
    @test other.classical_register !== state.classical_register
    @test other.measurement_results !== state.measurement_results
end

@testset "shots do not concatenate invsparsity telemetry" begin
    # `run` takes the runtime by value, but its contents used to be shared: the
    # history vector was never reset, so shot N returned shot 1..N's telemetry.
    rt = SimRuntime()
    circuit = gadget_circuit()
    lengths = Int[]
    for _ in 1:3
        result_state = PBCCompiler.run(copy(circuit), rt, nothing)
        push!(lengths, length(result_state.runtime.invsparsity_history))
    end
    @test allequal(lengths)
    @test isempty(rt.invsparsity_history)
end
##
end

@testitem "execution loop is type stable" tags=[:runtime] begin
##
using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, SimRuntime, build_compilerstate, execute!,
    do_quantum_step, get_measurement_result, quantum_measurement
using QuantumClifford: @P_str

state = build_compilerstate(Circuit(CircuitOp.Type[
    CircuitOp.ExpEighPiPauli(P"Z", [1]),
    CircuitOp.Measurement(P"Z", 1, [1]),
]), SimRuntime(), nothing)

# `run` itself is deliberately not checked: whether a circuit has magic-state
# gadgets is a runtime property, so the compiled state's type genuinely varies
# and `run` is the function barrier in front of it. What must stay concrete is
# everything *inside* the per-measurement loop -- when these infer as `Any`,
# every stabilizer-row commutation check becomes a dynamic dispatch.
@testset "per-measurement loop infers concretely" begin
    S = typeof(state)
    @test isconcretetype(S)
    @test Base.return_types(execute!, (S,)) == [S]
    @test Base.return_types(do_quantum_step, (S, Vector{Int})) == [S]

    rt = typeof(state.runtime)
    @test isconcretetype(rt)
    @test Base.return_types(quantum_measurement, (rt, CircuitOp.Type, Int)) == [Tuple{rt, Bool}]

    measres = only(Base.return_types(get_measurement_result, (S, CircuitOp.Type)))
    @test !(measres === Any)
    @test S <: measres.b.parameters[2]  # Union{Nothing, Tuple{MeasurementResult, S}}
end
##
end

@testitem "compile-once sampling matches repeated run" tags=[:runtime] begin
##
using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, SimRuntime, get_distribution,
    build_compilerstate, execute!, get_bit_number
using QuantumClifford: @P_str
using Random: seed!

# Two T gates and two measurements: enough for gadget measurements and for
# dependent outcomes, so a desynced magic register would show up in the counts.
circuit = Circuit(CircuitOp.Type[
    CircuitOp.ExpEighPiPauli(P"XY", [1, 2]),
    CircuitOp.ExpQuatPiPauli(P"ZX", [1, 2]),
    CircuitOp.ExpEighPiPauli(P"YZ", [1, 2]),
    CircuitOp.Measurement(P"Z", 1, [1]),
    CircuitOp.Measurement(P"Z", 2, [2]),
])

"Sampling loop as it was before compilation was amortized across shots."
function run_per_shot(input_circuit, rt, num_shots)
    num_bits = get_bit_number(input_circuit)
    data = zeros(Int, num_shots)
    for i in 1:num_shots
        register = PBCCompiler.run(copy(input_circuit), rt, nothing).classical_register
        bits = join(Int.(register[1:num_bits]))
        data[i] = parse(Int, bits; base=2)
    end
    return data
end

@testset "same outcomes as compiling per shot" begin
    for shots in (25, 100)
        seed!(4242); expected = run_per_shot(circuit, SimRuntime(), shots)
        seed!(4242); (_, actual) = get_distribution(circuit, SimRuntime(), nothing, shots)
        @test actual == expected
    end
end

@testset "sampling is reproducible under a seed" begin
    seed!(7); first_run  = get_distribution(circuit, SimRuntime(), nothing, 50)
    seed!(7); second_run = get_distribution(circuit, SimRuntime(), nothing, 50)
    @test first_run == second_run
end

@testset "sampling does not mutate the input circuit" begin
    before = copy(circuit)
    get_distribution(circuit, SimRuntime(), nothing, 10)
    @test length(circuit) == length(before)
    @test all(circuit[i] == before[i] for i in eachindex(circuit))
end

@testset "shots off one compilation are independent" begin
    # Each shot must start from the compiled state, not from the previous
    # shot's collapsed one. If `copy` aliased the magic register every shot
    # after the first would be forced onto one branch.
    compiled = build_compilerstate(copy(circuit), SimRuntime(), nothing)
    seed!(31)
    function run_shot(compiled)
        s = copy(compiled)
        while !PBCCompiler._execution_complete(s)
            s = execute!(s)
        end
        return s.classical_register[1:2]
    end
    outcomes = [run_shot(compiled) for _ in 1:60]
    @test !allequal(outcomes)
    # The compiled state itself is untouched: nothing measured yet.
    @test all(isnothing, compiled.classical_register)
    @test !any(compiled.runtime.activated)
end
##
end

@testitem "HybridRuntime converts into HybridStabilizerRuntime at its threshold" tags=[:runtime] begin
##
using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, HybridRuntime, HybridStabilizerRuntime,
    SimRuntime, StabilizerRuntime, to_result, num_gadget_qubits
using PBCCompiler: MeasurementResult, _magic_block_qpu_load
using .MeasurementResult: QuantumRes
using QuantumClifford: @P_str, nqubits
using Random: seed!

# Three independent pi/8 rotations -> three magic-state gadgets, so a
# threshold of 2 is crossed partway through the run (one gadget resolves
# pre-transition, two post-transition) rather than on the very first or last step.
three_gadget_circuit = Circuit(CircuitOp.Type[
    CircuitOp.ExpEighPiPauli(P"Z", [1]),
    CircuitOp.ExpEighPiPauli(P"Z", [2]),
    CircuitOp.ExpEighPiPauli(P"Z", [3]),
    CircuitOp.Measurement(P"Z", 1, [1]),
    CircuitOp.Measurement(P"Z", 2, [2]),
    CircuitOp.Measurement(P"Z", 3, [3]),
])

@testset "converts once total qubit support reaches the threshold" begin
    result = PBCCompiler.run(copy(three_gadget_circuit), HybridRuntime(2), nothing)
    @test result.runtime isa HybridStabilizerRuntime
    @test count(result.runtime.activated) >= 1

    out = to_result(result)
    num_input_qubits = nqubits(result.stabilizer_group) - num_gadget_qubits(result.runtime)
    expected_width = num_input_qubits + count(result.runtime.activated_at_transition)
    @test nqubits(out.stabilizer_group) == expected_width
    # activated_at_transition is a snapshot: it must not grow as `activated`
    # keeps being mutated for the rest of the run after conversion.
    @test count(result.runtime.activated_at_transition) <= count(result.runtime.activated)
    @test out.QPUDuration == length(out.QPU_workload)

    # Regression check for the pre/post-transition QPU_workload width
    # mismatch. `three_gadget_circuit`'s own pre-transition measurement only
    # ever touches one magic qubit, so `_magic_block_qpu_load` folds it away
    # entirely and the natural run above never actually exercises the padding
    # path. Force a pre-transition measurement with joint support on two
    # magic qubits (so it survives folding) to confirm it still comes out at
    # the same width as the post-transition entries, with the data-qubit
    # segment (qubit 1) forced to identity rather than leaking its Z support.
    num_qubits = nqubits(result.stabilizer_group)
    @test num_qubits == 6   # 3 data qubits + 3 gadgets, per three_gadget_circuit above
    padded_state = copy(result)
    padded_state.measurement_results[1] = QuantumRes(P"Z__XZ_", true)
    padded_out = to_result(padded_state)
    @test all(nqubits(entry.pauli) == num_qubits for entry in padded_out.QPU_workload)
    kept = padded_out.QPU_workload[1]
    @test all(kept.pauli[i] == (false, false) for i in 1:num_input_qubits)
end

@testset "never converts when maximum_measurement_support is nothing" begin
    seed!(1)
    result = PBCCompiler.run(copy(three_gadget_circuit), HybridRuntime(), nothing)
    @test result.runtime isa HybridRuntime

    # Regression check: a HybridRuntime that never converts must keep using
    # the generic to_result method, unaffected by the new dispatch. Seed both
    # runs identically so their random outcomes (and therefore QPU_workload)
    # line up.
    hybrid_out = to_result(result)
    seed!(1)
    sim_result = PBCCompiler.run(copy(three_gadget_circuit), SimRuntime(), nothing)
    sim_out = to_result(sim_result)
    @test nqubits(hybrid_out.stabilizer_group) == nqubits(sim_out.stabilizer_group)
    @test length(hybrid_out.QPU_workload) == length(sim_out.QPU_workload)
end

@testset "_magic_block_qpu_load embed_width pads with identity on the complement" begin
    magicqubits = 3:5
    embed_width = 6
    quantum = MeasurementResult.Type[
        QuantumRes(P"IIXZII", true),
        QuantumRes(P"IIIZYI", false),
    ]
    narrow = _magic_block_qpu_load(quantum, magicqubits)
    padded = _magic_block_qpu_load(quantum, magicqubits, embed_width)

    @test length(padded) == length(narrow)
    @test all(nqubits(entry.pauli) == embed_width for entry in padded)
    # Restricted back to magicqubits, padded entries agree with the narrow ones.
    @test all(padded[i].pauli[magicqubits] == narrow[i].pauli for i in eachindex(narrow))
    # The complement of magicqubits (the data-qubit segment) is identity.
    complement = setdiff(1:embed_width, magicqubits)
    @test all(entry.pauli[i] == (false, false) for entry in padded, i in complement)
end
##
end

@testitem "DummyHybridRuntime converts into DummyHybridStabilizerRuntime at its threshold" tags=[:runtime] begin
##
using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, DummyHybridRuntime, DummyHybridStabilizerRuntime,
    DummyRuntime, to_result, num_gadget_qubits
using QuantumClifford: @P_str, nqubits
using Random: seed!

# Same shape as the HybridRuntime testitem: three independent pi/8 rotations
# -> three magic-state gadgets, so a threshold of 2 is crossed partway through
# the run rather than on the very first or last step.
three_gadget_circuit = Circuit(CircuitOp.Type[
    CircuitOp.ExpEighPiPauli(P"Z", [1]),
    CircuitOp.ExpEighPiPauli(P"Z", [2]),
    CircuitOp.ExpEighPiPauli(P"Z", [3]),
    CircuitOp.Measurement(P"Z", 1, [1]),
    CircuitOp.Measurement(P"Z", 2, [2]),
    CircuitOp.Measurement(P"Z", 3, [3]),
])

@testset "converts once total qubit support reaches the threshold" begin
    result = PBCCompiler.run(copy(three_gadget_circuit), DummyHybridRuntime(2), nothing)
    @test result.runtime isa DummyHybridStabilizerRuntime
    @test count(result.runtime.activated) >= 1

    out = to_result(result)
    num_input_qubits = nqubits(result.stabilizer_group) - num_gadget_qubits(result.runtime)
    expected_width = num_input_qubits + count(result.runtime.activated_at_transition)
    @test nqubits(out.stabilizer_group) == expected_width
    # activated_at_transition is a snapshot: it must not grow as `activated`
    # keeps being mutated for the rest of the run after conversion.
    @test count(result.runtime.activated_at_transition) <= count(result.runtime.activated)
    @test out.QPUDuration == length(out.QPU_workload)
end

@testset "never converts when maximum_measurement_support is nothing" begin
    seed!(1)
    result = PBCCompiler.run(copy(three_gadget_circuit), DummyHybridRuntime(), nothing)
    @test result.runtime isa DummyHybridRuntime

    # Regression check: a DummyHybridRuntime that never converts must keep
    # using the generic to_result method, unaffected by the new dispatch.
    # Seed both runs identically so their random outcomes (and therefore
    # QPU_workload) line up.
    hybrid_out = to_result(result)
    seed!(1)
    dummy_result = PBCCompiler.run(copy(three_gadget_circuit), DummyRuntime(), nothing)
    dummy_out = to_result(dummy_result)
    @test nqubits(hybrid_out.stabilizer_group) == nqubits(dummy_out.stabilizer_group)
    @test length(hybrid_out.QPU_workload) == length(dummy_out.QPU_workload)
end
##
end

@testitem "collapseRuntime sizes QPU_workload from num_gadget_qubits, not the QuantumRes count" tags=[:runtime] begin
##
using PBCCompiler
using PBCCompiler: CircuitOp, collapseRuntime, CompilerState, MeasurementResult, to_result,
    _magic_block_qpu_load
using QuantumClifford: @P_str, MixedDestabilizer, Stabilizer, one

# The generic `to_result(state::CompilerState)` infers the magic-qubit window
# from `length(quantum)` (the QuantumRes count). That's correct for
# SimRuntime, where every gadget measurement is a QuantumRes, but
# collapseRuntime reclassifies some gadget measurements as ClassicalBiasedRes
# once their support collapses -- so the QuantumRes count can undercount the
# true number of gadget/magic qubits. Build a state by hand (2 magic qubits:
# qubits 3 and 4 of a 4-qubit register) with one QuantumRes spanning both
# magic qubits and one ClassicalBiasedRes, and check the magic-qubit window
# used to restrict/embed the QuantumRes Pauli is sized from `activated`
# (2 gadget qubits), not from the QuantumRes count (which is only 1).
@testset "to_result keeps both magic qubits' support on a joint QuantumRes" begin
    num_qubits = 4
    stab = MixedDestabilizer(one(Stabilizer, num_qubits; basis=:Z))
    rt = collapseRuntime(nothing, falses(2), falses(2), Int[])
    meas = MeasurementResult.Type[
        MeasurementResult.QuantumRes(P"__XX", false),
        MeasurementResult.ClassicalBiasedRes(P"__X_", true),
    ]
    state = CompilerState(; measurement_results=meas, stabilizer_group=stab,
        classical_register=Union{Nothing,Bool}[], circuit=CircuitOp.Type[],
        instruction_pointer=1, runtime=rt)

    result = to_result(state)
    # Before collapseRuntime had its own to_result, the generic method sized
    # the magic-qubit window from `length(quantum) == 1` instead of the true
    # `num_gadget_qubits(rt) == 2`, silently dropping the X on qubit 3 and
    # returning an empty QPU_workload.
    @test length(result.QPU_workload) == 1
    # Restricted to the 2-qubit magic-block width (qubits 3:4), not embedded
    # back to the full 4-qubit register.
    @test result.QPU_workload[1].pauli == P"XX"
    @test result.QPUDuration == 1
end
##
end

@testitem "collapseRuntime end-to-end: isolated gadget touches classify as ClassicalBiasedRes, not QuantumRes" tags=[:runtime] begin
##
using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, collapseRuntime, MeasurementResult, build_compilerstate,
    _execution_complete, execute!, to_result
using QuantumClifford: @P_str
using Random: seed!
using Moshi.Data: isa_variant

# A single pi/8 rotation measured by itself never shares its magic qubit with
# another gadget, so `_mark_collapsed!` finds exactly one remaining
# non-identity (not-yet-activated) magic qubit on the very first touch and
# marks it collapsed -- the measurement is classified `ClassicalBiasedRes`
# instead of `QuantumRes`. This is intended: an isolated T-gadget's statistics
# don't depend on entangling it with the live magic register, so it never
# needs to be actual QPU work.
circuit = Circuit(CircuitOp.Type[
    CircuitOp.ExpEighPiPauli(P"Z", [1]),
    CircuitOp.Measurement(P"Z", 1, [1]),
])

@testset "isolated gadget produces ClassicalBiasedRes and empty QPU_workload" begin
    seed!(1)
    state = build_compilerstate(circuit, collapseRuntime(), nothing)
    while !_execution_complete(state)
        state = execute!(state)
    end
    assigned = [state.measurement_results[i] for i in 1:length(state.measurement_results)
                if isassigned(state.measurement_results, i)]
    @test any(mr -> isa_variant(mr, MeasurementResult.ClassicalBiasedRes), assigned)
    @test !any(mr -> isa_variant(mr, MeasurementResult.QuantumRes), assigned)
    @test any(state.runtime.collapsed)
    @test any(state.runtime.activated)

    result = to_result(state)
    @test isempty(result.QPU_workload)
    @test result.QPUDuration == 0
end
##
end
