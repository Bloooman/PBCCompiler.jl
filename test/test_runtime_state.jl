@testitem "CompilerState copy isolates runtime state" tags=[:runtime] begin
##
using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, Measurement, ExpEighPiPauli, SimRuntime,
    DummyRuntime, DummyStabilizerRuntime, build_compilerstate, do_quantum_step
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
