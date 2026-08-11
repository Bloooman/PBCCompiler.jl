@testitem "to_result edge cases" tags=[:to_result] begin

using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, to_result, run, random_test_circuit,
    CompilerState, CompilationResult, MeasurementResult, SimRuntime, DummyRuntime,
    TraversalRuntime, get_circuit_width, make_stabilizer_list, validate_input
using .CircuitOp: Measurement, ExpEighPiPauli, ExpQuatPiPauli
using .MeasurementResult: ClassicalDetermRes, ClassicalRandomRes, QuantumRes
using Moshi.Match: isa_variant
using QuantumClifford: @P_str, @S_str, Stabilizer, MixedDestabilizer, nqubits, one
using Random

# Build a CompilerState directly so the magic register is known independently of
# how many measurements happened to come back quantum.
function _state(measurement_results, num_qubits)
    md = MixedDestabilizer(one(Stabilizer, num_qubits; basis=:Z))
    CompilerState(
        Vector{MeasurementResult.Type}(measurement_results),
        md,
        Union{Nothing,Bool}[nothing for _ in measurement_results],
        Circuit(),
        length(measurement_results) + 1,
        DummyRuntime(),
    )
end

# to_result recovers the magic-qubit columns as
#     num_qubits - count(QuantumRes) + 1 : num_qubits
# while build_rt_data derives the same register structurally as
#     get_circuit_width(preprocessed) - num_input_qubits.
# The two agree only while every gadget measurement returns QuantumRes and no
# other measurement does. That holds for a full-rank input state (a data-only
# Pauli commuting with every row of a full-rank group is in the group, hence
# ClassicalDetermRes), which is why the scan below passes -- but it is an
# assumption, not an enforced invariant. This testset pins it: if a future change
# breaks it, to_result starts slicing the wrong columns silently.
@testset "magic-register invariant holds for compiled circuits" begin
    for (rt_factory, seeds) in ((SimRuntime, 1:20), (DummyRuntime, 101:160),
                                (TraversalRuntime, 201:260))
        for seed in seeds
            Random.seed!(seed)
            input = random_test_circuit(20, 3)
            width_in = get_circuit_width(input)
            state = run(copy(input), rt_factory())
            isempty(state.circuit) && continue
            num_qubits = nqubits(state.stabilizer_group)
            num_magic = num_qubits - width_in
            quantum = filter(mr -> isa_variant(mr, QuantumRes), state.measurement_results)
            # If this ever fails, to_result silently slices the wrong columns
            @test length(quantum) == num_magic
        end
    end
end

# Same measurement Paulis, same 4-qubit register (data 1:2, magic 3:4). The only
# difference is that the first measurement resolved classically instead of
# quantumly -- something the runtime decides per shot. That must not change how
# the magic support of the *other* measurement is read.
#
# Reachability: this is masked today by the make_stabilizer_list bug in the last
# testset. Fix that (use nqubits(s), not length(s)) and the rank-deficient
# circuit below compiles to num_magic=1 with TWO QuantumRes -- to_result then
# infers the magic range as 2:3 instead of 3:3 and analyses data qubit 2 as if it
# were a magic qubit. The two defects hide each other.
@testset "magic-register range must not depend on the QuantumRes count" begin
    data_only = P"ZZII"
    across_magic = P"IIZZ"   # weight 2 on the true magic register -> genuine QPU load

    both_quantum = to_result(_state(
        [QuantumRes(data_only, false), QuantumRes(across_magic, false)], 4))
    one_classical = to_result(_state(
        [ClassicalDetermRes(data_only, false), QuantumRes(across_magic, false)], 4))

    # With 2 QuantumRes the inferred range is 3:4 (correct) and the weight-2
    # measurement is real QPU load.
    @test any(mr -> mr.pauli == across_magic[3:4], both_quantum.QPU_workload)

    # With 1 QuantumRes the inferred range collapses to 4:4, the measurement
    # looks like a single-qubit local, and the QPU round disappears.
    @test_broken any(mr -> mr.pauli == across_magic[3:4], one_classical.QPU_workload)
    @test_broken one_classical.QPUDuration == 1
end

# A QuantumRes whose support on the magic register is entirely identity needs no
# magic qubit, so it is not a QPU round. The `length(non_id) == 1` test in
# to_result excludes single-qubit locals but lets the empty case fall through to
# the else branch, where it is pushed onto QPU_workload and counted.
@testset "identity support on the magic register is not QPU load" begin
    # 3 qubits, one magic qubit (data 1:2, magic 3) and one quantum measurement,
    # so the inferred range 3:3 IS the true magic register -- this isolates the
    # empty-support branch from the range-inference issue above.
    r = to_result(_state([QuantumRes(P"ZZI", false)], 3))
    # P"ZZI" is identity on qubit 3: no magic qubit is consumed, so there is
    # nothing for the QPU to measure jointly
    @test_broken isempty(r.QPU_workload)
    @test_broken r.QPUDuration == 0
end

# `show(::CompilationResult)` deliberately handles undefined measurement slots
# (see test_bugfixes.jl), so a partially filled measurement_results is a state
# the package expects to exist. to_result's `filter` over the raw vector cannot
# survive it.
@testset "to_result tolerates a partially filled measurement_results" begin
    mres = Vector{MeasurementResult.Type}(undef, 2)
    mres[1] = QuantumRes(P"ZZ", false)
    state = CompilerState(mres, MixedDestabilizer(one(Stabilizer, 2; basis=:Z)),
                          Union{Nothing,Bool}[nothing, nothing], Circuit(), 2, DummyRuntime())
    @test_broken to_result(state) isa CompilationResult
end

# make_stabilizer_list still takes num_pauli_qubits from the generator COUNT
# (length(s)) rather than nqubits(s), so a stabilizer with fewer generators
# than qubits still embeds wrong internally. That is no longer reachable
# through the normal compilation path, though: validate_input now rejects any
# input whose rank is below the qubit count (mixed/underdetermined states)
# before make_stabilizer_list ever sees it, replacing the old opaque
# QuantumClifford ArgumentError with a clear, documented one.
@testset "rank-deficient input states are rejected by validate_input" begin
    circuit = Circuit([
        ExpEighPiPauli(P"Z", [1]),
        Measurement(P"Z", 1, [1]),
        Measurement(P"Z", 2, [2]),
    ])
    full_rank = Stabilizer(S"ZI IZ")
    rank_deficient = Stabilizer(S"ZI")      # 1 generator, 2 qubits: qubit 2 unstabilized

    @test nqubits(rank_deficient) == get_circuit_width(circuit)
    # Full rank: generator count == qubit count, so the embed happens to line up
    @test length(make_stabilizer_list(full_rank, circuit)) == 2
    # Rank deficient: caught up front now, with a message naming the cause
    @test_throws ArgumentError validate_input(circuit, rank_deficient)
    @test_throws ArgumentError run(copy(circuit), DummyRuntime(), rank_deficient)
end

end
