@testitem "Bugfix regressions" tags=[:bugfixes] begin

using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, Pauli, Measurement, ExpHalfPiPauli, ExpQuatPiPauli,
    ExpEighPiPauli, PauliConditional, BitConditional, traversal, preprocess_circuit,
    validate_CircuitOp, PauliQubitMismatchError, find_variant_indices, group_nonclifford,
    quantum_measurement, resolve_conditionals, get_distribution, get_graph,
    weight_std_graph, to_result, run, CompilerState, CompilationResult, MeasurementResult,
    SimRuntime, DummyRuntime, StabilizerRuntime
using .MeasurementResult: ClassicalDetermRes
using QuantumClifford: @P_str, @S_str, MixedDestabilizer, nqubits, Stabilizer
using Moshi.Derive: @derive
using Moshi.Match: isa_variant
using Graphs: nv

@derive CircuitOp[Eq, Show]

@testset "traversal with end_index left of start is a no-op" begin
    circuit = Circuit([Pauli(P"X", [1]), Pauli(P"Y", [1]), Pauli(P"Z", [1])])
    snapshot = copy(circuit)
    traversal(circuit, (a, b) -> (b, a), :left, 1, 0)
    @test circuit == snapshot
    traversal(circuit, (a, b) -> (b, a), :right, 2, 1)
    @test circuit == snapshot
end

@testset "group_nonclifford with a non-Clifford at index 1" begin
    e1 = ExpEighPiPauli(P"Z", [1])
    c  = ExpQuatPiPauli(P"X", [1])
    e2 = ExpEighPiPauli(P"Z", [1])
    circuit = Circuit([e1, c, e2])
    group_nonclifford(circuit)
    # e1 has nothing to its left and must stay put; e2 commutes past c (conjugated)
    @test circuit[1] == e1
    @test isa_variant(circuit[2], CircuitOp.ExpEighPiPauli)
    @test circuit[3] == c
end

@testset "validate_CircuitOp on PauliConditional with empty qubits" begin
    op = PauliConditional(P"X", Int[], P"Z", [2])
    @test_throws PauliQubitMismatchError validate_CircuitOp(op)
end

@testset "preprocess_circuit on a measurement-free circuit" begin
    circuit = Circuit([ExpQuatPiPauli(P"X", [1]), ExpHalfPiPauli(P"Z", [1])])
    preprocess_circuit(circuit)
    @test length(circuit) == 2
end

@testset "gadget classical bits never collide with user bits" begin
    # User bit index (2) exceeds the qubit count (1); gadget bits must go above it
    circuit = Circuit([ExpEighPiPauli(P"Z", [1]), Measurement(P"Z", 2, [1])])
    preprocess_circuit(circuit)
    bits = [circuit[i].bit for i in find_variant_indices(circuit, Measurement)]
    @test allunique(bits)
    @test 2 in bits
end

@testset "quantum_measurement without magic state throws ArgumentError" begin
    op = Measurement(P"X", 1, [1])
    @test_throws ArgumentError quantum_measurement(SimRuntime(), op, 1)
end

@testset "resolve_conditionals resolves every determined conditional" begin
    circuit = Circuit([
        BitConditional(ExpHalfPiPauli(P"X", [1]), 1),
        BitConditional(ExpHalfPiPauli(P"X", [1]), 2),
        Measurement(P"Z", 3, [1]),
    ])
    state = CompilerState(
        measurement_results = MeasurementResult.Type[],
        stabilizer_group = MixedDestabilizer(S"Z"),
        classical_register = Union{Nothing,Bool}[true, false, nothing],
        circuit = circuit,
        instruction_pointer = 1,
        runtime = DummyRuntime(),
    )
    resolve_conditionals(state)
    @test isempty(find_variant_indices(state.circuit, BitConditional))
end

@testset "get_distribution performs exactly num_shots shots" begin
    circuit = Circuit([Measurement(P"Z", 1, [1])])
    (distribution, data) = get_distribution(circuit, DummyRuntime(), nothing, 50)
    @test sum(distribution) == 50
    @test length(data) == 50
    # Z measurement on the default |0> state is deterministic
    @test distribution == [50, 0]
end

@testset "show on partially filled CompilationResult" begin
    mres = Vector{MeasurementResult.Type}(undef, 2)
    mres[1] = ClassicalDetermRes(P"Z", true)
    r = CompilationResult(mres, MeasurementResult.Type[], S"Z", 0)
    @test occursin("undefined", sprint(show, r))
end

@testset "get_graph and weight_std_graph on a compiled circuit" begin
    circuit = Circuit([
        ExpQuatPiPauli(P"X", [1]),
        ExpEighPiPauli(P"Z", [2]),
        Measurement(P"Z", 1, [1]),
        Measurement(P"Z", 2, [2]),
    ])
    state = run(copy(circuit), DummyRuntime())
    result = to_result(state)
    g_all = get_graph(result, false)
    @test nv(g_all) == size(result.stabilizer_group, 2)
    g_qpu = get_graph(result, true)
    @test nv(g_qpu) >= 0  # smoke: quantum-only graph builds without error

    g_std = weight_std_graph(circuit, DummyRuntime(); num_shots = 5)
    @test nv(g_std) == nv(g_all)
end

@testset "quantum measurements act on absolute qubit positions" begin
    # A leading T gate keeps its partial qubit list ([2]) through preprocessing,
    # so its gadget measurement Pauli must be embedded before magic-qubit slicing;
    # unembedded, the slice silently read out-of-bounds as identity
    circuit = Circuit([
        ExpEighPiPauli(P"Z", [2]),
        ExpEighPiPauli(P"Z", [1]),
        Measurement(P"Z", 1, [1]),
        Measurement(P"Z", 2, [2]),
    ])
    state = run(copy(circuit), SimRuntime())
    width = nqubits(state.stabilizer_group)
    quantum = filter(mr -> isa_variant(mr, MeasurementResult.QuantumRes),
                     state.measurement_results)
    @test !isempty(quantum)
    @test all(length(mr.pauli) == width for mr in quantum)
    # Each gadget measurement must touch its magic qubit (never sliced to identity)
    magicqubits = width - length(quantum) + 1 : width
    for mr in quantum
        magic_p = mr.pauli[magicqubits]
        @test any(magic_p[i] != (false, false) for i in 1:length(magic_p))
    end
end

@testset "imaginary-phase Paulis are rejected" begin
    @test_throws ArgumentError validate_CircuitOp(ExpQuatPiPauli(im * P"X", [1]))
    @test_throws ArgumentError validate_CircuitOp(Measurement(-im * P"Z", 1, [1]))
    # Real phases remain valid
    validate_CircuitOp(ExpQuatPiPauli(-P"X", [1]))
    @test true
end

@testset "random_test_circuit only generates Hermitian Paulis" begin
    circuit = random_test_circuit(50, 3)
    for op in circuit
        if !isa_variant(op, CircuitOp.PauliConditional)
            @test iseven(op.pauli.phase[])
        end
    end
end

@testset "traversal deletes both ops on empty-tuple result" begin
    cancel_all(op1, op2) = ()
    circuit = Circuit([Pauli(P"X", [1]), Pauli(P"Y", [1])])
    traversal(circuit, cancel_all, :right)
    @test isempty(circuit)
    # :left direction, deletion at the tail must not read out of bounds
    circuit = Circuit([Pauli(P"X", [1]), Pauli(P"Y", [1]), Pauli(P"Z", [1])])
    delete_yz(op1, op2) = (op1.pauli == P"Y" && op2.pauli == P"Z") ? () : nothing
    traversal(circuit, delete_yz, :left)
    @test length(circuit) == 1 && circuit[1].pauli == P"X"
end

@testset "inverse rotations cancel instead of leaving a placeholder" begin
    merge_ops = PBCCompiler.merge_ops
    # pi/8 followed by -pi/8 about the same axis is the identity
    circuit = Circuit([ExpEighPiPauli(P"Z", [1]), ExpEighPiPauli(-P"Z", [1])])
    merge_ops(circuit)
    @test isempty(circuit)
    # two pi/2 rotations about the same axis are a global phase
    circuit = Circuit([ExpHalfPiPauli(P"Z", [1]), ExpHalfPiPauli(P"Z", [1])])
    merge_ops(circuit)
    @test isempty(circuit)
    # same-sign pi/8 pair still merges into a pi/4 rotation
    circuit = Circuit([ExpEighPiPauli(P"Z", [1]), ExpEighPiPauli(P"Z", [1])])
    merge_ops(circuit)
    @test length(circuit) == 1
    @test isa_variant(circuit[1], CircuitOp.ExpQuatPiPauli)
end

@testset "get_hypergraph keeps isolated qubits in the vertex set" begin
    circuit = Circuit([
        ExpQuatPiPauli(P"X", [1]),
        ExpEighPiPauli(P"Z", [2]),
        Measurement(P"Z", 1, [1]),
        Measurement(P"Z", 2, [2]),
    ])
    state = run(copy(circuit), DummyRuntime())
    result = to_result(state)
    (A, h) = PBCCompiler.get_hypergraph(result)
    @test size(A, 1) == size(result.stabilizer_group, 2)
end

@testset "input validation rejects malformed operations" begin
    validate_circuit = PBCCompiler.validate_circuit
    # PrepMagic is not supported by the pipeline
    @test_throws ArgumentError validate_CircuitOp(PBCCompiler.PrepMagic(2, [1]))
    # Overlapping control/target registers in PauliConditional
    @test_throws ArgumentError validate_CircuitOp(PauliConditional(P"Z", [1], P"X", [1]))
    # Non-positive qubit indices
    @test_throws ArgumentError validate_CircuitOp(Measurement(P"Z", 1, [0]))
    @test_throws ArgumentError validate_CircuitOp(Measurement(P"Z", 1, [-2]))
    # Non-positive classical bit indices
    @test_throws ArgumentError validate_CircuitOp(Measurement(P"Z", 0, [1]))
    @test_throws ArgumentError validate_CircuitOp(BitConditional(ExpHalfPiPauli(P"X", [1]), 0))
    # Duplicate measurement bits would silently drop a measurement downstream
    @test_throws ArgumentError validate_circuit(Circuit([
        Measurement(P"Z", 1, [1]), Measurement(P"X", 1, [1])]))
    # BitConditional controlled by a bit no measurement writes
    @test_throws ArgumentError validate_circuit(Circuit([
        BitConditional(ExpHalfPiPauli(P"X", [1]), 7), Measurement(P"Z", 1, [1])]))
    # A well-formed circuit still validates
    validate_circuit(Circuit([
        BitConditional(ExpHalfPiPauli(P"X", [1]), 1),
        PauliConditional(P"Z", [1], P"X", [2]),
        Measurement(P"Z", 1, [1]), Measurement(P"Z", 2, [2])]))
    @test true
    # ... and run() surfaces the validation error
    @test_throws ArgumentError PBCCompiler.run(
        Circuit([Measurement(P"Z", 1, [1]), Measurement(P"X", 1, [1])]), DummyRuntime())
end

@testset "single-op circuits go through the full pipeline" begin
    circuit = Circuit([ExpEighPiPauli(P"Z", [1])])
    preprocess_circuit(circuit)
    @test isempty(find_variant_indices(circuit, ExpEighPiPauli))
    @test !isempty(find_variant_indices(circuit, Measurement))

    circuit = Circuit([PauliConditional(P"Z", [1], P"X", [2])])
    preprocess_circuit(circuit)
    @test isempty(find_variant_indices(circuit, PauliConditional))
end

@testset "parse_QuantumClifford expands .inc gate definitions" begin
    using QuantumClifford: sT, sCNOT, sHadamard, sInvPhase
    fixture = joinpath(@__DIR__, "fixtures", "toffoli3.qasm")
    qc_circuit = PBCCompiler.parse_QuantumClifford(fixture)
    # Standard qelib1 ccx: 2 h, 6 cx, 4 t, 3 tdg (each tdg -> sT + sInvPhase)
    @test PBCCompiler.get_T_count(qc_circuit) == 7
    @test count(op -> op isa sCNOT, qc_circuit) == 6
    @test count(op -> op isa sHadamard, qc_circuit) == 2
    @test count(op -> op isa sInvPhase, qc_circuit) == 3
    # Unknown gates now error instead of being silently dropped
    mktemp() do path, io
        write(io, "OPENQASM 2.0;\nqreg q[2];\nmystery q[0],q[1];\n"); close(io)
        @test_throws ErrorException PBCCompiler.parse_QuantumClifford(path)
    end
end

@testset "circuit that cancels to empty after preprocessing" begin
    circuit = Circuit([ExpEighPiPauli(P"Z", [1]), ExpEighPiPauli(-P"Z", [1])])
    state = run(circuit, SimRuntime())
    @test isempty(state.classical_register)
    @test isempty(state.measurement_results)
end

@testset "magic-state count matches gadgets in the executed circuit" begin
    # Leading and trailing T gates: exactly two gadgets, two magic qubits
    circuit = Circuit([
        ExpEighPiPauli(P"Z", [2]),
        ExpQuatPiPauli(P"X", [1]),
        ExpEighPiPauli(P"Z", [1]),
        Measurement(P"Z", 1, [1]),
        Measurement(P"Z", 2, [2]),
    ])
    state = run(copy(circuit), SimRuntime())
    input_width = 2
    @test nqubits(state.stabilizer_group) - input_width == 2
end

@testset "deferred T fires when a gadget's support is not the full register" begin
    # A pi/8 rotation with no Clifford ahead of it reaches `gadgetize` with its
    # original (partial) qubit list, so the gadget measurement's support is
    # [1, magic] rather than the full 1:width. The lazy-T activation must still
    # find the magic qubit. Selecting it by position within `op.qubits` instead
    # of by value yields an empty selection here and silently skips the T,
    # leaving the magic qubit in |+>. Circuits whose rotations all get conjugated
    # to full width (toffoli3, adder_n4) mask this: there position and value
    # coincide, so those fixtures cannot catch it.
    input = Stabilizer([P"X_", P"_Z"])          # |+>|0>
    circuit() = Circuit([ExpEighPiPauli(P"Z", [1]),
                         Measurement(P"X", 1, [1]),
                         Measurement(P"Z", 2, [2])])

    for rt in (SimRuntime, StabilizerRuntime)
        state = run(circuit(), rt(), copy(input))
        @test all(state.runtime.activated)
    end

    # T|+> measured along X gives P(-1) = sin^2(pi/8) ~ 0.1464. With the T
    # skipped the magic qubit stays |+> and the frequency comes out ~0.25, far
    # outside the tolerance below (~5 sigma at this shot count).
    nshots = 1000
    for rt in (SimRuntime, StabilizerRuntime)
        ones = count(_ -> run(circuit(), rt(), copy(input)).classical_register[1] === true,
                     1:nshots)
        @test isapprox(ones / nshots, sin(pi/8)^2; atol = 0.06)
    end
end

@testset "gadgetize keeps Pauli letters attached to their qubits" begin
    # `paulis` is positional in `op.qubits`; `affectedqubits` sorts. For an
    # unsorted qubit list the two disagree, and the gadget silently measures the
    # wrong observable. Z on qubit 2 and X on qubit 1 must stay that way.
    op = ExpEighPiPauli(P"ZX", [2, 1])
    gadget = PBCCompiler.gadgetize(op, 2, 2, 1)

    meas = gadget[1]
    @test meas.qubits == [1, 2, 3]
    @test meas.pauli == P"XZZ"          # X on q1, Z on q2, Z on the magic qubit

    # The conditional corrections carry the same Pauli and must agree with it
    for correction in (gadget[3], gadget[4])
        @test correction.op.qubits == [1, 2]
        @test correction.op.pauli == P"XZ"
    end

    # Sorted input is unaffected
    sorted_gadget = PBCCompiler.gadgetize(ExpEighPiPauli(P"ZX", [1, 2]), 2, 2, 1)
    @test sorted_gadget[1].pauli == P"ZXZ"

    # A negative phase must survive the permutation
    @test PBCCompiler.gadgetize(ExpEighPiPauli(-P"ZX", [2, 1]), 2, 2, 1)[1].pauli == -P"XZZ"
end

@testset "unsorted-qubit pi/8 rotation compiles to the right observable" begin
    # exp(-i*pi/8*P) with P = X on qubit 1, Z on qubit 2, written unsorted as
    # P"ZX" on [2,1]. On |00> that rotates within span{|00>, |10>}, so bit 1
    # picks up P(1) = sin^2(pi/8) while bit 2 stays 0. Under the letter/qubit
    # mismatch the observable becomes Z1*X2 instead, which rotates within
    # span{|00>, |01>} -- the randomness moves to bit 2 and bit 1 goes
    # deterministic, so both assertions below flip.
    circuit() = Circuit([ExpEighPiPauli(P"ZX", [2, 1]),
                         Measurement(P"Z", 1, [1]),
                         Measurement(P"Z", 2, [2])])
    nshots = 1000
    for rt in (SimRuntime, StabilizerRuntime)
        states = [run(circuit(), rt()) for _ in 1:nshots]
        ones = count(s -> s.classical_register[1] === true, states)
        @test isapprox(ones / nshots, sin(pi/8)^2; atol = 0.06)
        @test all(s -> s.classical_register[2] === false, states)
    end
end

@testset "input states must be fully stabilized" begin
    # The gadget path factors each measurement into a data part and a magic part
    # and needs the data part to have a definite eigenvalue, which only holds at
    # full rank. Underdetermined states used to get through validation and fail
    # much later -- either inside `embed` while padding the generators, or with
    # "no definite eigenvalue" from within a gadget measurement.
    circuit = Circuit([ExpEighPiPauli(P"Z", [2]),
                       Measurement(P"X", 1, [2]),
                       Measurement(P"Z", 2, [1])])

    # Too few generators
    @test_throws ArgumentError run(copy(circuit), SimRuntime(), Stabilizer([P"Z_"]))
    # Right number of generators, but linearly dependent -> still rank 1
    @test_throws ArgumentError run(copy(circuit), SimRuntime(), Stabilizer([P"Z_", P"Z_"]))
    # Full rank is accepted
    @test run(copy(circuit), SimRuntime(), Stabilizer([P"Z_", P"_X"])) isa CompilerState

    # The message must name the cause rather than surfacing an embed/eigenvalue error
    err = try
        run(copy(circuit), SimRuntime(), Stabilizer([P"Z_"]))
    catch e
        e
    end
    @test occursin("fully stabilized", err.msg)

    # The width check still fires, and no longer needs a non-empty tableau to run
    @test_throws ArgumentError PBCCompiler.validate_input(
        Circuit([Measurement(P"Z", 1, [1])]), Stabilizer([P"Z_", P"_Z"]))
end

@testset "get_distribution on circuits with no classical bits" begin
    (distribution, data) = get_distribution(Circuit(), DummyRuntime(), nothing, 3)
    @test distribution == [3]
    @test data == [0, 0, 0]
end

end
