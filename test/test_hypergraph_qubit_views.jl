@testitem "get_hypergraph/get_graph qubits and variant views" tags=[:statistics] begin

using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, ExpEighPiPauli, Measurement, DummyRuntime,
    DummyStabilizerRuntime, run, to_result, get_graph, get_hypergraph, weight_std_graph
using QuantumClifford: @P_str, nqubits
using Graphs: nv

# 2 input qubits, one T gate -> one gadget, so the register is
# input(2) + gadget(>=1). `state.stabilizer_group` (pre-`to_result`) is always
# full-register-width regardless of runtime -- unlike the `CompilationResult`
# it produces, whose tableau is sliced down to the input qubits for
# `AbstractStabilizerRuntime` results (see `to_result` in `logic.jl`) -- so use
# it, not `result.stabilizer_group`, to compute the expected register width.
circuit() = Circuit(CircuitOp.Type[
    CircuitOp.ExpEighPiPauli(P"Z", [1]),
    CircuitOp.Measurement(P"Z", 1, [1]),
    CircuitOp.Measurement(P"Z", 2, [2]),
])
const N_INPUT = 2

@testset "qubits=:gadget sizes to the gadget block, not the full register, for AbstractStabilizerRuntime" begin
    state = run(copy(circuit()), DummyStabilizerRuntime())
    register_n = Int(nqubits(state.stabilizer_group))
    gadget_n = register_n - N_INPUT
    @test gadget_n > 0  # sanity: there is a gadget block to test the fix against

    result = to_result(state)
    # The bug being fixed: result.stabilizer_group is sliced to input-only for
    # AbstractStabilizerRuntime, so this must NOT be mistaken for register_n.
    @test size(result.stabilizer_group, 2) == N_INPUT

    g = get_graph(result; qubits=:gadget, n_input=N_INPUT)
    @test nv(g) == gadget_n
end

@testset "get_hypergraph on AbstractStabilizerRuntime keeps only the input-qubit block" begin
    # get_hypergraph has no qubits=/n_input= kwargs (unlike get_graph): it
    # always sources hyperedges from QPU_workload, capped at
    # num_input = nqubits(result.stabilizer_group). For AbstractStabilizerRuntime
    # that field is sliced to the input qubits alone, so the cap intentionally
    # drops the gadget-qubit part of each measurement's support.
    state = run(copy(circuit()), DummyStabilizerRuntime())
    result = to_result(state)
    @test size(result.stabilizer_group, 2) == N_INPUT

    (A, h) = get_hypergraph(result)
    @test size(A, 1) == N_INPUT
    @test Int(h.n_vertices) == N_INPUT
end

@testset "get_graph qubits=:all sizes to the full register for both runtime families" begin
    for rt in (DummyRuntime(), DummyStabilizerRuntime())
        state = run(copy(circuit()), rt)
        register_n = Int(nqubits(state.stabilizer_group))
        result = to_result(state)

        g = get_graph(result; qubits=:all)
        @test nv(g) == register_n
    end
end

@testset "get_hypergraph on SimRuntime/DummyRuntime sizes to the full register (cap is a no-op)" begin
    # For this runtime family, QPU_workload Paulis are already restricted to
    # (and locally re-indexed over) the gadget-qubit block, while
    # stabilizer_group is full-register width -- so num_input is the full
    # register size and the `x <= num_input` cap never drops anything.
    for rt in (DummyRuntime(),)
        state = run(copy(circuit()), rt)
        register_n = Int(nqubits(state.stabilizer_group))
        result = to_result(state)

        (A, h) = get_hypergraph(result)
        @test size(A, 1) == register_n
    end
end

@testset "get_graph qubits=:input sizes to the input block alone, dropping the gadget block" begin
    state = run(copy(circuit()), DummyStabilizerRuntime())
    result = to_result(state)

    g = get_graph(result; qubits=:input, n_input=N_INPUT)
    @test nv(g) == N_INPUT
end

@testset "get_graph :input/:gadget without n_input raise ArgumentError" begin
    state = run(copy(circuit()), DummyRuntime())
    result = to_result(state)
    @test_throws ArgumentError get_graph(result; qubits=:input)
end

@testset "get_hypergraph on an empty QPU_workload yields an empty edge set, no error" begin
    # A circuit with no measurements at all has nothing to gate the gadget
    # register through, so QPU_workload is empty and every collected_edges
    # row is empty -- I/J/V stay empty, exercising the `sparse(..., num_input,
    # i - 1)` sizing on the zero-edge boundary rather than the removed
    # `maximum(I)`/`maximum(J)`, which throws on an empty collection.
    empty_circuit = Circuit(CircuitOp.Type[])
    state = run(copy(empty_circuit), DummyRuntime())
    result = to_result(state)
    @test isempty(result.QPU_workload)

    (A, h) = get_hypergraph(result)
    @test size(A) == (nqubits(result.stabilizer_group), 0)
end

@testset "weight_std_graph runs end-to-end with the new kwargs" begin
    state = run(copy(circuit()), DummyStabilizerRuntime())
    register_n = Int(nqubits(state.stabilizer_group))
    g = weight_std_graph(circuit(), DummyStabilizerRuntime(); qubits=:gadget, n_input=N_INPUT, num_shots=3)
    @test nv(g) == register_n - N_INPUT
end

end
