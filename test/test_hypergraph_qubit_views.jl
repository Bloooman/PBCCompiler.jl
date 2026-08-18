@testitem "get_hypergraph/get_graph qubits and variant views" tags=[:statistics] begin

using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, ExpEighPiPauli, Measurement, DummyRuntime,
    DummyStabilizerRuntime, run, to_result, get_graph, get_hypergraph, weight_std_graph
using QuantumClifford: @P_str, nqubits
using Graphs: nv

# 2 data qubits, one T gate -> one magic-state gadget, so the register is
# data(2) + magic(>=1). `state.stabilizer_group` (pre-`to_result`) is always
# full-register-width regardless of runtime -- unlike the `CompilationResult`
# it produces, whose tableau is sliced down to the data qubits for
# `AbstractStabilizerRuntime` results (see `to_result` in `logic.jl`) -- so use
# it, not `result.stabilizer_group`, to compute the expected register width.
circuit() = Circuit(CircuitOp.Type[
    CircuitOp.ExpEighPiPauli(P"Z", [1]),
    CircuitOp.Measurement(P"Z", 1, [1]),
    CircuitOp.Measurement(P"Z", 2, [2]),
])
const N_INPUT = 2

@testset "qubits=:magic sizes to the magic block, not the full register, for AbstractStabilizerRuntime" begin
    state = run(copy(circuit()), DummyStabilizerRuntime())
    register_n = Int(nqubits(state.stabilizer_group))
    magic_n = register_n - N_INPUT
    @test magic_n > 0  # sanity: there is a magic block to test the fix against

    result = to_result(state)
    # The bug being fixed: result.stabilizer_group is sliced to data-only for
    # AbstractStabilizerRuntime, so this must NOT be mistaken for register_n.
    @test size(result.stabilizer_group, 2) == N_INPUT

    g = get_graph(result; qubits=:magic, n_input=N_INPUT)
    @test nv(g) == magic_n

    (A, h) = get_hypergraph(result; qubits=:magic, n_input=N_INPUT)
    @test size(A, 1) == magic_n
    @test Int(h.n_vertices) == magic_n
end

@testset "qubits=:all sizes to the full register for both runtime families" begin
    for rt in (DummyRuntime(), DummyStabilizerRuntime())
        state = run(copy(circuit()), rt)
        register_n = Int(nqubits(state.stabilizer_group))
        result = to_result(state)

        g = get_graph(result; qubits=:all)
        @test nv(g) == register_n

        (A, h) = get_hypergraph(result; qubits=:all)
        @test size(A, 1) == register_n
    end
end

@testset "qubits=:data sizes to the data block alone, dropping the magic block" begin
    state = run(copy(circuit()), DummyStabilizerRuntime())
    result = to_result(state)

    g = get_graph(result; qubits=:data, n_input=N_INPUT)
    @test nv(g) == N_INPUT

    (A, h) = get_hypergraph(result; qubits=:data, n_input=N_INPUT)
    @test size(A, 1) == N_INPUT
    @test Int(h.n_vertices) == N_INPUT
end

@testset ":data/:magic without n_input raise ArgumentError" begin
    state = run(copy(circuit()), DummyRuntime())
    result = to_result(state)
    @test_throws ArgumentError get_hypergraph(result; qubits=:magic)
    @test_throws ArgumentError get_graph(result; qubits=:data)
end

@testset "variant filter composes with the qubits view" begin
    state = run(copy(circuit()), DummyStabilizerRuntime())
    result = to_result(state)
    (A_all, _) = get_hypergraph(result; qubits=:magic, n_input=N_INPUT, variant=:all)
    (A_q, _) = get_hypergraph(result; qubits=:magic, n_input=N_INPUT, variant=:quantum)
    @test size(A_q, 2) <= size(A_all, 2)
end

@testset "weight_std_graph runs end-to-end with the new kwargs" begin
    state = run(copy(circuit()), DummyStabilizerRuntime())
    register_n = Int(nqubits(state.stabilizer_group))
    g = weight_std_graph(circuit(), DummyStabilizerRuntime(); qubits=:magic, n_input=N_INPUT, num_shots=3)
    @test nv(g) == register_n - N_INPUT
end

end
