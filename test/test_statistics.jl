@testitem "Statistics graph/hypergraph utilities" tags=[:statistics] begin

using PBCCompiler
using PBCCompiler: CompilationResult, MeasurementResult, find_depth
using .MeasurementResult: QuantumRes
using SparseArrays
using KaHyPar
using Graphs: add_edge!, nv, ne
using SimpleWeightedGraphs: SimpleWeightedGraph, get_weight
using QuantumClifford: @P_str, @S_str, PauliOperator, nqubits
using PBCCompiler: Circuit, CircuitOp, ExpEighPiPauli, Measurement, DummyRuntime,
    DummyStabilizerRuntime, run, to_result

# Two tiny hypergraphs over 3 vertices (incidence matrices are 1-based;
# KaHyPar stores vertices 0-based internally):
#   h1: edges {1,2} and {2,3}
#   h2: edge {1,2}
A1 = sparse([1, 2, 2, 3], [1, 1, 2, 2], ones(Int, 4), 3, 2)
h1 = KaHyPar.HyperGraph(A1, ones(Int, 3), [1, 1])
A2 = sparse([1, 2], [1, 1], ones(Int, 2), 3, 1)
h2 = KaHyPar.HyperGraph(A2, ones(Int, 3), [1])

@testset "hyperedge_frequency" begin
    freq = PBCCompiler.hyperedge_frequency([h1, h2])
    @test freq == Dict([1, 2] => 2, [2, 3] => 1)
end

@testset "variant_hypergraph" begin
    vh = PBCCompiler.variant_hypergraph([h1, h2])
    # {1,2} appears in both graphs (weight +1); {2,3} only in h1 (weight -1)
    @test Int(vh.n_vertices) == 3
    @test sort(Int.(vh.e_weights)) == [-1, 1]
    verts_by_weight = Dict(
        sort(Int.(vh.hyperedges[vh.edge_indices[j]+1:vh.edge_indices[j+1]]) .+ 1) => Int(vh.e_weights[j])
        for j in 1:length(vh.edge_indices)-1)
    @test verts_by_weight == Dict([1, 2] => 1, [2, 3] => -1)
end

@testset "hyperedge_size_distribution" begin
    @test PBCCompiler.hyperedge_size_distribution(h1) == Dict(2 => 2)
end

@testset "HyperedgeCut and HyperedgeConnectivity" begin
    # parts[v+1] is the block of 0-based vertex v: blocks (0,0,1)
    parts = Int64[0, 0, 1]
    # {1,2} uncut; {2,3} spans blocks 0 and 1
    @test PBCCompiler.HyperedgeCut(h1, parts) == 1
    @test PBCCompiler.HyperedgeConnectivity(h1, parts) == 1
    all_one_block = Int64[0, 0, 0]
    @test PBCCompiler.HyperedgeCut(h1, all_one_block) == 0
    @test PBCCompiler.HyperedgeConnectivity(h1, all_one_block) == 0
end

@testset "variant_graph" begin
    g1 = SimpleWeightedGraph(3)
    add_edge!(g1, 1, 2, 1); add_edge!(g1, 2, 3, 1)
    g2 = SimpleWeightedGraph(3)
    add_edge!(g2, 1, 2, 5)
    vg = PBCCompiler.variant_graph([g1, g2])
    @test nv(vg) == 3 && ne(vg) == 2
    @test get_weight(vg, 1, 2) == 1    # present in both
    @test get_weight(vg, 2, 3) == -1   # missing from g2
end

@testset "find_depth" begin
    # find_depth only reads result.QPU_workload; stabilizer_group and
    # measurement_results are irrelevant, so a minimal dummy suffices.
    make_result(paulis) = CompilationResult(
        MeasurementResult.Type[],
        MeasurementResult.Type[QuantumRes(p, true) for p in paulis],
        S"Z",
        length(paulis),
    )

    @testset "empty QPU_workload" begin
        @test find_depth(make_result(PauliOperator[])) == []
    end

    @testset "single measurement is not dropped" begin
        layers = find_depth(make_result([P"XII"]))
        @test layers == [[Set([1])]]
    end

    @testset "all-disjoint measurements merge into one layer (previously infinite-looped)" begin
        layers = find_depth(make_result([P"XII", P"IXI", P"IIX"]))
        @test length(layers) == 1
        @test layers[1] == [Set([1]), Set([2]), Set([3])]
    end

    @testset "all-mutually-conflicting measurements each get their own layer, none dropped" begin
        layers = find_depth(make_result([P"XX", P"XX", P"XX"]))
        @test layers == [[Set([1, 2])], [Set([1, 2])], [Set([1, 2])]]
    end

    @testset "mixed: disjoint pair followed by a conflicting one starts a new layer" begin
        layers = find_depth(make_result([P"XII", P"IXI", P"XII"]))
        @test length(layers) == 2
        @test layers[1] == [Set([1]), Set([2])]
        @test layers[2] == [Set([1])]
    end

    @testset "layer count never exceeds QPUDuration" begin
        result = make_result([P"XII", P"IXI", P"XII", P"IIX"])
        @test length(find_depth(result)) <= result.QPUDuration
    end

    @testset "terminates and stays within QPUDuration on a real compiled circuit" begin
        circuit() = Circuit(CircuitOp.Type[
            ExpEighPiPauli(P"Z", [1]),
            Measurement(P"Z", 1, [1]),
            Measurement(P"Z", 2, [2]),
        ])
        for rt in (DummyRuntime(), DummyStabilizerRuntime())
            state = run(copy(circuit()), rt)
            result = to_result(state)
            layers = find_depth(result)
            @test length(layers) <= result.QPUDuration
            @test sum(length, layers; init=0) == length(result.QPU_workload)
        end
    end
end

end
