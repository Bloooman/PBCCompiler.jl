@testitem "Statistics graph/hypergraph utilities" tags=[:statistics] begin

using PBCCompiler
using SparseArrays
using KaHyPar
using Graphs: add_edge!, nv, ne
using SimpleWeightedGraphs: SimpleWeightedGraph, get_weight

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

end
