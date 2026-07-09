#This file contains functions for analyzing computed&compiled circuits
##
using Graphs, SimpleWeightedGraphs, KaHyPar, SparseArrays, LinearAlgebra
using QuantumClifford: xbit, zbit
using Statistics
using StatsBase

"""
    get_distribution(input_circuit::Circuit, input_state::Stabilizer, num_shots=1000)

Get measurement result distribution of given input circuit and input state by running compute&compile function for desired number of shots

Return calculated result distribution and raw result count
"""
function get_distribution(input_circuit::Circuit, rt::R, input_state::Union{Stabilizer, Nothing}=nothing, num_shots::Int=1000) where R<:AbstractRuntime
    num_bits = get_circuit_width(input_circuit)
    len = 2^num_bits
    distribution = zeros(Int, len)
    data = zeros(Int, num_shots)
    i=1
    while i<num_shots
        circuit = copy(input_circuit)
        result_i=run(circuit, rt, input_state)
        register=result_i.classical_register
        final_measurement_results=register[1:num_bits]
        bit_str = join(Int.(final_measurement_results))
        index = parse(Int, bit_str; base=2) + 1
        distribution[index]+=1
        data[i]=index-1
        i+=1
    end
    return (distribution, data)
end

"""
    get_graph(result::S, type=nothing) where S <: AbstractRuntime

Extract qubit interaction graph from result.compiler_state

When quantum_only is false, plot interaction hypergraph among all qubits using all Pauli Product Measurement
When quantum_only is true, plot interaction hypergraph of injected qubits that actually live on QPU
"""
function get_graph(result::CompilationResult, quantum_only::Bool=false)
    num_nodes=length(result.compiler_state.stabilizer_group[1])
    g=SimpleWeightedGraph{Int64, Int64}(Int64(num_nodes))
    measurements = quantum_only ? result.measurement_results : result.QPU_workload
    for m in measurements
        p=m.pauli
        for i in 1:length(p)
            @debug(i)
            if xbit(p)[i] || zbit(p)[i]
                for j in i+1:length(p)
                    @debug(j)
                    if xbit(p)[j] || zbit(p)[j]
                        current_w = get_weight(g, i, j)
                        add_edge!(g,i,j,current_w + 1)
                    end
                end
            end
        end
    end
    if quantum_only
        for h in 1:length(1 - result.compiler_state.num_gadgets)
            rem_vertex!(g,1)
        end
    end
    return g
end

"""
    variant_graph(graphs::Vector{<:SimpleWeightedGraph}) -> SimpleWeightedGraph

Merge all edges from `graphs` into a single weighted graph.

Hyperedges are identified by their vertex set (order-independent). An edge
that appears in every input graph receives weight `+1`; an edge absent
from at least one graph receives weight `-1`.

# Arguments
- `graphs`: non-empty vector of graphs sharing the same vertex count

# Returns
A `SimpleWeightedGraph` over the same vertex set containing the union of all
unique edges, with `±1` edge weights as described above.
"""
function variant_graph(graphs::Vector{<:SimpleWeightedGraph})::SimpleWeightedGraph
    N = length(graphs)
    N == 0 && throw(ArgumentError("Input vector must be non-empty"))

    n = nv(graphs[1])

    edge_count = Dict{Tuple{Int,Int}, Int}()
    for g in graphs
        for e in edges(g)
            key = minmax(Int(src(e)), Int(dst(e)))
            edge_count[key] = get(edge_count, key, 0) + 1
        end
    end

    I_vals = Int[]
    J_vals = Int[]
    V_vals = Int[]
    for ((u, v), count) in edge_count
        w = count == N ? 1 : -1
        push!(I_vals, u, v)
        push!(J_vals, v, u)
        push!(V_vals, w, w)
    end

    s = isempty(I_vals) ? spzeros(Int, n, n) :
                          sparse(I_vals, J_vals, V_vals, n, n)
    return SimpleWeightedGraph(s)
end
##
"""
    weight_std_graph(graphs::Vector{<:SimpleWeightedGraph})

Given a list of `SimpleWeightedGraph` objects with the same number of vertices,
return a new `SimpleWeightedGraph` whose edge weights are the standard deviation
of each edge's weight across the input graphs.

Edges absent from a graph in the list contribute weight 0 to the computation.
"""
function weight_std_graph(input_circuit::Circuit, input_state::Stabilizer; quantum_only::Bool=false, num_shots::Int=1000)
    graphs =  Vector{SimpleWeightedGraph}(undef, num_shots)
    i=1
    while i<num_shots+1
        circuit = copy(input_circuit)
        result_i=run(circuit, input_state)
        graphs[i]=get_graph(result_i, quantum_only)
        i+=1
    end
    n = nv(first(graphs))
    all_edges = Set{Tuple{Int,Int}}()
    for g in graphs
        for e in edges(g)
            push!(all_edges, (min(src(e), dst(e)), max(src(e), dst(e))))
        end
    end
    result = SimpleWeightedGraph(n)
    for (u, v) in all_edges
        ws = [get_weight(g, u, v) for g in graphs]
        add_edge!(result, u, v, std(ws))
    end
    return result
end
##
"""
    get_hypergraph(result::S, quantum_only=false) where S<:AbstractRuntime

Extract qubit interaction hypergraph from resulted CompilerState.

When quantum_only is false, plot interaction hypergraph among all qubits using all Pauli Product Measurement
When quantum_only is true, plot interaction hypergraph of injected qubits that actually live on QPU
"""
function get_hypergraph(result::CompilationResult, quantum_only::Bool=false)
    measurements = quantum_only ? result.QPU_workload : result.measurement_results
    paulis=[m.pauli for m in measurements]
    collected_edges = Vector{Vector{Int}}()
    for p in paulis
        push!(collected_edges, _qubit_coverage(p))
    end
    w=countmap(collected_edges)
    I = Int[]
    J = Int[]
    i=1
    for row in keys(w)
        append!(I,row)
        col=fill(i,length(row))
        append!(J,col)
        i+=1
    end
    V = Int.(ones(length(I)))
    edge_weights=collect(values(w))
    A = sparse(I, J, V)
    num_v, num_e = size(A)
    h = KaHyPar.HyperGraph(A,ones(Int, num_v),edge_weights)
    return (A, h)
end

function _qubit_coverage(p::PauliOperator)
    bool_vec = [p[i] for i in 1:nqubits(p)]
    idx=findall(x -> x !== (false,false), bool_vec)
    return idx
end

##
"""
    edgeweight_variant_hypergraph(hypergraphs::Vector{KaHyPar.HyperGraph}) -> KaHyPar.HyperGraph

Merge all hyperedges from `hypergraphs` into a single weighted hypergraph.

Hyperedges are identified by their vertex set (order-independent). An edge
that appears in every input hypergraph receives weight `+1`; an edge absent
from at least one hypergraph receives weight `-1`.

# Arguments
- `hypergraphs`: non-empty vector of hypergraphs sharing the same vertex count

# Returns
A `KaHyPar.HyperGraph` over the same vertex set containing the union of all
unique hyperedges, with `±1` edge weights as described above.
"""
function variant_hypergraph(hypergraphs::Vector{KaHyPar.HyperGraph})::KaHyPar.HyperGraph
    N = length(hypergraphs)
    N == 0 && throw(ArgumentError("Input vector must be non-empty"))

    n_vertices = Int(hypergraphs[1].n_vertices)

    # Count occurrences of each unique hyperedge (normalized to sorted 0-based vertex vector)
    edge_count = Dict{Vector{Int}, Int}()
    for hg in hypergraphs
        n_edges = length(hg.edge_indices) - 1
        for j in 1:n_edges
            verts = sort(Int.(hg.hyperedges[hg.edge_indices[j]+1 : hg.edge_indices[j+1]]))
            edge_count[verts] = get(edge_count, verts, 0) + 1
        end
    end

    # Build sparse incidence matrix (n_vertices × n_unique_edges)
    unique_edges = collect(keys(edge_count))
    n_out_edges = length(unique_edges)
    I_vals = Int[]
    J_vals = Int[]
    for (j, verts) in enumerate(unique_edges)
        for v in verts
            push!(I_vals, v + 1)  # 0-based vertex index → 1-based row
            push!(J_vals, j)
        end
    end
    A = sparse(I_vals, J_vals, ones(Int, length(I_vals)), n_vertices, n_out_edges)

    vertex_weights = ones(Int, n_vertices)
    edge_weights = [edge_count[e] == N ? 1 : -1 for e in unique_edges]

    return KaHyPar.HyperGraph(A, vertex_weights, edge_weights)
end
##
"""
    HyperedgeCut(h::KaHyPar.HyperGraph, parts::Vector{Int64}) -> Int64

Count the number of cut hyperedges in a partitioned hypergraph.

# Arguments
- `h`: a KaHyPar hypergraph with CSR-encoded edge structure
- `parts`: partition assignment vector (1-indexed); `parts[v+1]` gives the
  0-based block ID of 0-indexed vertex `v`; length must equal `h.n_vertices`

# Returns
Number of hyperedges whose vertices span more than one partition block.
"""
function HyperedgeCut(h::KaHyPar.HyperGraph, parts::Vector{Int64})::Int64
    n_edges = length(h.edge_indices) - 1
    count = 0
    for e in 1:n_edges
        start = Int(h.edge_indices[e]) + 1     # 0-based C offset → 1-based Julia index
        stop  = Int(h.edge_indices[e + 1])     # 0-based exclusive end = 1-based inclusive end
        first_part = parts[Int(h.hyperedges[start]) + 1]
        is_cut = any(parts[Int(h.hyperedges[i]) + 1] != first_part for i in (start + 1):stop)
        count += is_cut
    end
    return count
end

"""
    edgecut(g::SimpleWeightedGraph, parts::Vector{Int64}) -> Int64

Count the number of cut edges in a partitioned weighted graph.

# Arguments
- `g`: an undirected weighted graph; vertices are 1-indexed
- `parts`: partition assignment vector (1-indexed); values are 0-based block IDs;
  length must equal the number of vertices in `g`

# Returns
Number of edges whose endpoints belong to different partition blocks.
"""
function edgecut(g::SimpleWeightedGraph, parts::Vector{Int64})::Int64
    return count(e -> parts[src(e)] != parts[dst(e)], edges(g))
end
##
"""
    hyperedge_frequency(hypergraphs::Vector{KaHyPar.HyperGraph}) -> Dict{Vector{Int}, Int}

Count how many times each unique hyperedge appears across all input hypergraphs.

# Arguments
- `hypergraphs`: vector of `KaHyPar.HyperGraph` instances to analyze

# Returns
A `Dict` mapping each unique hyperedge — represented as a sorted, 1-indexed vertex
vector — to its total occurrence count across all input hypergraphs.
"""
function hyperedge_frequency(hypergraphs::Vector{KaHyPar.HyperGraph})::Dict{Vector{Int}, Int}
    counts = Dict{Vector{Int}, Int}()
    for hg in hypergraphs
        n_edges = length(hg.edge_indices) - 1
        for j in 1:n_edges
            verts = sort!(Int.(hg.hyperedges[hg.edge_indices[j]+1 : hg.edge_indices[j+1]]) .+ 1)
            counts[verts] = get(counts, verts, 0) + 1
        end
    end
    return counts
end

function num_nonclifford(circuit::Circuit)
    num = length(find_variant_indices(circuit, ExpEighPiPauli))
    return num
end
##
"""
    hyperedge_size_distribution(h::KaHyPar.HyperGraph) -> Dict{Int,Int}

Count the number of hyperedges at each size across the entire hypergraph.

# Arguments
- `h`: the hypergraph to analyse

# Returns
A `Dict` mapping each hyperedge size to its count.
"""
function hyperedge_size_distribution(h::KaHyPar.HyperGraph)::Dict{Int,Int}
    dist = Dict{Int,Int}()
    for i in 1:(length(h.edge_indices) - 1)
        s = Int(h.edge_indices[i + 1] - h.edge_indices[i])
        dist[s] = get(dist, s, 0) + 1
    end
    return dist
end
##
"""
    HyperedgeConnectivity(h::KaHyPar.HyperGraph, parts::Vector{Int64}) -> Int64

Compute the connectivity (λ-1) metric: sum over hyperedges of
(number of distinct partition blocks spanned - 1).

# Arguments
- `h`: a KaHyPar hypergraph with CSR-encoded edge structure
- `parts`: partition assignment vector (1-indexed); `parts[v+1]` gives the
  0-based block ID of 0-indexed vertex `v`; length must equal `h.n_vertices`

# Returns
Sum of (λ_e - 1) over all hyperedges e, where λ_e is the number of blocks
that hyperedge e touches. Equals `HyperedgeCut` only when every cut edge
spans exactly 2 blocks.
"""
function HyperedgeConnectivity(h::KaHyPar.HyperGraph, parts::Vector{Int64})::Int64
    n_edges = length(h.edge_indices) - 1
    seen = Set{Int64}()
    total = 0
    for e in 1:n_edges
        start = Int(h.edge_indices[e]) + 1
        stop  = Int(h.edge_indices[e + 1])
        empty!(seen)
        for i in start:stop
            push!(seen, parts[Int(h.hyperedges[i]) + 1])
        end
        total += length(seen) - 1
    end
    return total
end
