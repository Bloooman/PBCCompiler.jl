##
using Graphs, SimpleWeightedGraphs, KaHyPar, SparseArrays, LinearAlgebra
using Statistics
using StatsBase

"""
    get_distribution(input_circuit::Circuit, input_state::Stabilizer, num_shots=1000)

Get measurement result distribution of given input circuit and input state by running compute&compile function for desired number of shots

Return calculated result distribution and raw result count
"""
function get_distribution(input_circuit::Circuit, rt::R, input_state::Union{Stabilizer, Nothing}=nothing, num_shots::Int=1000) where R<:AbstractRuntime
    num_bits = get_bit_number(input_circuit)
    num_bits > 30 && throw(ArgumentError("Circuit uses $num_bits classical bits; the 2^$num_bits-entry distribution would not fit in memory"))
    len = 2^num_bits
    distribution = zeros(Int, len)
    data = zeros(Int, num_shots)
    # Compilation is shot-independent -- only the measurement outcomes differ --
    # so it is done once and each shot runs off a copy of the compiled state
    compiled = build_compilerstate(input_circuit, rt, input_state)
    for i in 1:num_shots
        result_i=copy(compiled)
        while !_execution_complete(result_i)
            result_i=execute!(result_i)
        end
        register=result_i.classical_register
        final_measurement_results=register[1:num_bits]
        any(isnothing, final_measurement_results) &&
            error("Classical bits 1:$num_bits were not all measured; register: $register")
        bit_str = join(Int.(final_measurement_results))
        # A circuit without measurements has an empty bit string; count it as outcome 0
        index = isempty(bit_str) ? 1 : parse(Int, bit_str; base=2) + 1
        distribution[index]+=1
        data[i]=index-1
    end
    return (distribution, data)
end

"""Valid values for the `qubits` keyword of [`get_graph`](@ref)/[`get_hypergraph`](@ref)."""
const QUBIT_VIEWS = (:all, :data, :magic)

"""Valid values for the `variant` keyword of [`get_graph`](@ref)/[`get_hypergraph`](@ref)."""
const RESULT_VARIANTS = (:all, :quantum, :determ, :random, :biased)

# `variant` value -> the MeasurementResult variant it keeps
const _VARIANT_NAME = Dict(
    :quantum => :QuantumRes,
    :determ  => :ClassicalDetermRes,
    :random  => :ClassicalRandomRes,
    :biased  => :ClassicalBiasedRes,
)

"""
`result.measurement_results` skipping unassigned slots.
"""
assigned_measurements(result::CompilationResult) =
    [result.measurement_results[i] for i in eachindex(result.measurement_results)
     if isassigned(result.measurement_results, i)]

"""
Keep only the measurements whose `MeasurementResult` variant matches `variant`
(`:all` keeps everything).
"""
filter_variant(measurements, variant::Symbol) =
    variant === :all ? measurements :
    filter(m -> variant_name(m) === _VARIANT_NAME[variant], measurements)

"""
    qubit_row_map(view, n_input, register_n) -> (row_of, nrows)

`row_of(q)` maps a register qubit index to its vertex row, or `nothing` when
the qubit is not included. `nrows` is the resulting number of vertices.

`:data` and `:magic` are complementary and symmetric: each keeps its own
block's qubits as individual vertices (re-indexed from 1) and drops the
other block entirely -- neither collapses multiple qubits onto one vertex.
(The playground's plotting-only `--qubits data` view does collapse the magic
block into a single lane for display; that collapsing is specific to
plotting and intentionally not reproduced here, since it would conflate the
magic block's resource footprint with whichever data qubit absorbed the
lane.)
"""
function qubit_row_map(view::Symbol, n_input::Int, register_n::Int)
    view === :all   && return (q -> (1 <= q <= register_n ? q : nothing), register_n)
    view === :magic && return (q -> (q > n_input ? q - n_input : nothing), register_n - n_input)
    view === :data  && return (q -> (q <= n_input ? q : nothing), n_input)
    error("unknown qubit view $view; must be one of $(join(QUBIT_VIEWS, ", "))")
end

"""
Vertex rows the measurement's Pauli support occupies.
"""
function measurement_rows(meas, row_of, register_n::Int)
    rows = Set{Int}()
    p = meas.pauli
    for q in 1:min(Int(nqubits(p)), register_n)
        p[q] == (false, false) && continue
        r = row_of(q)
        r === nothing || push!(rows, r)
    end
    return rows
end

"""
The support rule: keep only measurements with support in the selected rows.
"""
select_measurements(measurements, row_of, register_n::Int) =
    [m for m in measurements if !isempty(measurement_rows(m, row_of, register_n))]

"""
Select the measurements to build an interaction graph/hypergraph from, and the
vertex count to size it over.

`qubits` picks the vertex rows (`:all` every register qubit, `:data` the data
qubits alone, `:magic` the magic-state block alone -- both re-indexed from 1);
`:data`/`:magic` require `n_input` (the data-qubit count), since
`CompilationResult` does not record the data/magic split. A measurement
becomes an edge/hyperedge only if its Pauli has support in the selected rows.

The register width used to size and bound this is `max(size(stabilizer_group,
2), widest assigned measurement Pauli)`, not `size(stabilizer_group, 2)`
alone: [`AbstractStabilizerRuntime`](@ref)'s `to_result`
(`logic.jl`) slices `stabilizer_group` down to the data qubits before storing
it in `CompilationResult`, so its tableau width alone is the *data*-qubit
count, not the full register -- while `measurement_results` there keeps the
original full-register Paulis. For [`SimRuntime`](@ref)/[`DummyRuntime`](@ref)
results `stabilizer_group` is already full-register-width and this reduces to
that width. Either way this sizes off `measurement_results`/`stabilizer_group`
rather than trusting `QPU_workload`'s Pauli width, which is only
magic-register-wide for `SimRuntime`/`DummyRuntime` results and
full-register-wide for `AbstractStabilizerRuntime` results.

`variant` independently restricts which `MeasurementResult` variant is
eligible to become an edge/hyperedge (`:all`, `:quantum`, `:determ`,
`:random`). Note this measurement-support view is drawn from
`measurement_results`, not `QPU_workload`, so `qubits=:magic, variant=:quantum`
approximates but does not exactly reproduce the old `quantum_only=true`
selection -- see [`get_hypergraph`](@ref) for the precise difference.

Returns `(measurements, num_vertices, row_of, register_n)`: `row_of` must be
used to map each measurement's raw (full-register) qubit indices to vertex
rows before building edges/hyperedges, since `:data`/`:magic` renumber
qubits rather than using them directly.

Shared by `get_graph` and `get_hypergraph`.
"""
function _select_measurement_results(result::CompilationResult;
        qubits::Symbol=:all, n_input::Union{Int,Nothing}=nothing, variant::Symbol=:all)
    qubits in QUBIT_VIEWS ||
        error("qubits must be one of $(join(QUBIT_VIEWS, ", ")); got $qubits")
    variant in RESULT_VARIANTS ||
        error("variant must be one of $(join(RESULT_VARIANTS, ", ")); got $variant")
    (qubits === :data || qubits === :magic) && n_input === nothing &&
        throw(ArgumentError("qubits=$qubits requires n_input (the data-qubit count)"))
    all_m = assigned_measurements(result)
    tableau_n = size(result.stabilizer_group, 2)
    register_n = max(tableau_n,
        isempty(all_m) ? 0 : maximum(Int(nqubits(m.pauli)) for m in all_m))
    (row_of, num_vertices) = qubit_row_map(qubits, something(n_input, 0), register_n)
    measurements = select_measurements(filter_variant(all_m, variant), row_of, register_n)
    return (measurements, num_vertices, row_of, register_n)
end

"""
    get_graph(result::CompilationResult; qubits=:all, n_input=nothing, variant=:all) -> SimpleWeightedGraph

Extract the qubit interaction graph from a `CompilationResult`.

`qubits` selects which qubits become vertices: `:all` (default) every register
qubit, `:data` the data qubits alone, `:magic` the
magic-state block alone (re-indexed from 1). `:data`/`:magic` require
`n_input`, the number of data qubits, since `CompilationResult` does not
record the data/magic split. A measurement becomes an edge only if its Pauli
has support among the selected qubits. `variant` independently restricts
which `MeasurementResult` variant is eligible (`:all`, `:quantum`, `:determ`,
`:random`).

`qubits=:all` (the default) matches the old `quantum_only=false` default,
sized off the register width computed by [`_select_measurement_results`](@ref)
(not simply `size(result.stabilizer_group, 2)` -- see there for why).
`qubits=:magic, variant=:quantum` is the closest equivalent to the old
`quantum_only=true`,
but is not identical to it: unlike `QPU_workload`, this view is drawn from
`measurement_results` filtered by qubit support, so for
[`SimRuntime`](@ref)/[`DummyRuntime`](@ref) results it is a superset (weight-1
magic Paulis that `to_result` absorbed locally are not in `QPU_workload` but
do show up here), and for [`AbstractStabilizerRuntime`](@ref) results the two
are unrelated (`QPU_workload` there keeps full-register Paulis, including ones
with no magic support at all, which this view drops).
"""
function get_graph(result::CompilationResult;
        qubits::Symbol=:all, n_input::Union{Int,Nothing}=nothing, variant::Symbol=:all)
    (measurements, num_nodes, row_of, register_n) =
        _select_measurement_results(result; qubits, n_input, variant)
    g=SimpleWeightedGraph{Int64, Int64}(Int64(num_nodes))
    for m in measurements
        p=m.pauli
        rows = sort!(collect(measurement_rows(m, row_of, register_n)))
        for a in 1:length(rows)
            for b in a+1:length(rows)
                current_w = get_weight(g, rows[a], rows[b])
                add_edge!(g, rows[a], rows[b], current_w + 1)
            end
        end
    end
    return g
end

"""
Edge/hyperedge weight when merging `N` graphs/hypergraphs into one: `+1` if
the edge appears in every input, `-1` otherwise. Shared by `variant_graph` and
`variant_hypergraph`.
"""
_edge_occurrence_weight(count::Int, N::Int) = count == N ? 1 : -1

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

    edge_count = countmap([minmax(Int(src(e)), Int(dst(e))) for g in graphs for e in edges(g)])

    I_vals = Int[]
    J_vals = Int[]
    V_vals = Int[]
    for ((u, v), count) in edge_count
        w = _edge_occurrence_weight(count, N)
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
    weight_std_graph(input_circuit::Circuit, rt::AbstractRuntime, input_state=nothing;
                      qubits=:all, n_input=nothing, variant=:all, num_shots=1000)

Run the circuit `num_shots` times, extract the interaction graph of each shot with
[`get_graph`](@ref), and return a new `SimpleWeightedGraph` whose edge weights are the
standard deviation of each edge's weight across the shots.

`qubits`, `n_input`, and `variant` are forwarded to `get_graph`.

Edges absent from a shot's graph contribute weight 0 to the computation.
"""
function weight_std_graph(input_circuit::Circuit, rt::R, input_state::Union{Stabilizer, Nothing}=nothing;
        qubits::Symbol=:all, n_input::Union{Int,Nothing}=nothing, variant::Symbol=:all, num_shots::Int=1000) where R<:AbstractRuntime
    graphs =  Vector{SimpleWeightedGraph}(undef, num_shots)
    # See `get_distribution`: compile once, run each shot off a copy
    compiled = build_compilerstate(input_circuit, rt, input_state)
    for i in 1:num_shots
        state_i=copy(compiled)
        while !_execution_complete(state_i)
            state_i=execute!(state_i)
        end
        graphs[i]=get_graph(to_result(state_i); qubits, n_input, variant)
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
    get_hypergraph(result::CompilationResult; qubits=:all, n_input=nothing, variant=:all) -> (A, h)

Extract the qubit interaction hypergraph from a `CompilationResult`.

`qubits` selects which qubits become vertices: `:all` (default) every register
qubit, `:data` the data qubits alone, `:magic` the
magic-state block alone (re-indexed from 1). `:data`/`:magic` require
`n_input`, the number of data qubits, since `CompilationResult` does not
record the data/magic split. A measurement becomes a hyperedge only if its
Pauli has support among the selected qubits. `variant` independently
restricts which `MeasurementResult` variant is eligible (`:all`, `:quantum`,
`:determ`, `:random`).

`qubits=:all` (the default) matches the old `quantum_only=false` default,
sized off the register width computed by [`_select_measurement_results`](@ref)
(not simply `size(result.stabilizer_group, 2)` -- see there for why).
`qubits=:magic, variant=:quantum` is the closest equivalent to the old
`quantum_only=true`,
but is not identical to it: unlike `QPU_workload`, this view is drawn from
`measurement_results` filtered by qubit support, so for
[`SimRuntime`](@ref)/[`DummyRuntime`](@ref) results it is a superset (weight-1
magic Paulis that `to_result` absorbed locally are not in `QPU_workload` but
do show up here), and for [`AbstractStabilizerRuntime`](@ref) results the two
are unrelated (`QPU_workload` there keeps full-register Paulis, including ones
with no magic support at all, which this view drops).
"""
function get_hypergraph(result::CompilationResult;
        qubits::Symbol=:all, n_input::Union{Int,Nothing}=nothing, variant::Symbol=:all)
    # Fix the vertex count explicitly; inferring it from the sparse indices would
    # silently drop qubits that no measurement touches
    (measurements, num_v, row_of, register_n) =
        _select_measurement_results(result; qubits, n_input, variant)
    paulis=[m.pauli for m in measurements]
    collected_edges = Vector{Vector{Int}}()
    for p in paulis
        push!(collected_edges, _qubit_coverage(p, row_of, register_n))
    end
    w=countmap(collected_edges)
    I = Int[]
    J = Int[]
    i=1
    cleaned_w = filter(p-> length(p.first) > 1, w)
    for row in keys(cleaned_w)
        append!(I,row)
        col=fill(i,length(row))
        append!(J,col)
        i+=1
    end
    V = Int.(ones(length(I)))
    edge_weights=collect(values(cleaned_w))
    A = sparse(I, J, V, num_v, length(edge_weights))
    h = KaHyPar.HyperGraph(A,ones(Int, num_v),edge_weights)
    return (A, h)
end

"""
Vertex rows `p`'s support maps to under `row_of`, sorted and deduplicated.
"""
function _qubit_coverage(p::PauliOperator, row_of, register_n::Int)
    bool_vec = [p[i] for i in 1:min(Int(nqubits(p)), register_n)]
    idx = findall(x -> x !== (false,false), bool_vec)
    rows = Set{Int}()
    for i in idx
        r = row_of(i)
        r === nothing || push!(rows, r)
    end
    return sort!(collect(rows))
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

    edge_count = _hyperedge_counts(hypergraphs)

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
    edge_weights = [_edge_occurrence_weight(edge_count[e], N) for e in unique_edges]

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
    return count(verts -> !allequal(parts[v+1] for v in verts), _hyperedge_vertices(h))
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
    return Dict(verts .+ 1 => c for (verts, c) in _hyperedge_counts(hypergraphs))
end
##
"""
    _hyperedge_vertices(hg::KaHyPar.HyperGraph) -> Vector{Vector{Int}}

Decode the CSR edge structure of a KaHyPar hypergraph into one vector of
0-based vertex indices per hyperedge.
"""
function _hyperedge_vertices(hg::KaHyPar.HyperGraph)
    return [Int.(hg.hyperedges[hg.edge_indices[j]+1 : hg.edge_indices[j+1]])
            for j in 1:length(hg.edge_indices)-1]
end

"""
Count occurrences of each unique hyperedge across `hypergraphs`, keyed by the
sorted 0-based vertex vector.
"""
function _hyperedge_counts(hypergraphs)
    return countmap([sort!(verts) for hg in hypergraphs for verts in _hyperedge_vertices(hg)])
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
    return countmap(length.(_hyperedge_vertices(h)))
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
    return sum((length(unique(parts[v+1] for v in verts)) - 1
                for verts in _hyperedge_vertices(h)); init=0)
end
