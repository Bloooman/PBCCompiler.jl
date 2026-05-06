using Graphs, SimpleWeightedGraphs
using Metis

"""
    METIS_partition(g::SimpleWeightedGraph, num_partition::Int) -> Vector{Int32}

Partition the nodes of `g` into `num_partition` balanced parts using METIS k-way partitioning.

Edge weights are passed to METIS and must be positive integers.

# Arguments
- `g`: weighted undirected graph to partition
- `num_partition`: number of parts

# Returns
A `Vector{Int32}` of length `nv(g)` where element `i` is the 1-based partition index
(in `1:num_partition`) assigned to node `i`.
"""
function METIS_partition(g::SimpleWeightedGraph, num_partition::Int)::Vector{Int32}
    return Metis.partition(g, num_partition)
end
