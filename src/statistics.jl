    """
    This file contains functions for analyzing computed&compiled circuits
    """

##
using Graphs, SimpleWeightedGraphs
using QuantumClifford: xbit, zbit
using Statistics

"""
    get_distribution(input_circuit::Circuit, input_state::Stabilizer, num_shots=1000)

Get measurement result distribution of given input circuit and input state by running compute&compile function for desired number of shots

Return calculated result distribution and raw result count
"""
function get_distribution(input_circuit::Circuit, input_state::Stabilizer, num_shots::Int=1000)
    num_bits = _get_circuit_width(input_circuit)
    len = 2^num_bits
    distribution = zeros(Int, len)
    data = zeros(Int, num_shots)
    i=1
    while i<num_shots
        circuit = copy(input_circuit)
        result_i=run(circuit, input_state)
        register=result_i.memory_state.classical_register
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
    get_graph(result::ComputerState, type=nothing)

Extract qubit interaction graph from resulted ComputerState

When type is nothing(default), plot interaction graph among all qubits using all Pauli Product Measurement
When type is a MeasurementResultType, plot interaction graph using only MeasurementResult of given type
When type is QuantumRes, plot interaction graph of magic qubits denoted in ComputerState only
"""
function get_graph(result::ComputerState, type::Union{MeasurementResultType.Type,Nothing}=nothing)
    num_nodes=maximum(vcat(result.memory_state.pauli_qubits,result.memory_state.magic_qubits))
    g=SimpleWeightedGraph(Int64(num_nodes))
    measurements = type === nothing ?
        result.memory_state.measurement_results :
        filter(mr -> mr.result_type == type, result.memory_state.measurement_results)
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
    if type == QuantumRes()
        for h in 1:length(result.memory_state.pauli_qubits)
            rem_vertex!(g,1)
        end
    end
    return g
end
##
"""
    weight_std_graph(graphs::Vector{<:SimpleWeightedGraph})

Given a list of `SimpleWeightedGraph` objects with the same number of vertices,
return a new `SimpleWeightedGraph` whose edge weights are the standard deviation
of each edge's weight across the input graphs.

Edges absent from a graph in the list contribute weight 0 to the computation.
"""
function weight_std_graph(input_circuit::Circuit, input_state::Stabilizer; type::Union{MeasurementResultType.Type,Nothing}=nothing, num_shots::Int=1000)
    graphs =  Vector{SimpleWeightedGraph}(undef, num_shots)
    i=1
    while i<num_shots+1
        circuit = copy(input_circuit)
        result_i=run(circuit, input_state)
        graphs[i]=get_graph(result_i, type)
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
