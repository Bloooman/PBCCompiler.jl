using Graphs, SimpleWeightedGraphs
using QuantumClifford: xbit, zbit

function get_distribution(input_circuit::Circuit, input_state::Stabilizer, num_shots::Int=1000)
    num_bits = get_circuit_width(input_circuit)
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
        for v in result.memory_state.pauli_qubits
            rem_vertex!(g,v)
        end
    end
    return g
end
