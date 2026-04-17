function get_distribution(input_circuit::Circuit, input_state::Stabilizer, num_shots::Int)
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
