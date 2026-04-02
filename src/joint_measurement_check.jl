"""
Helper functions to check the first PPM in circuit, determine MeasurementResultType: ClassicalDetermRes, ClassicalRandomRes, QuantumRes
"""

function check_PPM(s::Stabilizer,op::CircuitOp.Type, num_qubits::Int)
    if !isa_variant(op,CircuitOp.Measurement)
        return nothing
    else
        Paulilen = num_qubits
        Pauli=embed(Paulilen, op.qubits, op.pauli)
        return project!(copy(s),Pauli)
    end
end

"""
0x00 denotes +1 eigenvalue, 0x02 denotes -1 eigenvalue
0 denotes +1 eigenvalue, 1 denotes -1 eigenvalue
false denotes +1 eigenvalue, true denotes -1 eigenvalue
"""

function get_measurement_result(s::Stabilizer, op::CircuitOp.Type, num_qubits::Int)
    len=length(s)
    projection = check_PPM(s, op, num_qubits)
    if projection === nothing
        return nothing
    else
        if projection[3] === nothing
            if projection[2]<=len
                result = rand(Bool[0,1])
                return (classical_random_result(result),projection[2])
            else
                result = Bool(quantum_measurement(op)>>1)
                return quantum_result(result)
            end
        else
            result = Bool(projection[3]>>1)
            return classical_deterministic_result(result)
        end
    end
end

function quantum_measurement(op::CircuitOp.Type)
    print("Enter quantum measurement result: ")
    measurement_result = parse(Int,readline())
    if abs(measurement_result) != 1
        throw(ArgumentError("Measurement Result can only be +1 or -1!!!"))
    else
        if measurement_result == 1
            return 0x00
        elseif measurement_result == -1
            return 0x02
        else
            return nothing
        end
    end
end
