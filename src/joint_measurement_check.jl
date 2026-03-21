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


function get_measurement_result(s::Stabilizer, op::CircuitOp.Type, num_qubits::Int)
    len=length(s)
    result = check_PPM(s, op, num_qubits)
    if result[3] === nothing
        if result[2]>len
            quantum_measurement(op)
        else
            return rand([-1,1])
        end
    else
        return result[3]
    end
end

function quantum_measurement(op::CircuitOp.Type)
    print("Enter quantum measurement result: ")
    measurement_result = parse(Int,readline())
    if abs(measurement_result) != 1
        throw(ArgumentError("Measurement Result can only be +1 or -1!!!"))
    else
        return measurement_result
    end
end
