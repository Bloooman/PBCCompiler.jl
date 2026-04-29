"""
Functions for compiling quantum circuits by moving measurement operations to the beginning of the circuit.
"""
struct PauliQubitMismatchError <: Exception
    msg::String
end

"""Function for checking if the pauli and qubits field denotes different number of qubits"""
function _validate_CircuitOp(op::CircuitOp.Type)
    p=_affectedpaulis(op)
    q=affectedqubits(op)
    name=variant_name(op)
    if length(p) != length(q)
        throw(PauliQubitMismatchError("$name($p, $q): The length of the Pauli string is not the same as the number of affected qubits. Please check the input operation."))
    end
    @match op begin
        CircuitOp.PauliConditional(cp, cq, tp, tq) => begin
            if cq == Int64[] || tq == Int64[]
                throw(PauliQubitMismatchError("$name($p, $q): Pauli String can't be empty"))
            end
        end
        _ => nothing
    end
end

"""Check every CircuitOp in a circuit"""
function _validate_circuit(circuit::Circuit)
    for op in circuit
        _validate_CircuitOp(op)
    end
end

function _get_circuit_width(circuit::Circuit)
    width=0
    for i in circuit
        width=max(width,maximum(affectedqubits(i)))
    end
    return width
end

function _get_bit_number(circuit::Circuit)
    num_bit=0
    for i in _find_measurement_indices(circuit)
        bits = @match circuit[i] begin
            CircuitOp.Measurement(pauli, bit, qubits) => bit
            _ => nothing
        end
        num_bit = max(num_bit,bits)
    end
    return num_bit
end

function _find_measurement_indices(circuit::Circuit)
    measurement_indices = []
    for (index, op) in enumerate(circuit)
        if isa_variant(op,CircuitOp.Measurement)
            push!(measurement_indices, index)
        else
            nothing
        end
    end
    return measurement_indices
end

function _find_nonclifford_indices(circuit::Circuit)
    nonclifford_indices = []
    for (index, op) in enumerate(circuit)
        if isa_variant(op, CircuitOp.ExpEighPiPauli)
            push!(nonclifford_indices, index)
        end
    end
    return nonclifford_indices
end

"""
Function that replace all non-Clifford circuit operations with BitConditional CircuitOps
Each BitConditional CircuitOp contains a gadget(a set of four consecutive CircuitOps) for pi/8 rotation implementation:
    Realize pi/8 rotation by consuming a |T ⟩ ancilla state
    perform a joint measurement P ⊗ Z between data and ancilla,
    then apply a conditional Clifford correction
"""
function _gadgetize(circuit::Circuit, index::Int, num_input_qubit::Int, num_magic_state::Int)
    op=circuit[index]
    num_bit=_get_bit_number(circuit)
            if isa_variant(op,CircuitOp.ExpEighPiPauli)
                P=_affectedpaulis(op)
                Q=affectedqubits(op)
                magic_state=[num_input_qubit+num_magic_state]
                Pauli=tensor(P,P"Z")
                Qubit=[Q;magic_state]
                magic_bit_1=num_bit+2*num_magic_state-1
                magic_bit_2=num_bit+2*num_magic_state
                Measurement_1=CircuitOp.Measurement(Pauli,magic_bit_1,Qubit)
                Measurement_2=CircuitOp.Measurement(P"X", magic_bit_2, magic_state)
                BitConditional_1=CircuitOp.BitConditional(CircuitOp.ExpQuatPiPauli(P,Q),magic_bit_1)
                BitConditional_2=CircuitOp.BitConditional(CircuitOp.ExpHalfPiPauli(P,Q),magic_bit_2)
                gadget=[Measurement_1, Measurement_2, BitConditional_1, BitConditional_2]
        splice!(circuit, index, gadget)
    else
        nothing
        end
end

"""
s is the stabilized part of input_state defined by user in the form of a stabilzier group
Function will expand the stabilizer group to cover the entire circuit width by adding Identities to each stabilizer
"""
function _make_stabilizer_list(s::Stabilizer, circuit::Circuit)
    paulilen=_get_circuit_width(circuit)
    num_pauli_qubits = length(s)
    new_s=PauliOperator[]
    pauli_qubits = collect(1:num_pauli_qubits)
    for i in s
        new_i = embed(paulilen, pauli_qubits, i)
        push!(new_s,new_i)
    end
    return Stabilizer(new_s)
end
