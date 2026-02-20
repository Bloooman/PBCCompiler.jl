    """
    This file contains the main iteration structures and functions for the PBCCompiler.
    It defines how to propagate Pauli Product Measurements (PPMs) through a quantum circuit
    Provides functions to handle commute condition and anticommute condition
    Provides functions to handle PBC List update and Check list update
    """


using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, Measurement, ExpHalfPiPauli, ExpQuatPiPauli, ExpEighPiPauli, PauliConditional, BitConditional, affectedqubits
using QuantumClifford: comm, embed
using Moshi.Match: @match


function affectedpaulis(op::CircuitOp.Type)
    qubits = @match op begin
        CircuitOp.Pauli(pauli, qubits) => pauli
        CircuitOp.Measurement(pauli, bit, qubits) => pauli
        CircuitOp.ExpHalfPiPauli(pauli, qubits) => pauli
        CircuitOp.ExpQuatPiPauli(pauli, qubits) => pauli
        CircuitOp.ExpEighPiPauli(pauli, qubits) => pauli
        CircuitOp.PauliConditional(cp, cq, tp, tq) => vcat(cp, tp)
    end
end

function complete_paulis(op1::CircuitOp.Type, op2::CircuitOp.Type)
    pu1=affectedpaulis(op1)
    pu2=affectedpaulis(op2)
    println("Affected Paulis of op1: ", pu1)
    println("Affected Paulis of op2: ", pu2)
    qu1=affectedqubits(op1)
    qu2=affectedqubits(op2)
    println("Affected qubits of op1: ", qu1)
    println("Affected qubits of op2: ", qu2)
    AffectedQubbits=sort(union(qu1,qu2))
    println("Affected qubits of both ops: ", AffectedQubbits)
    Paulilen=maximum(AffectedQubbits)
    println("Length of the affected Pauli string: ", Paulilen)
    Pauli1=embed(Paulilen, qu1, pu1)
    Pauli2=embed(Paulilen, qu2, pu2)
    println("New Pauli string of op1: ", Pauli1)
    println("New Pauli string of op2: ", Pauli2)
    return (Pauli1, Pauli2)
end

function check_commutation(op1::CircuitOp.Type, op2::CircuitOp.Type)
    
    @match (op1, op2) begin
        #scenario 1: One of them is classical controlled gate
        (op,CircuitOp.BitConditional(inner_op, bit)) || (CircuitOp.BitConditional(inner_op, bit), op) => begin
            println("Invaid input: Need to determine gate present first")
        end
        #scenario 2: One of them is Pauli Conditional gate
        (op, CircuitOp.PauliConditional(cp, cq, tp, tq)) || (CircuitOp.PauliConditional(cp, cq, tp, tq), op) => begin
            println("One of the operations is a Pauli conditional gate. ")
            cop=ExpQuatPiPauli(cp, cq)
            top=ExpQuatPiPauli(tp, tq)
            comm_cop=check_commutation(op, cop)
            comm_top=check_commutation(op, top)
            if comm_cop == 0 && comm_top == 0
                println("The operation commutes with both the control and target Paulis of the conditional gate.")
            elseif comm_cop != 0 && comm_top != 0
                println("The operation anticommutes with both the control and target Paulis of the conditional gate.")
            else
                println("The operation commutes with one of the control and target Paulis of the conditional gate, and anticommutes with the other.")
            end
            return (comm_cop, comm_top)
        end
        #scenario 3: Both are Pauli Product Rotations
        _ => begin
            (Pauli1,Pauli2)=complete_paulis(op1, op2)
            commutativity=comm(Pauli1,Pauli2)
            if commutativity == 0 
                println("The two operations commute.")
            else
                println("The two operations anticommute.")
            end
            return commutativity
        end
        #scenario 2: one is a Controlled gate
      
    end
end

