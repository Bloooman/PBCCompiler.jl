    """
    This file contains the main iteration structures and functions for the PBCCompiler.
    It defines how to propagate Pauli Product Measurements (PPMs) through a quantum circuit
    Provides functions to handle commute condition and anticommute condition
    Provides functions to handle PBC List update and Check list update
    """


using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, Measurement, ExpHalfPiPauli, ExpQuatPiPauli, ExpEighPiPauli, PauliConditional, BitConditional, affectedqubits
using QuantumClifford: comm, embed, ⊗
using Moshi.Match: @match
using Moshi.Data: variant_name


function affectedpaulis(op::CircuitOp.Type)
    pauli = @match op begin
        CircuitOp.Measurement(pauli, bit, qubits) => pauli
        CircuitOp.ExpHalfPiPauli(pauli, qubits) => pauli
        CircuitOp.ExpQuatPiPauli(pauli, qubits) => pauli
        CircuitOp.ExpEighPiPauli(pauli, qubits) => pauli
        CircuitOp.PauliConditional(cp, cq, tp, tq) => vcat(cp, tp)
    end
    return pauli
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
    Pauli1=embed(Paulilen, op1.qubits, pu1)
    Pauli2=embed(Paulilen, op2.qubits, pu2)
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
        #scenario 3: One of them is HalfPi Pauli
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

function conjugate(op1::CircuitOp.Type, op2::CircuitOp.Type) #first input is the one we conjugate by, second input is the one we want to conjugate
    conjugated_op=@match (op1, op2) begin
     #scenario 1: one is a BitControlled gate
        (op,CircuitOp.BitConditional(inner_op, bit)) || (CircuitOp.BitConditional(inner_op, bit), op) => begin
            println("Invaid input: Need to determine gate present first")
        end
    #scenario 2: Conjugated by a PauliControlled gate
        (CircuitOp.PauliConditional(cp, cq, tp, tq), op) => begin
            println("One of the operations is a Pauli conditional gate.")
            op_1=ExpQuatPiPauli(-cp, cq)
            println("First conjugation with the control Pauli of the conditional gate.")
            op_2=ExpQuatPiPauli(-tp, tq)
            println("Second conjugation with the target Pauli of the conditional gate.")
            op_3=ExpQuatPiPauli(cp⊗tp, sort(union(cq, tq)))
            println("Third conjugation with the combined control and target Paulis of the conditional gate.")
            conju_step1=conjugate(op_1, op2)[1]
            conju_step2=conjugate(op_2, conju_step1)[1]
            conju_final=conjugate(op_3, conju_step2)[1]
            conju_final
        end
    #scenario 3: Conjugated by a HalfPi Pauli
        (CircuitOp.ExpHalfPiPauli(p1,q1),op) => begin
        println("One of the operations is a HalfPi Pauli gate.")
            if check_commutation(op1,op2) == 0
                new_p=complete_paulis(op1, op2)[2]
                new_qm=maximum(sort(union(q1, affectedqubits(op2))))
                new_q=[x for x in 1:new_qm]
                println("The two operations commute, no change after conjugation.")
                println("The Pauli string of the conjugated operation is: ", new_p)
                println("The qubits affected by the conjugated operation are: ", new_q)
            else
                (pauli1,pauli2)=complete_paulis(op1, op2)
                new_p=-pauli2
                new_qm=maximum(sort(union(affectedqubits(op1), affectedqubits(op2))))
                new_q=[x for x in 1:new_qm]
                println("The two operations anticommute, the Pauli string of the conjugated operation will be changed after conjugation.")
                println("The Pauli string of the conjugated operation is: ", new_p)
                println("The qubits affected by the conjugated operation are: ", new_q)
            end
            typeofp=variant_name(op2)
            if typeofp == :Measurement
                constructor=getproperty(CircuitOp, typeofp)
                new_op=constructor(new_p, op2.bit, new_q)
            else
            constructor=getproperty(CircuitOp, typeofp)
            new_op=constructor(new_p, new_q)
            end

        end
    #scenario 4: PPM Conjugated by a ExpQuatPi Pauli
        (CircuitOp.ExpQuatPiPauli(p1,q1), CircuitOp.Measurement(p2,b,q2)) => begin
            if check_commutation(op1,op2) == 0
                new_p=complete_paulis(op1, op2)[2]
                new_qm=maximum(sort(union(q1, q2)))
                new_q=[x for x in 1:new_qm]
                println("The two operations commute, no change after conjugation.")
                println("The Pauli string of the conjugated operation is: ", new_p)
                println("The qubits affected by the conjugated operation are: ", new_q)
            else
                (pauli1,pauli2)=complete_paulis(op1, op2)
                new_p=1im*pauli1*pauli2
                new_qm=maximum(sort(union(q1, q2)))
                new_q=[x for x in 1:new_qm]
                println("The two operations anticommute, the Pauli string of the conjugated operation will be changed after conjugation.")
                println("The Pauli string of the conjugated operation is: ", new_p)
                println("The qubits affected by the conjugated operation are: ", new_q)
            end
            Measurement(new_p,b,new_q)
        end
    #scenario 5: PPR Conjugated by a ExpQuatPi Pauli
        (CircuitOp.ExpQuatPiPauli(p1,q1), op) => begin
            if check_commutation(op1,op2) == 0
                new_p=complete_paulis(op1, op2)[2]
                new_qm=maximum(sort(union(q1, affectedqubits(op2))))
                new_q=[x for x in 1:new_qm]
                println("The two operations commute, no change after conjugation.")
                println("The Pauli string of the conjugated operation is: ", new_p)
                println("The qubits affected by the conjugated operation are: ", new_q)
            else
                (pauli1,pauli2)=complete_paulis(op1, op2)
                new_p=1im*pauli1*pauli2
                new_qm=maximum(sort(union(affectedqubits(op1), affectedqubits(op2))))
                new_q=[x for x in 1:new_qm]
                println("The two operations anticommute, the Pauli string of the conjugated operation will be changed after conjugation.")
                println("The Pauli string of the conjugated operation is: ", new_p)
                println("The qubits affected by the conjugated operation are: ", new_q)
            end
            typeofp=variant_name(op2)
            constructor=getproperty(CircuitOp, typeofp)
            new_op=constructor(new_p, new_q)

        end
        _=> begin
            println("Invalid input: Conjugation between the given types of CircuitOps is not supported.")
        end

    end
    return (conjugated_op,op1)
end
