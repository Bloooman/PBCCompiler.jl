    """
    This file contains the main iteration structures and functions for the PBCCompiler.
    It defines how to propagate Pauli Product Measurements (PPMs) through a quantum circuit
    Provides functions to handle commute condition and anticommute condition
    Provides functions to handle PBC List update and Check list update
    """


    using PBCCompiler: CircuitOp
    using QuantumClifford: comm, embed

    #I need a function to check commutation between the Pauli string of two CircuitOps
    function check_commutation(op1::CircuitOp.Type, op2::CircuitOp.Type)
        #Qubits affected by both CircuitOps
        AffectedQubits=sort(union(op1.qubits,op2.qubits))
        #Total number of qubits affected
        Paulilen=maximum(AffectedQubits)
        
        Pauli1=embed(Paulilen, op1.qubits, op1.pauli)
        Pauli2=embed(Paulilen, op2.qubits, op2.pauli)

        commutativity=comm(Pauli1,Pauli2)

        return commutativity   
    end

