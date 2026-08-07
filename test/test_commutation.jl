@testitem "check_commutation" tags=[:check_commutation] begin

using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, Measurement, ExpHalfPiPauli, ExpQuatPiPauli, ExpEighPiPauli, PauliConditional, BitConditional, affectedqubits
using Test
using PBCCompiler: check_commutation
using QuantumClifford: @P_str, comm
@testset "Test unordered input" begin
    op1=ExpQuatPiPauli(P"XY", [1, 3])
    op2=ExpQuatPiPauli(P"ZXY", [3, 1, 2])

    t_1=comm(op1.pauli,op2.pauli)
    @test t_1 == 0x00

    t_2=check_commutation(op1, op2)
    @test t_2 == 0x01

end

@testset "Test non-overlapping input" begin
    op1=ExpQuatPiPauli(P"X", [1])
    op2=ExpQuatPiPauli(P"Z", [2])

    t_1=comm(op1.pauli,op2.pauli)
    @test t_1 == 0x01

    t_2=check_commutation(op1, op2)
    @test t_2 == 0x00

end

@testset "Test partially overlapping input" begin
    op1=ExpQuatPiPauli(P"XX", [1, 3])
    op2=ExpQuatPiPauli(P"ZY", [3, 4])

    t_1=comm(op1.pauli,op2.pauli)
    @test t_1 == 0x00

    t_2=check_commutation(op1, op2)
    @test t_2 == 0x01

end

@testset "Test Pauli Product Rotation/Measurement input" begin
    op1=ExpHalfPiPauli(P"X", [1])
    op2=ExpQuatPiPauli(P"Z", [1])
    op3=ExpEighPiPauli(P"Y", [1])
    M_Z=CircuitOp.Measurement(P"Z", 1, [1])

    t_1=check_commutation(op1,op2)
    @test t_1 == 0x01

    t_2=check_commutation(op1, op2)
    @test t_2 == 0x01

    t_3=check_commutation(op1,op3)
    @test t_3 == 0x01

    t_4=check_commutation(op1,M_Z)
    @test t_4 == 0x01

    t_5=check_commutation(op2,M_Z)
    @test t_5 == 0x00

    t_6=check_commutation(op3,M_Z)
    @test t_6 == 0x01
end

@testset "Test invalid input" begin
    CNOT=PauliConditional(P"Z", [1], P"X", [2])
    Con_Z=BitConditional(ExpHalfPiPauli(P"Z", [1]), 1)
    nonclifford_op=CircuitOp.ExpEighPiPauli(P"Z", [2])
    M_Z=CircuitOp.Measurement(P"Z", 1, [2])

    t_1=check_commutation(CNOT,nonclifford_op)
    @test t_1 === nothing

    t_2=check_commutation(CNOT,M_Z)
    @test t_2 === nothing

    t_3=check_commutation(Con_Z,nonclifford_op)
    @test t_3 === nothing

    t_4=check_commutation(Con_Z,M_Z)
    @test t_4 === nothing
end

end

@testitem "paulis_commute matches check_commutation" tags=[:check_commutation] begin
##
using PBCCompiler
using PBCCompiler: CircuitOp, check_commutation, paulis_commute
using QuantumClifford: PauliOperator, @P_str

# All Paulis of length n as (x, z) bit patterns.
allpaulis(n) = [PauliOperator(0x00,
                              [b & (1 << (k - 1)) != 0 for k in 1:n],
                              [b & (1 << (n + k - 1)) != 0 for k in 1:n])
                for b in 0:(4^n - 1)]

mk(kind, p, q) = kind == 1 ? CircuitOp.ExpHalfPiPauli(p, q) :
                 kind == 2 ? CircuitOp.ExpQuatPiPauli(p, q) :
                 kind == 3 ? CircuitOp.ExpEighPiPauli(p, q) :
                             CircuitOp.Measurement(p, 1, q)

# Supports that overlap fully, partially and not at all, sorted and unsorted --
# unsorted lists are the case where pairing letters with qubits by position
# matters.
supports = [[1], [2], [1, 2], [2, 1], [2, 3], [1, 2, 3], [3, 2, 1], [1, 3], [4]]

@testset "agrees with the embedding-based check" begin
    mismatches = 0
    for qa in supports, qb in supports
        for pa in allpaulis(length(qa)), pb in allpaulis(length(qb))
            (iszero(pa.xz) || iszero(pb.xz)) && continue
            for ka in 1:4, kb in 1:4
                op1 = mk(ka, pa, qa)
                op2 = mk(kb, pb, qb)
                mismatches += (check_commutation(op1, op2) == 0x00) != paulis_commute(op1, op2)
            end
        end
    end
    @test mismatches == 0
end

@testset "does not allocate" begin
    op1 = CircuitOp.ExpQuatPiPauli(P"XY", [1, 3])
    op2 = CircuitOp.Measurement(P"ZXY", 1, [3, 1, 2])
    paulis_commute(op1, op2)  # compile
    @test (@allocated paulis_commute(op1, op2)) == 0
end
##
end
