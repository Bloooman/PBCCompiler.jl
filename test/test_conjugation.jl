using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, Measurement, ExpHalfPiPauli, ExpQuatPiPauli, ExpEighPiPauli, PauliConditional, BitConditional, affectedqubits
using Test
using PBCCompiler: conjugate
using QuantumClifford: @P_str
using Moshi.Derive: @derive

@derive CircuitOp[Eq, Show]

@testset "Test unordered input" begin
    op1=ExpQuatPiPauli(P"XY", [1, 3])
    op2=ExpQuatPiPauli(P"ZXY", [3, 1, 2])
    conjugated_op=ExpQuatPiPauli(P"-_YX", [1, 2, 3])

    t_1=conjugate(op1,op2)
    @test t_1 == (conjugated_op, op1)

end

@testset "Test non-overlapping input" begin
    op1=ExpQuatPiPauli(P"X", [1])
    op2=ExpQuatPiPauli(P"Z", [2])
    conjugated_op=ExpQuatPiPauli(P"_Z", [1, 2])

    t_1=conjugate(op1,op2)
    @test t_1 == (conjugated_op, op1)

end
