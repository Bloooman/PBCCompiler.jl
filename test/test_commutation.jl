using PBCCompiler: MainIteration
using Test

@testset "Test check_commutation function" begin
    a=P"XY"
    b=P"ZXY"

    op1=CircuitOp.Pauli(a,[1,3])
    op2=CircuitOp.Pauli(b,[3,1,2])

    t_1=comm(op1.pauli,op2.pauli)
    @test t_1 == 0x00

    t_2=MainIteration.check_commutation(op1,op2)
    @test t_2 == 0x01

end