@testitem "affectedqubits" tags=[:affectedqubits] begin

using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, Measurement, ExpHalfPiPauli, ExpQuatPiPauli, ExpEighPiPauli, PauliConditional, BitConditional, affectedqubits
using QuantumClifford: @P_str

@testset "Single operation qubits" begin
    # Pauli gate
    op = ExpQuatPiPauli(P"X", [1])
    @test affectedqubits(op) == [1]

    op = ExpQuatPiPauli(P"XY", [1, 3])
    @test affectedqubits(op) == [1, 3]

    # Measurement
    op = Measurement(P"Z", 0, [2])
    @test affectedqubits(op) == [2]

    op = Measurement(P"ZZ", 1, [1, 4])
    @test affectedqubits(op) == [1, 4]

    # ExpHalfPiPauli
    op = ExpHalfPiPauli(P"Y", [3])
    @test affectedqubits(op) == [3]

    # ExpQuatPiPauli
    op = ExpQuatPiPauli(P"X", [2])
    @test affectedqubits(op) == [2]

    # ExpEighPiPauli
    op = ExpEighPiPauli(P"Z", [5])
    @test affectedqubits(op) == [5]
end


@testset "PauliConditional qubits" begin
    # PauliConditional has control and target qubits
    op = PauliConditional(P"X", [1], P"Z", [2])
    @test affectedqubits(op) == [1, 2]

    op = PauliConditional(P"XX", [1, 2], P"ZZ", [3, 4])
    @test affectedqubits(op) == [1, 2, 3, 4]

    # Overlapping qubits should be deduplicated
    op = PauliConditional(P"X", [1], P"Z", [1])
    @test affectedqubits(op) == [1]
end

@testset "BitConditional qubits" begin
    # BitConditional wraps another operation
    inner = ExpQuatPiPauli(P"XY", [2, 3])
    op = BitConditional(inner, 0)
    @test affectedqubits(op) == [2, 3]

    inner = Measurement(P"Z", 1, [5])
    op = BitConditional(inner, 0)
    @test affectedqubits(op) == [5]
end

@testset "Circuit qubits" begin
    # Empty circuit
    circuit = Circuit()
    @test affectedqubits(circuit) == Int[]

    # Single operation circuit
    circuit = Circuit([ExpQuatPiPauli(P"X", [1])])
    @test affectedqubits(circuit) == [1]

    # Multiple operations
    circuit = Circuit([
        ExpQuatPiPauli(P"X", [1]),
        Measurement(P"Z", 0, [2]),
        ExpHalfPiPauli(P"Y", [3])
    ])
    @test affectedqubits(circuit) == [1, 2, 3]

    # Operations with overlapping qubits
    circuit = Circuit([
        ExpQuatPiPauli(P"XX", [1, 2]),
        ExpQuatPiPauli(P"YY", [2, 3])
    ])
    @test affectedqubits(circuit) == [1, 2, 3]

    # Non-contiguous qubits
    circuit = Circuit([
        ExpQuatPiPauli(P"X", [1]),
        ExpQuatPiPauli(P"Z", [5]),
        ExpQuatPiPauli(P"Y", [3])
    ])
    @test affectedqubits(circuit) == [1, 3, 5]
end

end
