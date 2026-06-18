@testitem "io" tags=[:io] begin

# test_parse.jl — unit tests for parse_input
using JLD2
using PBCCompiler: Circuit, CircuitOp, Measurement, ExpHalfPiPauli, ExpQuatPiPauli, ExpEighPiPauli, PauliConditional
using PBCCompiler: MeasurementResult, CompilerState, SimRuntime, CompilationResult
using .MeasurementResult: ClassicalDetermRes, ClassicalRandomRes, QuantumRes
using Moshi.Data: isa_variant
using QuantumClifford: @P_str, Stabilizer
using PBCCompiler: parse_input, save, load

"""
    with_qasm(f, gates)

Write a minimal OpenQASM 2.0 file containing `gates` to a temp path,
call `f(path)`, then delete the file. The header (`OPENQASM 2.0`,
`include`, `qreg q[3]`, `creg c[3]`) is prepended automatically.
"""
function with_qasm(f::Function, gates::String)
    path = tempname() * ".qasm"
    write(path, """OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
$gates
""")
    try
        f(path)
    finally
        rm(path; force=true)
    end
end

const FIXTURES = joinpath(@__DIR__, "fixtures")

@testset "Header lines produce no ops" begin
    with_qasm("") do path
        @test isempty(parse_input(path))
    end
end

@testset "Hadamard gate (h)" begin
    with_qasm("h q[0];") do path
        c = parse_input(path)
        # H → three consecutive ExpQuatPiPauli: Z, X, Z
        @test length(c) == 3
        @test all(isa_variant(op, CircuitOp.ExpQuatPiPauli) for op in c)
        @test c[1].pauli == P"Z" && c[1].qubits == [1]
        @test c[2].pauli == P"X" && c[2].qubits == [1]
        @test c[3].pauli == P"Z" && c[3].qubits == [1]
    end
end

@testset "Phase gate (s)" begin
    with_qasm("s q[1];") do path
        c = parse_input(path)
        @test length(c) == 1
        @test isa_variant(c[1], CircuitOp.ExpQuatPiPauli)
        @test c[1].pauli == P"Z" && c[1].qubits == [2]
    end
end

@testset "Phase-dagger gate (sdg)" begin
    with_qasm("sdg q[0];") do path
        c = parse_input(path)
        @test length(c) == 1
        @test isa_variant(c[1], CircuitOp.ExpQuatPiPauli)
        @test c[1].pauli == -P"Z" && c[1].qubits == [1]
    end
end

@testset "T gate (t)" begin
    with_qasm("t q[2];") do path
        c = parse_input(path)
        @test length(c) == 1
        @test isa_variant(c[1], CircuitOp.ExpEighPiPauli)
        @test c[1].pauli == P"Z" && c[1].qubits == [3]
    end
end

@testset "T-dagger gate (tdg)" begin
    with_qasm("tdg q[1];") do path
        c = parse_input(path)
        @test length(c) == 1
        @test isa_variant(c[1], CircuitOp.ExpEighPiPauli)
        @test c[1].pauli == -P"Z" && c[1].qubits == [2]
    end
end

@testset "Pauli X gate (x)" begin
    with_qasm("x q[0];") do path
        c = parse_input(path)
        @test length(c) == 1
        @test isa_variant(c[1], CircuitOp.ExpHalfPiPauli)
        @test c[1].pauli == P"X" && c[1].qubits == [1]
    end
end

@testset "Pauli Y gate (y)" begin
    with_qasm("y q[0];") do path
        c = parse_input(path)
        @test length(c) == 1
        @test isa_variant(c[1], CircuitOp.ExpHalfPiPauli)
        @test c[1].pauli == P"Y" && c[1].qubits == [1]
    end
end

@testset "Pauli Z gate (z)" begin
    with_qasm("z q[0];") do path
        c = parse_input(path)
        @test length(c) == 1
        @test isa_variant(c[1], CircuitOp.ExpHalfPiPauli)
        @test c[1].pauli == P"Z" && c[1].qubits == [1]
    end
end

@testset "CNOT gate (cx)" begin
    with_qasm("cx q[0],q[1];") do path
        c = parse_input(path)
        @test length(c) == 1
        @test isa_variant(c[1], CircuitOp.PauliConditional)
        @test c[1].control_pauli  == P"Z"  && c[1].control_qubits == [1]
        @test c[1].target_pauli   == P"X"  && c[1].target_qubits  == [2]
    end
end

@testset "Measurement" begin
    with_qasm("measure q[2] -> c[1];") do path
        c = parse_input(path)
        @test length(c) == 1
        @test isa_variant(c[1], CircuitOp.Measurement)
        @test c[1].pauli  == P"Z"
        @test c[1].qubits == [3]
        @test c[1].bit    == 2
    end
end

@testset "ccx (Toffoli) gate" begin
    c = parse_input(joinpath(FIXTURES, "toffoli3.qasm"))
    # h×2 (3 ops each) + 4 cx + t×4 + tdg×3 = 19 gate ops
    @test length(c) == 19

    # First 3 ops: h q[2] → ExpQuatPiPauli Z,X,Z on qubit 3
    @test isa_variant(c[1], CircuitOp.ExpQuatPiPauli) && c[1].pauli == P"Z" && c[1].qubits == [3]
    @test isa_variant(c[2], CircuitOp.ExpQuatPiPauli) && c[2].pauli == P"X" && c[2].qubits == [3]
    @test isa_variant(c[3], CircuitOp.ExpQuatPiPauli) && c[3].pauli == P"Z" && c[3].qubits == [3]

    # Op 4: cx q[1],q[2] → ctrl=2, tgt=3
    @test isa_variant(c[4], CircuitOp.PauliConditional)
    @test c[4].control_qubits == [2] && c[4].target_qubits == [3]

    # Op 19: last cx q[0],q[1] → ctrl=1, tgt=2
    @test isa_variant(c[19], CircuitOp.PauliConditional)
    @test c[19].control_qubits == [1] && c[19].target_qubits == [2]
end

@testset "barrier is skipped" begin
    with_qasm("barrier q;\nx q[0];") do path
        c = parse_input(path)
        @test length(c) == 1
        @test isa_variant(c[1], CircuitOp.ExpHalfPiPauli)
    end
end

@testset "Gate sequence ordering is preserved" begin
    with_qasm("s q[0];\nt q[0];\nsdg q[0];") do path
        c = parse_input(path)
        @test length(c) == 3
        @test isa_variant(c[1], CircuitOp.ExpQuatPiPauli) && c[1].pauli == P"Z"   # s
        @test isa_variant(c[2], CircuitOp.ExpEighPiPauli) && c[2].pauli == P"Z"   # t
        @test isa_variant(c[3], CircuitOp.ExpQuatPiPauli) && c[3].pauli == -P"Z"  # sdg
    end
end

@testset "Mixed multi-qubit sequence" begin
    with_qasm("h q[0];\ncx q[0],q[1];\nmeasure q[0] -> c[0];") do path
        c = parse_input(path)
        # h → 3 ops, cx → 1 op, measure → 1 op
        @test length(c) == 5
        @test all(isa_variant(c[i], CircuitOp.ExpQuatPiPauli) for i in 1:3)
        @test isa_variant(c[4], CircuitOp.PauliConditional)
        @test isa_variant(c[5], CircuitOp.Measurement)
    end
end

end
