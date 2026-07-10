@testitem "io" tags=[:io] begin

# test_parse.jl — unit tests for parse_input
using JLD2
using PBCCompiler: Circuit, CircuitOp, Measurement, ExpHalfPiPauli, ExpQuatPiPauli, ExpEighPiPauli, PauliConditional
using PBCCompiler: MeasurementResult, CompilerState, SimRuntime, CompilationResult
using .MeasurementResult: ClassicalDetermRes, ClassicalRandomRes, QuantumRes
using Moshi.Data: isa_variant
using QuantumClifford: @P_str, Stabilizer
using PBCCompiler: parse_input, save_result, load_result
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

@testset "Header lines produce no gate ops" begin
    with_qasm("") do path
        c = parse_input(path)
        # qreg q[3] triggers measure-all; no gate ops from header lines
        @test length(c) == 3
        @test all(isa_variant(op, CircuitOp.Measurement) for op in c)
    end
end

@testset "Hadamard gate (h)" begin
    with_qasm("h q[0];") do path
        c = parse_input(path)
        # H → three consecutive ExpQuatPiPauli: Z, X, Z; then 3 measure-all
        @test length(c) == 6
        @test all(isa_variant(c[i], CircuitOp.ExpQuatPiPauli) for i in 1:3)
        @test c[1].pauli == P"Z" && c[1].qubits == [1]
        @test c[2].pauli == P"X" && c[2].qubits == [1]
        @test c[3].pauli == P"Z" && c[3].qubits == [1]
    end
end

@testset "Phase gate (s)" begin
    with_qasm("s q[1];") do path
        c = parse_input(path)
        @test length(c) == 4  # 1 gate + 3 measure-all
        @test isa_variant(c[1], CircuitOp.ExpQuatPiPauli)
        @test c[1].pauli == P"Z" && c[1].qubits == [2]
    end
end

@testset "Phase-dagger gate (sdg)" begin
    with_qasm("sdg q[0];") do path
        c = parse_input(path)
        @test length(c) == 4
        @test isa_variant(c[1], CircuitOp.ExpQuatPiPauli)
        @test c[1].pauli == -P"Z" && c[1].qubits == [1]
    end
end

@testset "T gate (t)" begin
    with_qasm("t q[2];") do path
        c = parse_input(path)
        @test length(c) == 4
        @test isa_variant(c[1], CircuitOp.ExpEighPiPauli)
        @test c[1].pauli == P"Z" && c[1].qubits == [3]
    end
end

@testset "T-dagger gate (tdg)" begin
    with_qasm("tdg q[1];") do path
        c = parse_input(path)
        @test length(c) == 4
        @test isa_variant(c[1], CircuitOp.ExpEighPiPauli)
        @test c[1].pauli == -P"Z" && c[1].qubits == [2]
    end
end

@testset "Pauli X gate (x)" begin
    with_qasm("x q[0];") do path
        c = parse_input(path)
        @test length(c) == 4
        @test isa_variant(c[1], CircuitOp.ExpHalfPiPauli)
        @test c[1].pauli == P"X" && c[1].qubits == [1]
    end
end

@testset "Pauli Y gate (y)" begin
    with_qasm("y q[0];") do path
        c = parse_input(path)
        @test length(c) == 4
        @test isa_variant(c[1], CircuitOp.ExpHalfPiPauli)
        @test c[1].pauli == P"Y" && c[1].qubits == [1]
    end
end

@testset "Pauli Z gate (z)" begin
    with_qasm("z q[0];") do path
        c = parse_input(path)
        @test length(c) == 4
        @test isa_variant(c[1], CircuitOp.ExpHalfPiPauli)
        @test c[1].pauli == P"Z" && c[1].qubits == [1]
    end
end

@testset "CNOT gate (cx)" begin
    with_qasm("cx q[0],q[1];") do path
        c = parse_input(path)
        @test length(c) == 4  # 1 gate + 3 measure-all
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

@testset "multiple quantum and classical registers" begin
    path = tempname() * ".qasm"
    write(path, """OPENQASM 2.0;
    qreg a[2];
    qreg b[1];
    creg c[2];
    creg d[1];
    x a[1];
    cx a[0],b[0];
    measure a[0] -> c[0];
    measure b[0] -> d[0];
    """)
    try
        c = parse_input(path)
        # a[0]→1, a[1]→2, b[0]→3; c[0]→1, c[1]→2, d[0]→3
        @test length(c) == 4
        @test isa_variant(c[1], CircuitOp.ExpHalfPiPauli) && c[1].qubits == [2]
        @test isa_variant(c[2], CircuitOp.PauliConditional)
        @test c[2].control_qubits == [1] && c[2].target_qubits == [3]
        @test isa_variant(c[3], CircuitOp.Measurement) && c[3].qubits == [1] && c[3].bit == 1
        @test isa_variant(c[4], CircuitOp.Measurement) && c[4].qubits == [3] && c[4].bit == 3
    finally
        rm(path; force=true)
    end
end

@testset "OpenQASM 3.0 registers and measure" begin
    path = tempname() * ".qasm"
    write(path, """OPENQASM 3.0;
    qubit[2] q;
    bit[2] m;
    x q[1];
    m[0] = measure q[1];
    """)
    try
        c = parse_input(path)
        @test length(c) == 2
        @test isa_variant(c[1], CircuitOp.ExpHalfPiPauli) && c[1].qubits == [2]
        @test isa_variant(c[2], CircuitOp.Measurement) && c[2].qubits == [2] && c[2].bit == 1
    finally
        rm(path; force=true)
    end
end

@testset "undeclared register reference errors" begin
    with_qasm("x nosuch[0];") do path
        @test_throws ErrorException parse_input(path)
    end
end

@testset "ccx (Toffoli) gate" begin
    c = parse_input(joinpath(FIXTURES, "toffoli3.qasm"))
    # h×2 (3 ops each) + 4 cx + t×4 + tdg×3 = 19 gate ops + 3 measure-all = 22
    @test length(c) == 22

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

    # Ops 20–22: measure-all on qubits 1,2,3
    @test isa_variant(c[20], CircuitOp.Measurement) && c[20].bit == 1 && c[20].qubits == [1]
    @test isa_variant(c[21], CircuitOp.Measurement) && c[21].bit == 2 && c[21].qubits == [2]
    @test isa_variant(c[22], CircuitOp.Measurement) && c[22].bit == 3 && c[22].qubits == [3]
end

@testset "barrier is skipped" begin
    with_qasm("barrier q;\nx q[0];") do path
        c = parse_input(path)
        @test length(c) == 4  # 1 gate + 3 measure-all
        @test isa_variant(c[1], CircuitOp.ExpHalfPiPauli)
    end
end

@testset "Gate sequence ordering is preserved" begin
    with_qasm("s q[0];\nt q[0];\nsdg q[0];") do path
        c = parse_input(path)
        @test length(c) == 6  # 3 gates + 3 measure-all
        @test isa_variant(c[1], CircuitOp.ExpQuatPiPauli) && c[1].pauli == P"Z"   # s
        @test isa_variant(c[2], CircuitOp.ExpEighPiPauli) && c[2].pauli == P"Z"   # t
        @test isa_variant(c[3], CircuitOp.ExpQuatPiPauli) && c[3].pauli == -P"Z"  # sdg
    end
end

@testset "Measure-all" begin
    with_qasm("x q[0];") do path
        c = parse_input(path)
        # 1 gate op + 3 measure-all (qubit q → classical bit q)
        @test length(c) == 4
        for q in 1:3
            @test isa_variant(c[1+q], CircuitOp.Measurement)
            @test c[1+q].pauli == P"Z" && c[1+q].bit == q && c[1+q].qubits == [q]
        end
    end
    with_qasm("x q[0];\nmeasure q[0] -> c[0];") do path
        c = parse_input(path)
        # explicit measure suppresses measure-all
        @test length(c) == 2
        @test isa_variant(c[2], CircuitOp.Measurement) && c[2].bit == 1
    end
end

@testset "Sample input file (sample_input.qasm)" begin
    c = parse_input(joinpath(FIXTURES,"sample_input.qasm"))
    cx_ops      = filter(op -> isa_variant(op, CircuitOp.PauliConditional), c)
    measure_ops = filter(op -> isa_variant(op, CircuitOp.Measurement),      c)
    @test length(cx_ops)      == 3
    @test length(measure_ops) == 2
    # Measurements are in Z basis on qubits 0 and 1
    @test all(op.pauli == P"Z" for op in measure_ops)
    @test measure_ops[1].qubits == [1] && measure_ops[1].bit == 1
    @test measure_ops[2].qubits == [2] && measure_ops[2].bit == 2
end

##
# ---------------------------------------------------------------------------
# Helper: construct a minimal CompilerState for coverage testing
# ---------------------------------------------------------------------------
function _make_compiler_state()
    return CompilerState(
        measurement_results = MeasurementResult.Type[ClassicalDetermRes(P"XZ", true)],
        stabilizer_group    = Stabilizer([P"XX", P"ZZ"]),
        classical_register  = Union{Nothing,Bool}[true, nothing],
        circuit             = Circuit(),
        instruction_pointer = 2,
        runtime             = SimRuntime(),
    )
end

# ---------------------------------------------------------------------------
# Helper: construct a minimal CompilationResult for testing
# ---------------------------------------------------------------------------
function _make_result(;
    measurement_results = MeasurementResult.Type[
        ClassicalDetermRes(P"XZ", true),
        ClassicalRandomRes(P"ZX", false),
        QuantumRes(P"YI", true),
    ],
    qpu_workload = MeasurementResult.Type[
        QuantumRes(P"YI", true),
    ],
    stabilizer_group = Stabilizer([P"XX", P"ZZ"]),
    qpu_duration = 5,
)
    return CompilationResult(measurement_results, qpu_workload, stabilizer_group, qpu_duration)
end

# ---------------------------------------------------------------------------

# Tests that save_result(CompilerState) creates a file with the expected keys.
# No round-trip: CompilerState has no load_result counterpart.
@testset "save_result(CompilerState) creates file with correct keys" begin
    path = tempname() * ".jld2"
    save_result(_make_compiler_state(), path)
    @test isfile(path)
    data = JLD2.load(path)
    @test haskey(data, "measurement_results")
    @test haskey(data, "stabilizer_group")
    @test haskey(data, "classical_register")
    @test haskey(data, "instruction_pointer")
    @test !haskey(data, "circuit")
    @test !haskey(data, "runtime")
    rm(path)
end

# Tests that `save_result` creates a file at the exact path supplied.
@testset "save_result creates file at specified path" begin
    path = tempname() * ".jld2"
    save_result(_make_result(), path)
    @test isfile(path)
    rm(path)
end

# Tests that the saved file contains exactly the four expected top-level keys
# and that no CompilerState internals (circuit, classical_register, etc.) leak in.
@testset "save_result writes correct keys" begin
    path = tempname() * ".jld2"
    save_result(_make_result(), path)
    data = JLD2.load(path)
    @test haskey(data, "measurement_results")
    @test haskey(data, "QPU_workload")
    @test haskey(data, "stabilizer_group")
    @test haskey(data, "QPUDuration")
    @test !haskey(data, "classical_register")
    @test !haskey(data, "circuit")
    rm(path)
end

# Tests that load_result() reconstructs a CompilationResult identical to what was saved,
# covering measurement_results round-trip across all three variant types.
@testset "load_result round-trips measurement_results" begin
    path = tempname() * ".jld2"
    original = _make_result()
    save_result(original, path)
    loaded = load_result(path)

    @test length(loaded.measurement_results) == 3

    @test loaded.measurement_results[1].pauli  == P"XZ"
    @test loaded.measurement_results[1].result == true
    @test isa_variant(loaded.measurement_results[1], ClassicalDetermRes)

    @test loaded.measurement_results[2].pauli  == P"ZX"
    @test loaded.measurement_results[2].result == false
    @test isa_variant(loaded.measurement_results[2], ClassicalRandomRes)

    @test loaded.measurement_results[3].pauli  == P"YI"
    @test loaded.measurement_results[3].result == true
    @test isa_variant(loaded.measurement_results[3], QuantumRes)
    rm(path)
end

# Tests that load_result() correctly restores QPU_workload entries.
@testset "load_result round-trips QPU_workload" begin
    path = tempname() * ".jld2"
    save_result(_make_result(), path)
    loaded = load_result(path)

    @test length(loaded.QPU_workload) == 1
    @test loaded.QPU_workload[1].pauli  == P"YI"
    @test loaded.QPU_workload[1].result == true
    @test isa_variant(loaded.QPU_workload[1], QuantumRes)
    rm(path)
end

# Tests that the stabilizer_group round-trips through the JLD2 file,
# checking both the number of generators and their individual Pauli strings.
@testset "load_result round-trips stabilizer_group" begin
    stab = Stabilizer([P"XX", P"ZZ"])
    path = tempname() * ".jld2"
    save_result(_make_result(stabilizer_group=stab), path)
    loaded = load_result(path)

    @test length(loaded.stabilizer_group) == 2
    @test loaded.stabilizer_group[1] == P"XX"
    @test loaded.stabilizer_group[2] == P"ZZ"
    rm(path)
end

# Tests that QPUDuration round-trips as an exact integer value.
@testset "load_result round-trips QPUDuration" begin
    path = tempname() * ".jld2"
    save_result(_make_result(qpu_duration=42), path)
    loaded = load_result(path)

    @test loaded.QPUDuration == 42
    rm(path)
end

end
