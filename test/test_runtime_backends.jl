@testitem "Runtime backends: measurement classification" tags=[:deterministic] begin

# `StabilizerRuntime` classifies each measurement purely from the return value of
# `project!`, without the anticommutation pre-scan that the `SimRuntime` method
# performs. `project!(::MixedDestabilizer, p)` returns `(state, anticom, result)`
# with `result === nothing` in two physically different situations:
#
#   1. `p` anticommutes with a stabilizer row  -> genuinely random outcome
#   2. `p` commutes but is independent         -> rank grows (the gadget case)
#
# `result !== nothing` happens only when `anticom == 0`. These tests pin down that
# the two situations are told apart, and that the result-type sequence agrees with
# the faithful `SimRuntime` backend.

using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, Measurement, ExpHalfPiPauli, ExpEighPiPauli,
    PauliConditional, SimRuntime, StabilizerRuntime, MeasurementResult
using Moshi.Data: variant_name
using QuantumClifford: @P_str

"""Result-type sequence (`:QuantumRes`, `:ClassicalRandomRes`, ...) of a run."""
result_types(state) = [variant_name(m) for m in state.measurement_results]

# A Clifford-only circuit allocates no magic qubits, so any path that reaches the
# magic register is a misclassification: `X` anticommutes with the `+Z` stabilizer
# of |0>, so this is a coin-flip outcome, not a gadget measurement.
@testset "Clifford-only random measurement needs no magic state" begin
    circuit() = Circuit(CircuitOp.Type[Measurement(P"X", 1, [1])])
    @test result_types(PBCCompiler.run(circuit(), SimRuntime())) == [:ClassicalRandomRes]
    @test result_types(PBCCompiler.run(circuit(), StabilizerRuntime())) == [:ClassicalRandomRes]
end

# Same circuit as the "SimRuntime Correctness" testset in test_compilation_pass.jl.
# Measurement 2 anticommutes with the group and must stay a classical coin flip;
# only measurement 1 (the gadget) is a true quantum measurement.
@testset "result types agree with SimRuntime" begin
    circuit() = Circuit(CircuitOp.Type[
        PauliConditional(P"X", [1], P"Z", [2]),
        ExpHalfPiPauli(P"YZ", [1, 2]),
        ExpEighPiPauli(P"Z", [2]),
        Measurement(P"Z", 1, [1]),
        Measurement(P"Z", 2, [2]),
    ])
    expected = [:QuantumRes, :ClassicalRandomRes, :ClassicalDetermRes, :ClassicalDetermRes]
    @test result_types(PBCCompiler.run(circuit(), SimRuntime())) == expected
    @test result_types(PBCCompiler.run(circuit(), StabilizerRuntime())) == expected
end

end
