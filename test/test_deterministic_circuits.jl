@testitem "Deterministic non-Clifford circuits" tags=[:deterministic] begin

# Truth-table tests: deterministic circuits containing non-Clifford gates.
# The compiled circuit still takes random branches at every gadget measurement,
# so a bit-exact deterministic output on every shot verifies that the Pauli
# correction machinery compensates each random branch correctly. This is the
# regression test for the stabilizer-row sign bug: the recorded row must carry
# the measured (-1)^result sign or the target bit flips on random shots.

using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, ExpHalfPiPauli, parse_input, SimRuntime
using QuantumClifford: @P_str

const FIXTURES = joinpath(@__DIR__, "fixtures")

xgate(q) = ExpHalfPiPauli(P"X", [q])

# Only SimRuntime is a faithful simulator. DummyRuntime coins every quantum
# outcome at a fixed bias, but later gadget measurements have true conditional
# probabilities other than 1/2 (the magic memory is already partially
# projected), so it can take physically impossible branches that no Pauli
# correction can compensate — deterministic circuits do NOT come out
# deterministic under DummyRuntime, by design.
@testset "Toffoli truth table (SimRuntime)" begin
    ccx = parse_input(joinpath(FIXTURES, "toffoli3.qasm"))
    for a in 0:1, b in 0:1, c in 0:1
        prep = CircuitOp.Type[]
        a == 1 && push!(prep, xgate(1))
        b == 1 && push!(prep, xgate(2))
        c == 1 && push!(prep, xgate(3))
        expected = Bool[a, b, xor(c, a & b)]
        for shot in 1:3
            circuit = Circuit([prep; ccx])
            state = PBCCompiler.run(circuit, SimRuntime())
            got = collect(Bool, state.classical_register[1:3])
            @test got == expected
        end
    end
end

end
