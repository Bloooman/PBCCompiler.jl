using BenchmarkTools
using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, Pauli, Measurement, ExpHalfPiPauli, traversal,
    preprocess_circuit, random_test_circuit, get_distribution, SimRuntime
using QuantumClifford: @P_str, random_pauli
using Moshi.Data: isa_variant
using Random: MersenneTwister, seed!, randperm

const SUITE = BenchmarkGroup()

# Helper function to create a circuit with n Pauli gates
function make_pauli_circuit(n::Int)
    paulis = [P"X", P"Y", P"Z", P"I"]
    ops = [Pauli(paulis[mod1(i, 4)], [1]) for i in 1:n]
    return Circuit(ops)
end

# Helper function to create a mixed circuit with different gate types
function make_mixed_circuit(n::Int)
    paulis = [P"X", P"Y", P"Z"]
    ops = CircuitOp.Type[]
    for i in 1:n
        p = paulis[mod1(i, 3)]
        if i % 4 == 0
            push!(ops, Measurement(p, i, [1]))
        elseif i % 4 == 1
            push!(ops, Pauli(p, [1]))
        elseif i % 4 == 2
            push!(ops, ExpHalfPiPauli(p, [1]))
        else
            push!(ops, Pauli(p, [1, 2]))
        end
    end
    return Circuit(ops)
end

# Transformation functions for benchmarks
swap_transform(op1, op2) = (op2, op1)
noop_transform(op1, op2) = nothing
function combine_paulis(op1, op2)
    if isa_variant(op1, CircuitOp.Pauli) && isa_variant(op2, CircuitOp.Pauli)
        return Pauli(P"X", [1])
    end
    return nothing
end

# Traversal benchmarks
SUITE["traversal"] = BenchmarkGroup(["traversal"])

# Swap traversal - moves operations around without changing circuit length
SUITE["traversal"]["swap"] = BenchmarkGroup(["swap"])
SUITE["traversal"]["swap"]["10"] = @benchmarkable traversal(c, swap_transform) setup=(c=make_pauli_circuit(10)) evals=1
SUITE["traversal"]["swap"]["100"] = @benchmarkable traversal(c, swap_transform) setup=(c=make_pauli_circuit(100)) evals=1
SUITE["traversal"]["swap"]["1000"] = @benchmarkable traversal(c, swap_transform) setup=(c=make_pauli_circuit(1000)) evals=1

# No-op traversal - visits all pairs but makes no changes
SUITE["traversal"]["noop"] = BenchmarkGroup(["noop"])
SUITE["traversal"]["noop"]["10"] = @benchmarkable traversal(c, noop_transform) setup=(c=make_pauli_circuit(10)) evals=1
SUITE["traversal"]["noop"]["100"] = @benchmarkable traversal(c, noop_transform) setup=(c=make_pauli_circuit(100)) evals=1
SUITE["traversal"]["noop"]["1000"] = @benchmarkable traversal(c, noop_transform) setup=(c=make_pauli_circuit(1000)) evals=1

# Combine traversal - reduces circuit by combining adjacent Paulis
SUITE["traversal"]["combine"] = BenchmarkGroup(["combine"])
SUITE["traversal"]["combine"]["10"] = @benchmarkable traversal(c, combine_paulis) setup=(c=make_pauli_circuit(10)) evals=1
SUITE["traversal"]["combine"]["100"] = @benchmarkable traversal(c, combine_paulis) setup=(c=make_pauli_circuit(100)) evals=1
SUITE["traversal"]["combine"]["1000"] = @benchmarkable traversal(c, combine_paulis) setup=(c=make_pauli_circuit(1000)) evals=1

# Mixed circuit traversal
SUITE["traversal"]["mixed"] = BenchmarkGroup(["mixed"])
SUITE["traversal"]["mixed"]["swap_100"] = @benchmarkable traversal(c, swap_transform) setup=(c=make_mixed_circuit(100)) evals=1
SUITE["traversal"]["mixed"]["noop_100"] = @benchmarkable traversal(c, noop_transform) setup=(c=make_mixed_circuit(100)) evals=1

# Direction comparison
SUITE["traversal"]["direction"] = BenchmarkGroup(["direction"])
SUITE["traversal"]["direction"]["right_100"] = @benchmarkable traversal(c, swap_transform, :right) setup=(c=make_pauli_circuit(100)) evals=1
SUITE["traversal"]["direction"]["left_100"] = @benchmarkable traversal(c, swap_transform, :left) setup=(c=make_pauli_circuit(100)) evals=1

# ---------------------------------------------------------------------------
# Pipeline benchmarks
#
# The traversal benchmarks above drive synthetic stub kernels, so they are blind
# to the cost of the real conjugation kernels and of execution. These exercise
# the actual pipeline. Every circuit is generated from a seeded RNG, and the
# execution benchmarks seed the global RNG in `setup` as well: `run` takes a
# coin-flip branch that splices compensating rotations into the circuit, so
# without a seed the *amount of work* varies between samples, not just the
# timing noise.
# ---------------------------------------------------------------------------

const BENCH_SEED = 20260806

"""
Clifford-heavy circuit: `nops` rotations of which `nT` are pi/8. This is the
Clifford+T shape that comes out of QASM input — a long Clifford run with a
sparse scattering of T gates — and it is the regime where compilation, rather
than the chi-expansion of the magic register, dominates the per-shot cost.
`random_test_circuit` is T-heavy (a quarter of its ops are pi/8) and covers the
opposite regime.
"""
function make_clifford_heavy(nops::Int, nq::Int, nT::Int; seed::Int=BENCH_SEED)
    rng = MersenneTwister(seed)
    qubits = collect(1:nq)
    tpos = Set(randperm(rng, nops)[1:nT])
    ops = CircuitOp.Type[]
    for i in 1:nops
        p = random_pauli(rng, nq; nophase=true)
        while iszero(p.xz)
            p = random_pauli(rng, nq; nophase=true)
        end
        push!(ops, i in tpos ? CircuitOp.ExpEighPiPauli(p, qubits) :
                   isodd(i)  ? CircuitOp.ExpQuatPiPauli(p, qubits) :
                               CircuitOp.ExpHalfPiPauli(p, qubits))
    end
    for i in 1:nq
        push!(ops, Measurement(P"Z", i, [i]))
    end
    return Circuit(ops)
end

# Compilation only. Scaling here is dominated by `group_nonclifford`, which
# bubbles each pi/8 rotation left past every preceding Clifford.
SUITE["preprocess"] = BenchmarkGroup(["preprocess"])
for n in (160, 320, 640)
    SUITE["preprocess"]["random_$n"] =
        @benchmarkable preprocess_circuit(c) setup=(c=random_test_circuit($n, 4; rng=MersenneTwister(BENCH_SEED))) evals=1
end
for n in (200, 400, 800)
    SUITE["preprocess"]["cliffordT_$n"] =
        @benchmarkable preprocess_circuit(c) setup=(c=make_clifford_heavy($n, 4, 6)) evals=1
end

# Compile plus a single execution.
SUITE["run"] = BenchmarkGroup(["run"])
SUITE["run"]["cliffordT_400"] =
    @benchmarkable PBCCompiler.run(c, SimRuntime(), nothing) setup=(seed!(BENCH_SEED); c=make_clifford_heavy(400, 4, 6)) evals=1
SUITE["run"]["Theavy_200"] =
    @benchmarkable PBCCompiler.run(c, SimRuntime(), nothing) setup=(seed!(BENCH_SEED); c=make_clifford_heavy(200, 4, 20)) evals=1

# Shot sampling: this is what the research samplers actually spend their time
# in, and where amortizing compilation across shots shows up.
SUITE["sample"] = BenchmarkGroup(["sample"])
SUITE["sample"]["cliffordT_400_50shots"] =
    @benchmarkable get_distribution(c, SimRuntime(), nothing, 50) setup=(seed!(BENCH_SEED); c=make_clifford_heavy(400, 4, 6)) evals=1
SUITE["sample"]["Theavy_200_50shots"] =
    @benchmarkable get_distribution(c, SimRuntime(), nothing, 50) setup=(seed!(BENCH_SEED); c=make_clifford_heavy(200, 4, 20)) evals=1
