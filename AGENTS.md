# PBCCompiler.jl

Tools for Pauli Based Computation (PBC), a modality of quantum computation.

## Project Structure

- `src/PBCCompiler.jl` - Main module with circuit operations and compiler infrastructure
- `src/type.jl` - `CircuitOp`/`MeasurementResult` ADTs, runtimes, `CompilerState`
- `src/traversal.jl` - Circuit traversal utilities for gate simplifications
- `src/pair_transformation.jl` - The 2->2 kernels driven by `traversal` (conjugation, merging, commutation)
- `src/preprocess.jl` - The `preprocess_circuit` compilation pipeline
- `src/joint_measurement_check.jl` - Measurement outcomes against the tableau and magic register
- `src/logic.jl` - `build_compilerstate`, `execute!`, `run`, `to_result`
- `src/statistics.jl` - Shot sampling and interaction-graph extraction
- `src/io.jl` - QASM parsing, result save/load
- `src/random_circuit.jl` - Random circuit generation
- `src/affectedqubits.jl` - Query functions for qubit indices affected by operations
- `src/plotting.jl` - Plotting function stubs (scaffolding for extensions)
- `src/partition.jl` - Branch-enumeration search. **Not in the module's include
  list**, i.e. currently dead code
- `ext/PBCCompilerMakieExt/` - Makie extension for circuit visualization
- `test/` - Test suite using TestItemRunner.jl
- `benchmark/` - Performance benchmarks using BenchmarkTools.jl

## Dependencies

- **Moshi.jl** - Algebraic data types via `@data` macro
- **QuantumClifford.jl** - Pauli operators and stabilizer formalism

## Key Concepts

### Circuit Operations (`CircuitOp`)
Algebraic data type representing quantum circuit operations:
- `Measurement` - Pauli measurement with classical bit output
- `Pauli` - Pauli gate application
- `ExpHalfPiPauli`, `ExpQuatPiPauli`, `ExpEighPiPauli` - Rotation gates (pi/2, pi/4, pi/8)
- `PrepMagic` - Magic state preparation
- `PauliConditional` - Pauli gate conditioned on another Pauli
- `BitConditional` - Operation conditioned on classical bit

### Compilation Pipeline
The `preprocess_circuit` function transforms circuits through stages:
1. Remove Pauli conditionals
2. Commute non-Clifford gates to front
3. Group non-Clifford operations
4. Commute measurements to end
5. Remove non-Clifford gates (introduce magic states)
6. Remove post-measurement operations

### Runtime
- `AbstractRuntime` - Supertype of the measurement backends
- `SimRuntime` - Simulates the magic register with a `GeneralizedStabilizer`;
  outcomes that anticommute with the stabilizer group become coin flips resolved
  by splicing compensating rotations into the circuit
- `StabilizerRuntime` - Simulates the full register (data + magic) together, so
  it projects for every non-deterministic outcome and never yields a
  `ClassicalRandomRes`. **This is intended design — do not "fix" it.** Every
  non-deterministic outcome becomes a `QuantumRes`, including ones that are
  random only because of the data register, so `QPU_workload` runs far larger
  than `SimRuntime`'s (~5x on random small circuits) and the two are not
  comparable. Outcomes themselves agree between the runtimes; only the
  classification differs
- `DummyRuntime` / `TraversalRuntime` - Replace quantum measurements with
  classical coin flips of a fixed bias. `DummyRuntime` still tracks
  `activated` the same way `SimRuntime` does, for parity/diagnostics — it
  isn't consumed by `to_result`, which derives its magic block from the
  `QuantumRes` count instead
- `DummyStabilizerRuntime` - Cheap stand-in for `StabilizerRuntime`: same
  control flow (no anticommuting coin-flip branch) and same `activated`
  bookkeeping of which magic qubits a measurement touched, but coin-flips
  every non-deterministic outcome instead of simulating the full register.
  `to_result`/`QPU_workload` extraction work the same way they do for
  `StabilizerRuntime`; only the outcome bias is not physically faithful
- `HybridRuntime` - Starts out simulating the magic register like
  `SimRuntime`, then converts in place into a `StabilizerRuntime` once
  `maximum_measurement_support` activated magic qubits have been reached.
  `PBCCompiler.run` drives the conversion by calling `transition` after every
  measurement step. With `maximum_measurement_support = nothing` (the
  default) it never converts and behaves exactly like `SimRuntime`
- `CompilerState{R,T}` - Tracks measurement results, tableau, classical
  register, circuit, instruction pointer and runtime. Parameterized on the
  runtime and tableau types to keep the execution loop type stable

Runtimes carry mutable contents (`quantum_memory`, `activated`,
`invsparsity_history`), so `copy(::CompilerState)` copies them. Anything that
derives two states from one compilation must go through `copy`, never share a
runtime.

### Compiling once, running many shots
`PBCCompiler.run` is `build_compilerstate` followed by `execute!`. Compilation
is shot-independent, so sampling loops compile once and run each shot off a
copy:

```julia
compiled = build_compilerstate(circuit, SimRuntime())
shots = map(1:1000) do _
    s = copy(compiled)
    while !PBCCompiler._execution_complete(s)
        s = execute!(s)
    end
    s
end
```

`execute!` performs one measurement step per call (a no-op once every
measurement is done), so driving a shot to completion means looping it as
above. `s.circuit` itself never shrinks -- looping on `!isempty(s.circuit)`
instead of the measurement-completion check spins forever.

`get_distribution` and `weight_std_graph` already do this.

**Reproducibility**: outcomes come from the global RNG (`Random.seed!` makes a
run reproducible). There is no `rng` keyword on `run`: `QuantumClifford`'s
`projectrand!` takes no RNG argument, so one cannot be threaded all the way
through without an upstream change. `random_test_circuit` does accept `rng`.

### Circuit Traversal
The `traversal` function (`src/traversal.jl`) applies transformations to adjacent pairs of circuit operations:
```julia
traversal(circuit, pair_transformation, direction=:right, starting_index=1, end_index=:end)
```
- `pair_transformation(op1, op2)` returns:
  - `(new_op1, new_op2)` tuple to replace the pair
  - Single operation to combine the pair into one
  - `nothing` to keep unchanged
- Supports left-to-right (`:right`) or right-to-left (`:left`) traversal
- Used for gate commutation, simplification, and compilation passes

**Note on Moshi types**: Use `Moshi.Data.isa_variant(op, CircuitOp.Pauli)` instead of `op isa CircuitOp.Pauli` to check variant types.

### Affected Qubits
The `affectedqubits` function (`src/affectedqubits.jl`) returns the sorted list of qubit indices affected by an operation or circuit:
```julia
affectedqubits(op::CircuitOp.Type) -> Vector{Int}
affectedqubits(circuit::Circuit) -> Vector{Int}
```

### Circuit Visualization (Makie Extension)
When Makie is loaded, the `PBCCompilerMakieExt` extension provides circuit plotting:
```julia
using CairoMakie  # or GLMakie
using PBCCompiler

circuit = Circuit([...])
circuitplot(circuit)  # Create a plot
circuitplot!(ax, circuit)  # Add to existing axis
circuitplot_axis(fig[1,1], circuit)  # Create complete figure panel
```

**Plot features:**
- Gates shown as colored rectangles spanning affected qubits
- Horizontal qubit wire lines
- Measurement results marked with classical bit index (e.g., "c0")
- Conditional operations marked with dependency index (e.g., "?c0")
- PrepMagic gates not visualized (placeholder for future)

**Configurable attributes:**
- `gatewidth`, `qubitspacing` - Gate dimensions
- `wirecolor`, `wirelinewidth` - Wire appearance
- `paulicolor`, `measurementcolor`, etc. - Gate colors by variant

### Circuit Visualization (Quantikz Extension)
When Quantikz is loaded, the `PBCCompilerQuantikzExt` extension provides
publication-quality LaTeX diagrams:
```julia
using Quantikz
using PBCCompiler

circuitplot_quantikz(circuit)                      # rendered image
circuitplot_quantikz(circuit, "circuit.pdf")       # or .png / .tex
circuitstring_quantikz(circuit)                    # raw LaTeX, no toolchain needed
```

**Why this and not the Makie/Qiskit renderers:** both of those draw a multi-qubit
operation as *one rectangle spanning min→max qubit*, which covers wires the
operation does not act on and hides the per-qubit Pauli letter. Qiskit's mpl
drawer has no public API for anything else. The Quantikz renderer instead draws
one labelled box per supported qubit, joined by a vertical `\vqw` rail
(Litinski-style), so an operation's support is readable straight off the diagram.

**Implementation notes:**
- `PauliBoxOp <: Quantikz.QuantikzOp` emits the per-wire cells. Quantikz's own
  `MultiControlU`/`Measurement`/`ParityMeasurement` are deliberately *not* reused:
  the first two draw spanning rectangles and the latter two terminate wires via
  `deleteoutputs`, but PBC Pauli measurements are non-destructive (`\meter`-style
  boxes that keep the wire, not `\meterD`).
- The environment options must not include Quantikz.jl's default `transparent`
  key — it suppresses every `fill=` the renderer emits. `circuitstring_quantikz`
  therefore passes its own `quantikzoptions`.
- Column packing comes free from `Quantikz.circuit2table_compressed`; no
  hand-rolled layering pass is needed.
- `CircuitOp.Pauli` and `CircuitOp.PrepMagic` are not handled (same coverage as
  `circuitplot_qiskit`).

## Development

### Workflow
1. Always pull latest master: `git pull`
2. Create feature branches for new work
3. Commit often at each change
4. Update CLAUDE.md with new functionality
5. Run tests before creating PRs
6. Add benchmarks for new performance-critical functionality

### Docstring Guidelines
- Docstrings are for **users**, not developers
- Do not include implementation details (e.g., "uses pattern matching", "implemented via recursion")
- Focus on: what the function does, its arguments, return values, and usage examples
- Implementation notes belong in code comments, not docstrings

### Run tests
```bash
julia -tauto --project -e 'using Pkg; Pkg.test("PBCCompiler")'
```

### Benchmarks
Benchmarks are managed with BenchmarkTools.jl and run in CI via AirspeedVelocity.jl.

Run benchmarks locally:
```bash
julia -tauto --project=benchmark -e 'include("benchmark/benchmarks.jl"); run(SUITE)'
```

**When to add benchmarks:**
- New compilation passes or transformations
- New traversal operations
- Any function that processes circuits at scale
- Performance-critical code paths

**Benchmark file structure:**
- Add new benchmarks to `benchmark/benchmarks.jl`
- Use `evals=1` for functions that modify state in-place
- Use `setup=` to create fresh data for each evaluation
- Group related benchmarks using `BenchmarkGroup`

### Julia invocation
Always use the `-tauto` flag when launching Julia to utilize all available threads, which drastically speeds up compilation times:
```bash
julia -tauto --project
```

### Related source code
- QuantumClifford.jl source: `../QuantumClifford.jl`
- QuantumInterface.jl source: `../QuantumInterface.jl`

### Reference paper
- "Game of Surface Codes" - https://quantum-journal.org/papers/q-2019-03-05-128/pdf/

## Related Packages

- [QuantumClifford.jl](https://github.com/QuantumSavory/QuantumClifford.jl) - Stabilizer formalism
- [BPGates.jl](https://github.com/QuantumSavory/BPGates.jl) - Bell-preserving gates
