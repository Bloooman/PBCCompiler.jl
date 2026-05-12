# PBCCompiler.jl

```@meta
DocTestSetup = quote
    using QuantumClifford, PBCCompiler
end
```

Compiles arbitrary gate-based quantum circuits into the **Pauli-Based Computation** model,
Clifford operations into Pauli product rotations, inject gadgets with magic states for
non-Clifford gates and Meausurements into Pauli product measurements. Optionally performs simulation and returns measurement statistics alongside the compiled Pauli product measurement circuit.

## Overview

PBCCompiler.jl provides infrastructure for compiling and executing Pauli-based quantum circuits. The package includes:

- **Circuit Operations** - Algebraic data types for representing quantum operations
- **Compilation Pipeline** - Transform circuits through optimization stages
- **Runtime System** - Execute circuits on quantum backends

---

## Installation

```julia
] add PBCCompiler
```

Or explicitly:

```julia
using Pkg; Pkg.add("PBCCompiler")
```

Requires Julia ≥ 1.9.

---

## Circuit Operations

The `CircuitOp` type represents quantum operations:

- `Measurement` - Pauli measurement with classical bit output
- `Pauli` - Pauli gate application
- `ExpHalfPiPauli`, `ExpQuatPiPauli`, `ExpEighPiPauli` - Rotation gates
- `PrepMagic` - Magic state preparation
- `PauliConditional` - Pauli gate conditioned on another Pauli
- `BitConditional` - Operation conditioned on classical bit

---

## Quickstart

```julia
using PBCCompiler

input_circuit = Circuit([
ExpQuatPiPauli(P"Z", [1]),
ExpQuatPiPauli(P"X", [1]),
ExpQuatPiPauli(P"Z", [1]),
PauliConditional(P"Z", [1], P"X", [2]),
ExpEighPiPauli(P"Z", [2]),
Measurement(P"Z", 1, [1]),
Measurement(P"Z", 2, [2])
])

input_state = S"+Z_
                +_Z"

result = run(input_circuit,input_state)
```

---

## Background

### Pauli-Based Computation

PBC is a model of quantum computation in which all operations are Pauli product rotations

$$e^{-i\frac{\theta}{2} P}, \quad P \in \{I, X, Y, Z\}^{\otimes n}$$

Clifford gates correspond to $\theta \in \{\pi/2, \pi\}$ and are classically tractable
(Gottesman–Knill). Non-Clifford gates require **gadget injection**: an ancilla qubit
prepared in a magic state $|M\rangle = T|+\rangle$ is consumed via an adaptive
measurement, with a classically-controlled Pauli correction applied on the outcome.

### Compilation strategy

1. **Clifford extraction** — decompose the input into a Clifford layer (Pauli rotations
   at Clifford angles) and a residual non-Clifford part.
2. **Gadget insertion** — replace each non-Clifford gate with a measurement-based gadget
   consuming one magic state per $T$-gate equivalent.
3. **Output** — a `PBCCircuit` containing an ordered sequence of `PauliGadget` operations
   and ancilla resources.

---

## API Reference

### Functions

#### `compile`

```julia
compile(circuit::String; gadget_strategy::Symbol=:auto) → PBCCircuit
compile(circuit::GateCircuit; gadget_strategy::Symbol=:auto) → PBCCircuit
```

Compiles a gate-based circuit into the PBC model. Accepts either an OpenQASM 3 string
or a `GateCircuit`. Returns a `PBCCircuit`.

`gadget_strategy` controls non-Clifford gate handling:
- `:auto` — infers optimal gadget type per gate (default)
- `:t_gate` — forces T-gate gadget injection
- `:ccz` — forces CCZ gadget decomposition

---

#### `simulate`

```julia
simulate(circuit::String; shots::Int=1024) → MeasurementResult
simulate(circuit::GateCircuit; shots::Int=1024) → MeasurementResult
simulate(circuit::PBCCircuit; shots::Int=1024) → MeasurementResult
```

Compiles (if needed) and simulates the circuit. The third dispatch accepts a
pre-compiled `PBCCircuit` directly, skipping recompilation.

Returns a `MeasurementResult`.

---

#### `to_qasm`

```julia
to_qasm(circuit::PBCCircuit) → String
```

Serializes a compiled `PBCCircuit` to OpenQASM 3.

---

### Types

#### `GateCircuit`

```julia
GateCircuit(n_qubits::Int)
```

Native gate-based circuit IR. Build programmatically if not using QASM input.

| Field | Type | Description |
|---|---|---|
| `n_qubits` | `Int` | Number of logical qubits |
| `gates` | `Vector{Gate}` | Ordered gate sequence |

```julia
qc = GateCircuit(3)
push!(qc, HGate(1))
push!(qc, CXGate(1, 2))
push!(qc, TGate(3))
```

---

#### `PauliGadget`

```julia
PauliGadget(pauli_string::String, angle::Float64)
```

Represents $e^{-i\frac{\theta}{2}P}$. Users construct these directly only when
building `PBCCircuit`s manually; normally produced by `compile`.

| Field | Type | Description |
|---|---|---|
| `pauli_string` | `String` | Pauli operator string over `{I,X,Y,Z}`, e.g. `"XZIY"` |
| `angle` | `Float64` | Rotation angle $\theta$ in radians |

**Invariant**: `length(pauli_string)` must equal the circuit's `n_qubits`.

---

#### `PBCCircuit` *(returned by `compile`)*

Not constructed directly. Inspect via:

| Field | Type | Description |
|---|---|---|
| `n_qubits` | `Int` | Logical qubit count (excluding ancillae) |
| `n_ancilla` | `Int` | Ancilla qubits added for gadget injection |
| `gadgets` | `Vector{PauliGadget}` | Compiled Pauli product rotations, in order |
| `non_clifford_count` | `Int` | Number of injected non-Clifford gadgets |

---

#### `MeasurementResult` *(returned by `simulate`)*

Not constructed directly. Inspect via:

| Field | Type | Description |
|---|---|---|
| `counts` | `Dict{String,Int}` | Bitstring → count histogram |
| `shots` | `Int` | Total shots executed |
| `circuit` | `PBCCircuit` | Compiled circuit used in simulation |

---

## Output format

`to_dict(pbc::PBCCircuit)` returns a serializable structure:

```json
{
  "n_qubits": 3,
  "n_ancilla": 1,
  "gadgets": [
    {"pauli": "XXI", "angle": 1.5707963267948966},
    {"pauli": "ZII", "angle": 3.141592653589793}
  ],
  "non_clifford_gadgets": [
    {"type": "T_gadget", "target": 2, "ancilla": 3, "correction": "Z"}
  ]
}
```

---

## Examples

See [`examples/`](examples/) for notebooks covering:

- Bernstein–Vazirani compiled to PBC
- Toffoli gate gadget injection
- QFT gadget count scaling with $n$

---

## Contributing

```bash
git clone https://github.com/you/PBCCompiler.jl
cd PBCCompiler.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. test/runtests.jl
```

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for conventions and PR guidelines.

---

## Citation

```bibtex
@software{pbccompiler_jl,
  author  = {Your Name},
  title   = {{PBCCompiler.jl}: A gate-to-PBC quantum circuit compiler},
  year    = {2025},
  url     = {https://github.com/you/PBCCompiler.jl}
}
```

---

## License

GPL-3.0. See [LICENSE](LICENSE) for details.

## Related Packages

- [QuantumClifford.jl](https://github.com/QuantumSavory/QuantumClifford.jl) - Stabilizer formalism
- [BPGates.jl](https://github.com/QuantumSavory/BPGates.jl) - Bell-preserving gates
