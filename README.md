# PBCCompiler.jl

> Gate-based quantum circuit compiler targeting Pauli-Based Computation (PBC).

[![Build Status](https://github.com/you/PBCCompiler.jl/actions/workflows/CI.yml/badge.svg)]()
[![Coverage](https://codecov.io/gh/you/PBCCompiler.jl/branch/main/graph/badge.svg)]()
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://you.github.io/PBCCompiler.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://you.github.io/PBCCompiler.jl/dev)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Compiles arbitrary gate-based quantum circuits into the **Pauli-Based Computation** model,
Clifford operations into Pauli product rotations, inject gadgets with magic states for
non-Clifford gates and Meausurements into Pauli product measurements. Optionally performs simulation and returns measurement statistics alongside the compiled Pauli product measurement circuit.

Full documentation at **[you.github.io/PBCCompiler.jl/stable](https://you.github.io/PBCCompiler.jl/stable)**.

---

## Installation

```julia
] add PBCCompiler
```

Requires Julia ≥ 1.9.

---

## Quick example

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

## Citation

```bibtex
@software{pbccompiler_jl,
  author = {Yuhan Gao},
  title  = {{PBCCompiler.jl}: A gate-to-PBC quantum circuit compiler},
  year   = {2025},
  url    = {https://github.com/Bloooman/PBCCompiler.jl}
}
```

## License

GPL-3.0. See [LICENSE](LICENSE) for details.
