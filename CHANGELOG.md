# News

## v0.1.0 - dev

- **Fix**: `copy(::CompilerState)` now copies the runtime's mutable contents
  (magic-state memory, `activated` bits, `invsparsity_history`) instead of
  aliasing them. Two states derived from one compilation previously shared a
  magic register, so the second silently skipped its deferred T gate and read an
  already-collapsed qubit.
- **Fix**: `build_rt_data` installs a fresh `invsparsity_history`, so reusing one
  runtime across shots no longer concatenates every shot's telemetry.
- `execute!(state)` runs an already-compiled `CompilerState`, splitting
  compilation from execution. `get_distribution` and `weight_std_graph` now
  compile once and run each shot off a copy instead of recompiling per shot.
- Clifford conjugation leaves commuting operations untouched rather than
  re-embedding them over the union of both supports. **Behaviour change**:
  Paulis in preprocessed circuits keep their own support instead of carrying
  identity padding (`Measurement(P"X", 1, [1])` where it used to be
  `Measurement(P"X_", 1, [1, 2])`). The operators are equivalent.
- Runtimes and `CompilerState` are parameterized on their memory and tableau
  types, making the per-measurement execution loop type stable.
- `random_test_circuit` accepts an `rng` keyword for reproducible circuits.
- Benchmarks cover `preprocess_circuit`, `run` and `get_distribution`, not just
  `traversal` with synthetic kernels.
- Basic datastructures for representing circuits through ADT with Moshi.jl
- Basic circuit traversal routines
- Plotting extension through Makie.jl