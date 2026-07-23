using JLD2
using QuantumClifford
import Base: show

"""
Parsed gate definition from a `.inc` file or an inline `gate` block.

Fields:
- `name`: lowercase gate name
- `params`: formal angle parameter names (empty for non-parameterized gates)
- `qubits`: formal qubit argument names
- `body`: gate body split into individual statements (semicolons removed)
"""
struct GateDef
    name::String
    params::Vector{String}
    qubits::Vector{String}
    body::Vector{String}
end

"""
    parse_input(filepath::String) -> Circuit

Read an OpenQASM 2.0 or 3.0 file and return a `Circuit`.

The version is detected from the `OPENQASM X.Y;` header (defaults to 2.0).

`include "filename.inc";` directives are followed and gate definitions are loaded
from the referenced file. Gate definitions in the circuit file itself are also parsed.
Loaded gates are expanded recursively into the primitive Clifford+T gate set.

Primitive gates (mapped directly to `CircuitOp` variants):
- `h`   → three `ExpQuatPiPauli`: Z, X, Z
- `s`   → `ExpQuatPiPauli(P"Z")`
- `sdg` → `ExpQuatPiPauli(-P"Z")`
- `t`   → `ExpEighPiPauli(P"Z")`
- `tdg` → `ExpEighPiPauli(-P"Z")`
- `x/y/z` → `ExpHalfPiPauli`
- `cx q[c],q[t]` → `PauliConditional(P"Z", [c], P"X", [t])`
- `id` → no-op

Composite gates (`ccx`, `cz`, `cy`, `ch`, etc.) are expanded via their definitions
from a loaded `.inc` file.

Parameterized gate calls (`rx(θ)`, `u3(θ,φ,λ)`, etc.) raise an error — they are
not representable in the Clifford+T gate set.

Registers are laid out consecutively in declaration order: with `qreg a[2]; qreg b[1];`,
`a[0]`→qubit 1, `a[1]`→qubit 2, `b[0]`→qubit 3, and likewise for classical registers.

QASM 2.0 measure: `measure q[i] -> c[j]` → `Measurement(P"Z", bit(c[j]), [qubit(q[i])])`
QASM 3.0 measure: `meas[j] = measure q[i]` → `Measurement(P"Z", bit(meas[j]), [qubit(q[i])])`

References to undeclared registers raise an error.

If no `measure` statement is present, all qubits are measured at the end: qubit `q`
maps to classical bit `q`.
"""
function parse_input(filepath::String)::Circuit
    content = read(filepath, String)
    dir     = dirname(abspath(filepath))

    # gate_defs: includes first, then inline definitions (later wins)
    gate_defs = Dict{String,GateDef}()
    for m in eachmatch(r"include\s+\"([^\"]+)\"\s*;", content)
        inc_path = joinpath(dir, m[1])
        isfile(inc_path) && _load_gate_defs!(gate_defs, read(inc_path, String))
    end
    _load_gate_defs!(gate_defs, content)

    circuit     = Circuit()
    # Registers are laid out consecutively in declaration order, so files with
    # several qreg/creg declarations map to distinct global indices instead of
    # silently collapsing onto the same qubits/bits
    qubit_map    = Dict{String,Int}()   # "q[0]" -> 1-based global qubit index
    bit_map      = Dict{String,Int}()   # "c[0]" -> 1-based global bit index
    qubit_offset = 0
    bit_offset   = 0
    has_measure = false
    in_gate_def = false

    for raw in split(content, '\n')
        line = strip(replace(raw, r"\s*//.*$" => ""))
        isempty(line) && continue

        if in_gate_def
            contains(line, "}") && (in_gate_def = false)
            continue
        end
        if startswith(line, "gate ")
            in_gate_def = !contains(line, "}")
            continue
        end

        any(startswith(line, p) for p in ("OPENQASM", "barrier", "include")) && continue

        # quantum registers: v2 "qreg q[N];"  or  v3 "qubit[N] name;"
        m = match(r"^qreg\s+(\w+)\[(\d+)\]", line)
        m === nothing && (m = match(r"^qubit\[(\d+)\]\s+(\w+)", line);
                          m === nothing || (m = (m[2], m[1])))
        if m !== nothing
            qubit_offset = _declare_register!(qubit_map, String(m[1]), parse(Int, m[2]), qubit_offset)
            continue
        end

        # classical registers: v2 "creg c[N];"  or  v3 "bit[N] name;"
        m = match(r"^creg\s+(\w+)\[(\d+)\]", line)
        m === nothing && (m = match(r"^bit\[(\d+)\]\s+(\w+)", line);
                          m === nothing || (m = (m[2], m[1])))
        if m !== nothing
            bit_offset = _declare_register!(bit_map, String(m[1]), parse(Int, m[2]), bit_offset)
            continue
        end

        # measure: v2 "measure q[i] -> c[j];"
        m = match(r"^measure\s+(\w+\[\d+\])\s*->\s*(\w+\[\d+\]);$", line)
        if m !== nothing
            push!(circuit, CircuitOp.Measurement(P"Z", _resolve(bit_map, m[2], "classical bit", line),
                                                 [_resolve(qubit_map, m[1], "qubit", line)]))
            has_measure = true; continue
        end
        # measure: v3 "meas[j] = measure q[i];"
        m = match(r"^(\w+\[\d+\])\s*=\s*measure\s+(\w+\[\d+\]);$", line)
        if m !== nothing
            push!(circuit, CircuitOp.Measurement(P"Z", _resolve(bit_map, m[1], "classical bit", line),
                                                 [_resolve(qubit_map, m[2], "qubit", line)]))
            has_measure = true; continue
        end

        # reject parameterized gate calls
        m = match(r"^(\w+)\s*\([^)]*\)\s+\w+\[", line)
        m !== nothing && error("Parameterized gate '$(m[1])' is not supported — only Clifford+T gates are representable")

        # gate application: name q[i] or name q[i],q[j] or name q[i],q[j],q[k]
        m = match(r"^(\w+)\s+(.+);$", line)
        if m !== nothing
            refs = [String(r.match) for r in eachmatch(r"\w+\[\d+\]", m[2])]
            if !isempty(refs)
                qubits = [_resolve(qubit_map, ref, "qubit", line) for ref in refs]
                _apply_gate!(circuit, String(m[1]), qubits, gate_defs)
                continue
            end
        end
    end

    if !has_measure
        for q in 1:qubit_offset
            push!(circuit, CircuitOp.Measurement(P"Z", q, [q]))
        end
    end
    return circuit
end

"""
Map every element `name[i]` of a newly declared register onto consecutive
1-based global indices starting after `offset`; return the new offset.
"""
function _declare_register!(map::Dict{String,Int}, name::String, size::Int, offset::Int)
    for i in 0:size-1
        map["$name[$i]"] = offset + i + 1
    end
    return offset + size
end

"""Resolve a `name[i]` reference to its global index, erroring if the register was never declared."""
function _resolve(map::Dict{String,Int}, ref::AbstractString, kind::String, line::AbstractString)
    haskey(map, ref) ||
        error("Unknown $kind reference '$ref' in line '$line' — register not declared")
    return map[ref]
end

"""
    _load_gate_defs!(defs, content)

Parse all `gate name(params) qubits { body }` blocks in `content` (QASM or .inc
text) and store them in `defs`. Existing entries are overwritten, so later includes
and inline definitions take precedence.

Gate names are normalised to lowercase. Comments are stripped before parsing.
"""
function _load_gate_defs!(defs::Dict{String,GateDef}, content::AbstractString)
    cleaned = replace(content, r"//[^\n]*" => "")
    for m in eachmatch(r"gate\s+(\w+)\s*(?:\(([^)]*)\))?\s+([^{]+)\{([^}]*)\}", cleaned)
        name   = lowercase(strip(m[1]))
        params = m[2] !== nothing ? strip.(split(m[2], ',')) : String[]
        qubits = strip.(split(strip(m[3]), ','))
        body   = filter(!isempty, strip.(split(m[4], ';')))
        defs[name] = GateDef(name, params, qubits, body)
    end
end

"""
    _apply_gate!(circuit, name, qubit_indices, gate_defs)

Expand gate `name` onto `qubit_indices` (1-indexed) into `circuit`.

Primitive gates (`h`, `s`, `sdg`, `t`, `tdg`, `x`, `y`, `z`, `cx`, `id`) are
mapped directly to `CircuitOp` variants. All other gates are looked up in
`gate_defs` and expanded recursively. Errors if a gate is unknown or if its body
reaches a parameterized primitive (e.g., `U`, `u1`) that is outside the Clifford+T
gate set.
"""
function _apply_gate!(circuit::Circuit, name::String, qubit_indices::Vector{Int},
                      gate_defs::Dict{String,GateDef})
    lname = lowercase(name)

    if lname in ("h", "s", "sdg", "t", "tdg", "x", "y", "z")
        length(qubit_indices) == 1 || error("Gate '$name' takes 1 qubit, got $(length(qubit_indices))")
        append!(circuit, _single_qubit_ops(lname, qubit_indices[1]))
        return
    elseif lname == "id"
        return
    elseif lname == "cx"
        length(qubit_indices) == 2 || error("Gate 'cx' takes 2 qubits, got $(length(qubit_indices))")
        push!(circuit, CircuitOp.PauliConditional(P"Z", [qubit_indices[1]], P"X", [qubit_indices[2]]))
        return
    end

    haskey(gate_defs, lname) ||
        error("Unknown gate '$name' — not a Clifford+T primitive and not defined in any loaded .inc file")

    def = gate_defs[lname]
    length(qubit_indices) == length(def.qubits) ||
        error("Gate '$name' expects $(length(def.qubits)) qubits, got $(length(qubit_indices))")
    qubit_map = Dict(zip(def.qubits, qubit_indices))

    for stmt in def.body
        stmt = strip(stmt)
        isempty(stmt) && continue
        if match(r"^\w+\s*\(", stmt) !== nothing
            error("Gate '$(def.name)' expands to parameterized primitive '$(split(stmt,'(')[1])' — not representable in the Clifford+T gate set")
        end
        m = match(r"^(\w+)\s+(.+)$", stmt)
        m !== nothing || error("Cannot parse body statement in gate '$(def.name)': '$stmt'")
        _apply_gate!(circuit, String(m[1]),
                     [qubit_map[strip(String(s))] for s in split(m[2], ',')],
                     gate_defs)
    end
end

"""
    _single_qubit_ops(gate, q) -> Vector{CircuitOp.Type}

Translate a primitive single-qubit gate name to its `CircuitOp` sequence.
Supported: `h`, `s`, `sdg`, `t`, `tdg`, `x`, `y`, `z`.
"""
function _single_qubit_ops(gate::AbstractString, q::Int)::Vector{CircuitOp.Type}
    Q(p) = CircuitOp.ExpQuatPiPauli(p, [q])
    E(p) = CircuitOp.ExpEighPiPauli(p, [q])
    H(p) = CircuitOp.ExpHalfPiPauli(p, [q])
    gate == "h"   && return [Q(P"Z"), Q(P"X"), Q(P"Z")]
    gate == "s"   && return [Q(P"Z")]
    gate == "sdg" && return [Q(-P"Z")]
    gate == "t"   && return [E(P"Z")]
    gate == "tdg" && return [E(-P"Z")]
    gate == "x"   && return [H(P"X")]
    gate == "y"   && return [H(P"Y")]
    gate == "z"   && return [H(P"Z")]
    error("Unsupported primitive gate: $gate")
end

"""
    save_result(result::CompilerState, filepath::String)

Save transient compiler state fields to a `.jld2` file for debugging.

Saves `measurement_results`, `stabilizer_group`, `classical_register`, and
`instruction_pointer`. Omits `circuit` (fixed input) and `runtime` (ADT).

# Arguments
- `result`: compiler state to serialize
- `filepath`: destination path (should end in `.jld2`)

# Returns
Nothing.
"""
function save_result(result::CompilerState, filepath::String)
    mkpath(dirname(filepath))
    JLD2.jldsave(filepath;
        measurement_results  = result.measurement_results,
        # keep the on-disk representation a plain Stabilizer regardless of the
        # in-memory tableau type
        stabilizer_group     = copy(stabilizerview(result.stabilizer_group)),
        classical_register   = result.classical_register,
        instruction_pointer  = result.instruction_pointer,
    )
end

"""
    save_result(r::CompilationResult, filepath::String)

Save all fields of a `CompilationResult` to a `.jld2` file at `filepath`.

# Arguments
- `r`: compilation result to serialize
- `filepath`: destination path (should end in `.jld2`)

# Returns
Nothing.
"""
function save_result(r::CompilationResult, filepath::String)
    mkpath(dirname(filepath))
    JLD2.jldsave(filepath;
        measurement_results = r.measurement_results,
        QPU_workload        = r.QPU_workload,
        stabilizer_group    = r.stabilizer_group,
        QPUDuration         = r.QPUDuration,
    )
end

"""
    load_result(filepath::String) -> CompilationResult

Load a `CompilationResult` from a `.jld2` file previously written by `save_result`.

# Arguments
- `filepath`: path to a `.jld2` file

# Returns
A `CompilationResult` reconstructed from the saved fields.
"""
function load_result(filepath::String)::CompilationResult
    data = JLD2.load(filepath)
    return CompilationResult(
        data["measurement_results"],
        data["QPU_workload"],
        data["stabilizer_group"],
        data["QPUDuration"],
    )
end
##
"""
    parse_QuantumClifford(filepath::String) -> Vector{AbstractOperation}

Read an OpenQASM 2.0 file and return a circuit as a vector of QuantumClifford
operations.

Primitive Clifford+T gates (`x`, `y`, `z`, `h`, `cx`/`cnot`, `s`, `sdg`, `t`,
`tdg`) are mapped directly; T† is decomposed as T followed by S†. Composite
gates (e.g. `ccx`) are expanded recursively via `gate` definitions loaded from
`include`d `.inc` files and from inline `gate` blocks. Unknown gates and
parameterized gate calls raise an error. Qubit indices are 1-based; multiple
`qreg` declarations are laid out consecutively.
"""
function parse_QuantumClifford(filepath::String)
    content = read(filepath, String)
    dir     = dirname(abspath(filepath))

    # gate definitions: includes first, then inline definitions (later wins)
    gate_defs = Dict{String,GateDef}()
    for m in eachmatch(r"include\s+\"([^\"]+)\"\s*;", content)
        inc_path = joinpath(dir, m[1])
        isfile(inc_path) && _load_gate_defs!(gate_defs, read(inc_path, String))
    end
    _load_gate_defs!(gate_defs, content)

    circuit = QuantumClifford.AbstractOperation[]
    qubit_map = Dict{String,Int}()
    qubit_offset = 0
    in_gate_def = false

    for raw_line in split(content, '\n')
        line = strip(replace(raw_line, r"\s*//.*$" => ""))
        isempty(line) && continue

        # Skip gate definition blocks; their bodies were consumed by _load_gate_defs!
        if in_gate_def
            contains(line, "}") && (in_gate_def = false)
            continue
        end
        if startswith(line, "gate ")
            in_gate_def = !contains(line, "}")
            continue
        end

        # Skip directives that carry no gate information
        (startswith(line, "OPENQASM") || startswith(line, "include") ||
         startswith(line, "creg")     || startswith(line, "measure") ||
         startswith(line, "reset")    || startswith(line, "barrier") ||
         startswith(line, "opaque")) && continue

        # qreg declaration: build name[index] -> 1-based qubit number map
        m = match(r"^qreg\s+(\w+)\[(\d+)\]\s*;", line)
        if !isnothing(m)
            reg, sz = m.captures[1], parse(Int, m.captures[2])
            for i in 0:sz-1
                qubit_map["$(reg)[$(i)]"] = qubit_offset + i + 1
            end
            qubit_offset += sz
            continue
        end

        # Gate application — strip trailing semicolon then parse
        line = endswith(line, ";") ? line[1:end-1] : line
        m = match(r"^(\w+)\s*(.*)$", line)
        isnothing(m) && continue

        gate = lowercase(m.captures[1])
        args_str = strip(m.captures[2])
        contains(args_str, "(") &&
            error("Parameterized gate '$gate' is not supported — only Clifford+T gates are representable")

        qargs = [strip(q) for q in split(args_str, ",") if !isempty(strip(q))]
        isempty(qargs) && continue

        qs = [begin
                  haskey(qubit_map, q) || error("Unknown qubit reference '$q' in line '$line'")
                  qubit_map[q]
              end for q in qargs]
        _apply_gate_qc!(circuit, gate, qs, gate_defs)
    end

    return circuit
end

"""
    _apply_gate_qc!(circuit, name, qs, gate_defs)

Expand gate `name` onto 1-based qubit indices `qs` into a vector of
QuantumClifford operations. Primitive Clifford+T gates are mapped directly;
all other gates are looked up in `gate_defs` and expanded recursively.
Errors on unknown gates and parameterized primitives.
"""
function _apply_gate_qc!(circuit::Vector{QuantumClifford.AbstractOperation}, name::AbstractString,
                         qs::Vector{Int}, gate_defs::Dict{String,GateDef})
    gate = lowercase(name)
    single = Dict("x" => sX, "y" => sY, "z" => sZ, "h" => sHadamard,
                  "s" => sPhase, "sdg" => sInvPhase, "t" => sT)
    if gate == "id"
        return
    elseif haskey(single, gate)
        length(qs) == 1 || error("Gate '$name' takes 1 qubit, got $(length(qs))")
        push!(circuit, single[gate](qs[1]))
        return
    elseif gate == "tdg"
        length(qs) == 1 || error("Gate 'tdg' takes 1 qubit, got $(length(qs))")
        # T† = S†·T, so apply T first then S†
        push!(circuit, sT(qs[1]))
        push!(circuit, sInvPhase(qs[1]))
        return
    elseif gate == "cx" || gate == "cnot"
        length(qs) == 2 || error("Gate '$name' takes 2 qubits, got $(length(qs))")
        push!(circuit, sCNOT(qs[1], qs[2]))
        return
    end

    haskey(gate_defs, gate) ||
        error("Unknown gate '$name' — not a Clifford+T primitive and not defined in any loaded .inc file")

    def = gate_defs[gate]
    length(qs) == length(def.qubits) ||
        error("Gate '$name' expects $(length(def.qubits)) qubits, got $(length(qs))")
    qubit_map = Dict(zip(def.qubits, qs))

    for stmt in def.body
        stmt = strip(stmt)
        isempty(stmt) && continue
        if match(r"^\w+\s*\(", stmt) !== nothing
            error("Gate '$(def.name)' expands to parameterized primitive '$(split(stmt,'(')[1])' — not representable in the Clifford+T gate set")
        end
        m = match(r"^(\w+)\s+(.+)$", stmt)
        m !== nothing || error("Cannot parse body statement in gate '$(def.name)': '$stmt'")
        _apply_gate_qc!(circuit, String(m[1]),
                        [qubit_map[strip(String(s))] for s in split(m[2], ',')],
                        gate_defs)
    end
end

"""
    get_T_count(circuit::Vector{QuantumClifford.AbstractOperation}) -> Int

Return the number of `sT` gates in `circuit`.
"""
function get_T_count(circuit::Vector{QuantumClifford.AbstractOperation})
    return count(op -> op isa sT, circuit)
end
