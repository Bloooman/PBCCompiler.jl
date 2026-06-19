using JLD2
using QuantumClifford
import Base: show

"""Sequence of abstract circuit operations."""
const Circuit = Vector{CircuitOp.Type}

"""
    parse_input(filepath::String) -> Circuit

Read an OpenQASM 2.0 or 3.0 file and return a `Circuit`.

The version is detected from the `OPENQASM X.Y;` header (defaults to 2.0).

Supported gates: `h`, `s`, `sdg`, `t`, `tdg`, `x`, `y`, `z`, `cx`, `ccx`, `measure`.
Each gate is translated to `CircuitOp` variants:
- `h`   → three `ExpQuatPiPauli`: Z, X, Z
- `s`   → `ExpQuatPiPauli(P"Z")`
- `sdg` → `ExpQuatPiPauli(-P"Z")`
- `t`   → `ExpEighPiPauli(P"Z")`
- `tdg` → `ExpEighPiPauli(-P"Z")`
- `x/y/z` → `ExpHalfPiPauli`
- `cx q[c],q[t]` → `PauliConditional(P"Z", [c], P"X", [t])`
- `ccx q[c1],q[c2],q[t]` → 15-op Toffoli decomposition (H, T, Tdg, CX sequence)

QASM 2.0 measure: `measure q[i] -> c[j]` → `Measurement(P"Z", j, [i])`
QASM 3.0 measure: `meas[j] = measure q[i]` → `Measurement(P"Z", j, [i])`

Header/declaration lines (`OPENQASM`, `include`, `barrier`, `creg`/`bit[N]`) are skipped.
If no `measure` statement is present, all qubits are measured at the end: qubit `q`
maps to classical bit `q`.
"""
function parse_input(filepath::String)::Circuit
    circuit = Circuit()
    n_qubits    = 0
    has_measure = false
    version     = 2

    open(filepath) do file
        for raw in eachline(file)
            line = strip(raw)
            isempty(line) && continue
            startswith(line, "//") && continue
            any(startswith(line, pfx) for pfx in ("include", "barrier")) && continue

            # detect OPENQASM version from header
            m = match(r"^OPENQASM\s+(\d+)\.", line)
            if m !== nothing
                version = parse(Int, m[1])
                continue
            end

            if version == 2
                startswith(line, "creg") && continue

                # qreg q[N];
                m = match(r"^qreg\s+\w+\[(\d+)\];$", line)
                if m !== nothing
                    n_qubits = parse(Int, m[1])
                    continue
                end

                # measure q[i] -> c[j];
                m = match(r"^measure\s+\w+\[(\d+)\]\s*->\s*\w+\[(\d+)\];$", line)
                if m !== nothing
                    q = parse(Int, m[1])+1
                    c = parse(Int, m[2])+1
                    push!(circuit, CircuitOp.Measurement(P"Z", c, [q]))
                    has_measure = true
                    continue
                end
            else
                startswith(line, "bit[") && continue

                # qubit[N] name;
                m = match(r"^qubit\[(\d+)\]\s+\w+;$", line)
                if m !== nothing
                    n_qubits = parse(Int, m[1])
                    continue
                end

                # meas[j] = measure q[i];
                m = match(r"^\w+\[(\d+)\]\s*=\s*measure\s+\w+\[(\d+)\];$", line)
                if m !== nothing
                    c = parse(Int, m[1])+1
                    q = parse(Int, m[2])+1
                    push!(circuit, CircuitOp.Measurement(P"Z", c, [q]))
                    has_measure = true
                    continue
                end
            end

            # cx q[ctrl],q[tgt];
            m = match(r"^[Cc][Xx]\s+\w+\[(\d+)\]\s*,\s*\w+\[(\d+)\];$", line)
            if m !== nothing
                ctrl = parse(Int, m[1])+1
                tgt  = parse(Int, m[2])+1
                push!(circuit, CircuitOp.PauliConditional(P"Z", [ctrl], P"X", [tgt]))
                continue
            end

            # ccx Toffoli;
            m = match(r"^[Cc][Cc][Xx]\s+\w+\[(\d+)\]\s*,\s*\w+\[(\d+)\],\s*\w+\[(\d+)\];$", line)
            if m !== nothing
                ctrl_1 = parse(Int, m[1])+1
                ctrl_2 = parse(Int, m[2])+1
                tgt    = parse(Int, m[3])+1
                append!(circuit, _toffoli_ops(ctrl_1, ctrl_2, tgt))
                continue
            end

            # single-qubit gate: <name> q[i];
            m = match(r"^(\w+)\s+\w+\[(\d+)\];$", line)
            if m !== nothing
                gate = m[1]
                q    = parse(Int, m[2])+1
                append!(circuit, _single_qubit_ops(gate, q))
                continue
            end
        end
    end

    if !has_measure
        for q in 1:n_qubits
            push!(circuit, CircuitOp.Measurement(P"Z", q, [q]))
        end
    end

    return circuit
end

"""
    _toffoli_ops(ctrl_1::Int, ctrl_2::Int, tgt::Int) -> Vector{CircuitOp.Type}

Return the 15-op decomposition of a Toffoli (CCX) gate into H, T, Tdg, and CX ops.

# Arguments
- `ctrl_1`: 1-indexed first control qubit
- `ctrl_2`: 1-indexed second control qubit
- `tgt`: 1-indexed target qubit

# Returns
15-element `Vector{CircuitOp.Type}` implementing the standard Toffoli decomposition.
"""
function _toffoli_ops(ctrl_1::Int, ctrl_2::Int, tgt::Int)::Vector{CircuitOp.Type}
    cx(c, t) = CircuitOp.PauliConditional(P"Z", [c], P"X", [t])
    ops = CircuitOp.Type[]
    append!(ops, _single_qubit_ops("h",   tgt))
    push!(ops,   cx(ctrl_2, tgt))
    append!(ops, _single_qubit_ops("tdg", tgt))
    push!(ops,   cx(ctrl_1, tgt))
    append!(ops, _single_qubit_ops("t",   tgt))
    push!(ops,   cx(ctrl_2, tgt))
    append!(ops, _single_qubit_ops("tdg", tgt))
    push!(ops,   cx(ctrl_1, tgt))
    append!(ops, _single_qubit_ops("t",   ctrl_2))
    append!(ops, _single_qubit_ops("t",   tgt))
    append!(ops, _single_qubit_ops("h",   tgt))
    push!(ops,   cx(ctrl_1, ctrl_2))
    append!(ops, _single_qubit_ops("t",   ctrl_1))
    append!(ops, _single_qubit_ops("tdg", ctrl_2))
    push!(ops,   cx(ctrl_1, ctrl_2))
    return ops
end

# Internal: translate a single-qubit gate name to its CircuitOps.
function _single_qubit_ops(gate::AbstractString, q::Int)::Vector{CircuitOp.Type}
    if gate == "h"
        return CircuitOp.Type[
            CircuitOp.ExpQuatPiPauli(P"Z", [q]),
            CircuitOp.ExpQuatPiPauli(P"X", [q]),
            CircuitOp.ExpQuatPiPauli(P"Z", [q]),
        ]
    elseif gate == "s"
        return [CircuitOp.ExpQuatPiPauli(P"Z",  [q])]
    elseif gate == "sdg"
        return [CircuitOp.ExpQuatPiPauli(-P"Z", [q])]
    elseif gate == "t"
        return [CircuitOp.ExpEighPiPauli(P"Z",  [q])]
    elseif gate == "tdg"
        return [CircuitOp.ExpEighPiPauli(-P"Z", [q])]
    elseif gate == "x"
        return [CircuitOp.ExpHalfPiPauli(P"X",  [q])]
    elseif gate == "y"
        return [CircuitOp.ExpHalfPiPauli(P"Y",  [q])]
    elseif gate == "z"
        return [CircuitOp.ExpHalfPiPauli(P"Z",  [q])]
    else
        error("Unsupported gate: $gate")
    end
end


"""
    save(result::CompilerState, filepath::String)

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
        stabilizer_group     = result.stabilizer_group,
        classical_register   = result.classical_register,
        instruction_pointer  = result.instruction_pointer,
    )
end

"""
    save(r::CompilationResult, filepath::String)

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
    load(filepath::String) -> CompilationResult

Load a `CompilationResult` from a `.jld2` file previously written by `save`.

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
operations. Supports Clifford+T gates: `x`, `y`, `z`, `h`, `cx`, `s`, `sdg`,
`t`, `tdg`. Qubit indices are 1-based. T† is decomposed as T followed by S†.
"""
function parse_QuantumClifford(filepath::String)
    circuit = QuantumClifford.AbstractOperation[]
    qubit_map = Dict{String,Int}()
    qubit_offset = 0

    for raw_line in eachline(filepath)
        line = strip(raw_line)

        # Strip inline comments
        ci = findfirst("//", line)
        !isnothing(ci) && (line = strip(line[1:ci.start-1]))
        isempty(line) && continue

        # Skip directives that carry no gate information
        (startswith(line, "OPENQASM") || startswith(line, "include") ||
         startswith(line, "creg")     || startswith(line, "measure") ||
         startswith(line, "reset")    || startswith(line, "barrier") ||
         startswith(line, "gate")     || startswith(line, "opaque")) && continue

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

        # Drop any parenthesised parameter list (e.g. U(theta,phi,lambda))
        args_str = replace(args_str, r"\([^)]*\)" => "")
        args_str = strip(args_str)

        qargs = [strip(q) for q in split(args_str, ",") if !isempty(strip(q))]
        isempty(qargs) && continue

        qs = [get(qubit_map, q, nothing) for q in qargs]
        any(isnothing, qs) && continue  # skip unresolvable qubit references

        if gate == "x" && length(qs) == 1
            push!(circuit, sX(qs[1]))
        elseif gate == "y" && length(qs) == 1
            push!(circuit, sY(qs[1]))
        elseif gate == "z" && length(qs) == 1
            push!(circuit, sZ(qs[1]))
        elseif gate == "h" && length(qs) == 1
            push!(circuit, sHadamard(qs[1]))
        elseif (gate == "cx" || gate == "cnot") && length(qs) == 2
            push!(circuit, sCNOT(qs[1], qs[2]))
        elseif gate == "s" && length(qs) == 1
            push!(circuit, sPhase(qs[1]))
        elseif gate == "sdg" && length(qs) == 1
            push!(circuit, sInvPhase(qs[1]))
        elseif gate == "t" && length(qs) == 1
            push!(circuit, sT(qs[1]))
        elseif gate == "tdg" && length(qs) == 1
            # T† = S†·T, so apply T first then S†
            push!(circuit, sT(qs[1]))
            push!(circuit, sInvPhase(qs[1]))
        end
    end

    return circuit
end

"""
    get_T_count(circuit::Vector{QuantumClifford.AbstractOperation}) -> Int

Return the number of `sT` gates in `circuit`.
"""
function get_T_count(circuit::Vector{QuantumClifford.AbstractOperation})
    return count(op -> op isa sT, circuit)
end
