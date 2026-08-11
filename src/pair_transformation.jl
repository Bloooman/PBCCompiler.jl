"""
    paulis(op::CircuitOp.Type) -> PauliOperator

Return Pauli string that defines the collective eigen-axis of the operation; a Pauli Product Rotation rotates the multi-qubit state around this axis by an angle phi,
while a Pauli Product Measurement projects the system directly into its collective directional eigenstates

# Examples
```jldoctest
julia> op = PBCCompiler.Pauli(P"XY", [1, 2]);

julia> PBCCompiler.paulis(op)
+ XY
```
```jldoctest
julia> op = PBCCompiler.PauliConditional(P"X", [1], P"Z", [3]);

julia> PBCCompiler.paulis(op)
ERROR: pauli called on PauliConditional — decompose first
```
"""
function paulis(op::CircuitOp.Type)
    pauli = @match op begin
        CircuitOp.Pauli(pauli, qubits) => pauli
        CircuitOp.Measurement(pauli, bit, qubits) => pauli
        CircuitOp.ExpHalfPiPauli(pauli, qubits) => pauli
        CircuitOp.ExpQuatPiPauli(pauli, qubits) => pauli
        CircuitOp.ExpEighPiPauli(pauli, qubits) => pauli
        CircuitOp.PauliConditional(cp, cq, tp, tq) => error("pauli called on PauliConditional — decompose first")
        CircuitOp.BitConditional(inner_op, bit) => paulis(inner_op)
    end
    return pauli
end

"""
    complete_paulis(op1::CircuitOp.Type, op2::CircuitOp.Type) -> (PauliOperator, PauliOperator)

This helper function ensures that both operators are represented over the
union of their affected qubits. It reorders strings to a canonical qubit
ordering and pads missing sites with Identity operators ('_') to ensure
equal string length.

# Examples
```jldoctest
julia> op1 = PBCCompiler.ExpQuatPiPauli(P"XY", [1, 3]);

julia> op2 = PBCCompiler.ExpQuatPiPauli(P"ZXY",[3, 1, 2]);

julia> PBCCompiler.complete_paulis(op1,op2)
(+ X_Y, + XYZ)
```
```jldoctest
julia> op1 = PBCCompiler.ExpQuatPiPauli(P"XY", [1, 3]);

julia> op2 = PBCCompiler.ExpHalfPiPauli(P"Z", [5]);

julia> PBCCompiler.complete_paulis(op1,op2)
(+ X_Y__, + ____Z)
```
"""
function complete_paulis(op1::CircuitOp.Type, op2::CircuitOp.Type)
    pu1 = paulis(op1)
    pu2 = paulis(op2)
    paulilen = max(max_affected_qubit(op1), max_affected_qubit(op2))
    pauli1 = embed(paulilen, op1.qubits, pu1)
    pauli2 = embed(paulilen, op2.qubits, pu2)
    return (pauli1, pauli2)
end

"""
    check_commutation(op1::CircuitOp.Type, op2::CircuitOp.Type) -> Union{Int8, Nothing}

Return 0x00 if the two Pauli Product Rotations/Measurements commute, and return 0x01 if they anticommute.
For inputs containing PauliConditional or BitConditional, the function returns nothing.

# Examples
```jldoctest
julia> op1 = PBCCompiler.ExpQuatPiPauli(P"XY", [1, 3]);

julia> op2 = PBCCompiler.ExpQuatPiPauli(P"ZXY",[3, 1, 2]);

julia> PBCCompiler.check_commutation(op1, op2)
0x01
```
```jldoctest
julia> op1 = PBCCompiler.ExpQuatPiPauli(P"XY", [1, 3]);

julia> CNOT = PBCCompiler.PauliConditional(P"Z", [1], P"X", [2]);

julia> PBCCompiler.check_commutation(op1, CNOT)
```
"""
function check_commutation(op1::CircuitOp.Type, op2::CircuitOp.Type)
    @match (op1, op2) begin
        (op,CircuitOp.BitConditional(inner_op, bit)) || (CircuitOp.BitConditional(inner_op, bit), op) => begin
            return nothing
        end
        (op,CircuitOp.PauliConditional()) || (CircuitOp.PauliConditional(), op) => begin
            return nothing
        end
        _ => begin
            (pauli1,pauli2) = complete_paulis(op1, op2)
            commutativity = comm(pauli1,pauli2)
            return commutativity
        end
    end
end

"""
Largest `|q1| * |q2|` for which [`paulis_commute`](@ref) walks the supports
directly instead of falling back to the bit-parallel `comm`. See the comment at
the branch for how this was chosen.
"""
const PAULI_WALK_MAX_WORK = 256

"""
Test whether two Pauli-carrying operations commute, without allocating.

[`check_commutation`](@ref) answers the same question by embedding both Paulis
to a common width first, which costs two `PauliOperator` allocations per call.
This walks the two support lists instead: a qubit outside the shared support
carries an identity on one side and cannot contribute, and the pair
anticommutes exactly when an odd number of shared qubits carry anticommuting
single-qubit Paulis.

Both operands must be plain Pauli-carrying variants (with `pauli` and `qubits`
fields); conditionals are rejected by the callers before this is reached.
"""
function paulis_commute(op1::CircuitOp.Type, op2::CircuitOp.Type)
    q1 = op1.qubits
    q2 = op2.qubits
    # The support walk is O(|q1| * |q2|) scalar work, while `comm` on embedded
    # Paulis is bit-parallel (one word per 64 qubits) at the cost of two
    # allocations. Operations reaching here are usually narrow, where the walk
    # wins by a wide margin, but dense multi-qubit rotations on a large register
    # flip that -- measured crossover is around a product of a few hundred.
    if length(q1) * length(q2) > PAULI_WALK_MAX_WORK
        (embedded1, embedded2) = complete_paulis(op1, op2)
        return comm(embedded1, embedded2) == 0x00
    end
    p1 = paulis(op1)
    p2 = paulis(op2)
    anticommuting = false
    for i in eachindex(q1)
        j = findfirst(isequal(q1[i]), q2)
        j === nothing && continue
        (x1, z1) = p1[i]
        (x2, z2) = p2[j]
        anticommuting ⊻= (x1 & z2) ⊻ (z1 & x2)
    end
    return !anticommuting
end

"""
Conjugate op2's Pauli by the Clifford rotation op1 via `conjugated_by_clifford`
and rebuild op2 with `rebuild(pauli, qubits)`. A commuting pair swaps unchanged
(returns `(op2, op1)` directly) rather than rebuilding op2 at the union width;
otherwise returns `(rebuild(conjugated...), op1)`.
"""
function _conjugate_via(op1::CircuitOp.Type, op2::CircuitOp.Type, rebuild)
    conjugated = conjugated_by_clifford(op1,op2)
    conjugated === nothing && return (op2, op1)
    return (rebuild(conjugated[1], conjugated[2]), op1)
end

"""
    conjugate_noncliff(op1::CircuitOp.Type, op2::CircuitOp.Type) -> Union{Tuple{CircuitOp.Type, CircuitOp.Type}, Nothing}

Move a non-Clifford CircuitOp op2 pass a Clifford CircuitOp op1 and update op2 by conjugating its pauli string by op1's pauli string.
Will throw an error if op2 is not a ExpEighPiPauli CircuitOp. Return nothing if op1 is not a ExpHalfPiPauli or a ExpQuatPiPauli.

# Examples
```jldoctest
julia> op1 = PBCCompiler.ExpQuatPiPauli(P"XY", [1, 3]);

julia> op2 = PBCCompiler.ExpEighPiPauli(P"ZXY",[3, 1, 2]);

julia> PBCCompiler.conjugate_noncliff(op1, op2)
(CircuitOp.ExpEighPiPauli(pauli=- _YX, qubits=[1, 2, 3]), CircuitOp.ExpQuatPiPauli(pauli=+ XY, qubits=[1, 3]))
```
```jldoctest
julia> op1 = PBCCompiler.ExpQuatPiPauli(P"XY", [1, 3]);

julia> M_Z = PBCCompiler.Measurement(P"Z", 1, [2]);

julia> PBCCompiler.conjugate_noncliff(op1, M_Z)
```
"""
function conjugate_noncliff(op1::CircuitOp.Type, op2::CircuitOp.Type)
    @match (op1, op2) begin
        (CircuitOp.ExpHalfPiPauli(), CircuitOp.ExpEighPiPauli()) ||
        (CircuitOp.ExpQuatPiPauli(), CircuitOp.ExpEighPiPauli()) => begin
            _conjugate_via(op1, op2, CircuitOp.ExpEighPiPauli)
        end
        _=> nothing
    end
end
##

"""
    conjugate_measurement(op1::CircuitOp.Type, op2::CircuitOp.Type) -> Union{Tuple{CircuitOp.Type, CircuitOp.Type}, Nothing}

Move a Measurement CircuitOp op2 past a Clifford CircuitOp op1 and update op2 by conjugating its Pauli string by op1's Pauli string.
Will throw an error if op2 is not a Measurement CircuitOp. Return nothing if op1 is not a ExpHalfPiPauli or a ExpQuatPiPauli.

A measurement that commutes with `op1` crosses it unchanged and keeps its own
support; only an anticommuting one is rebuilt over the union of both supports.

# Examples
```jldoctest
julia> op1 = PBCCompiler.ExpQuatPiPauli(P"XY", [1, 3]);

julia> M_Z=PBCCompiler.Measurement(P"Z", 1, [2]);

julia> PBCCompiler.conjugate_measurement(op1, M_Z)
(CircuitOp.Measurement(pauli=+ Z, bit=1, qubits=[2]), CircuitOp.ExpQuatPiPauli(pauli=+ XY, qubits=[1, 3]))
```
```jldoctest
julia> op1 = PBCCompiler.ExpQuatPiPauli(P"X", [1]);

julia> M_Z=PBCCompiler.Measurement(P"Z", 1, [1]);

julia> PBCCompiler.conjugate_measurement(op1, M_Z)
(CircuitOp.Measurement(pauli=+ Y, bit=1, qubits=[1]), CircuitOp.ExpQuatPiPauli(pauli=+ X, qubits=[1]))
```
```jldoctest
julia> op1 = PBCCompiler.ExpQuatPiPauli(P"XY", [1, 3]);

julia> op2 = PBCCompiler.ExpEighPiPauli(P"ZXY",[3, 1, 2]);

julia> PBCCompiler.conjugate_measurement(op1, op2)
```
"""
function conjugate_measurement(op1::CircuitOp.Type, op2::CircuitOp.Type)
    @match (op1, op2) begin
        (CircuitOp.ExpHalfPiPauli(), CircuitOp.Measurement()) ||
        (CircuitOp.ExpQuatPiPauli(), CircuitOp.Measurement()) => begin
            # See `conjugate_noncliff`: a commuting measurement crosses unchanged
            _conjugate_via(op1, op2, (p, q) -> CircuitOp.Measurement(p, op2.bit, q))
        end
        _=> nothing
    end
end
##
"""
Conjugate op2's Pauli by the Clifford rotation op1 (op1 must be an
ExpHalfPiPauli or ExpQuatPiPauli).

Returns `nothing` when the two Paulis commute -- op2 is then unchanged and the
caller should reuse it. Otherwise returns the conjugated Pauli together with its
qubit list, expanded to `1:width` over the union of both ops' supports:
conjugation by exp(-iπ/2·P₁) gives -P₂ and by exp(-iπ/4·P₁) gives i·P₁·P₂.
"""
function conjugated_by_clifford(op1::CircuitOp.Type, op2::CircuitOp.Type)
    # Commuting pairs return `nothing`: op2 is unchanged, so the caller reuses it
    # as-is. Widening it to the union width -- which is all the old code did on
    # this branch -- would allocate twice and would stick, making every later
    # embed/comm on the op run over the full register instead of its support.
    #
    # The two branches differ only in how commutation is decided, and each
    # embeds at most once: on the wide path `complete_paulis` serves both the
    # test and the conjugation, so testing never costs an extra embedding.
    local pauli1, pauli2
    if length(op1.qubits) * length(op2.qubits) > PAULI_WALK_MAX_WORK
        (pauli1, pauli2) = complete_paulis(op1, op2)
        comm(pauli1, pauli2) == 0x00 && return nothing
    else
        paulis_commute(op1, op2) && return nothing
        (pauli1, pauli2) = complete_paulis(op1, op2)
    end
    # complete_paulis already embedded both Paulis to this width
    new_q = collect(1:length(pauli1))
    new_p = isa_variant(op1, CircuitOp.ExpHalfPiPauli) ? -pauli2 : 1im*pauli1*pauli2
    return (new_p, new_q)
end
##
"""
Merge two same-variant, same-axis rotations into one `mergedtype` rotation of
twice the angle, or cancel them to `()` if they're inverses. Shared by the
`ExpEighPiPauli`/`ExpEighPiPauli` and `ExpQuatPiPauli`/`ExpQuatPiPauli` arms of
`merge_rotations`, which differ only in which "double-angle" type the merge
produces.
"""
function _merge_same_axis(op1::CircuitOp.Type, op2::CircuitOp.Type, mergedtype)
    (p1,p2) = complete_paulis(op1,op2)
    # The qubit list is only needed on the merge branch; building it up
    # front pays for every non-matching pair, which is nearly all of them
    p1.xz == p2.xz || return nothing
    if xor(p1.phase[1], p2.phase[1]) == 0x02
        # exp(-i*theta*P) * exp(-i*theta*(-P)) = identity
        return ()
    elseif op1.pauli.phase == op2.pauli.phase
        return mergedtype(p1,collect(1:length(p1)))
    else
        return nothing
    end
end

"""
    merge_rotations(op1::CircuitOp.Type, op2::CircuitOp.Type)

Combine an adjacent pair of Pauli Product Rotations about the same axis.

Two rotations of the same angle about the same signed Pauli axis merge into one
rotation of twice the angle. A rotation followed by its inverse (opposite-sign
axis) — or any pair of pi/2 rotations about the same axis, which compose to a
global phase — cancels entirely: the empty tuple `()` is returned so `traversal`
deletes both operations. Returns `nothing` when the pair cannot be merged.
"""
function merge_rotations(op1::CircuitOp.Type, op2::CircuitOp.Type)
    @match (op1,op2) begin
        (ExpEighPiPauli(),ExpEighPiPauli()) => _merge_same_axis(op1, op2, ExpQuatPiPauli)
        (ExpQuatPiPauli(),ExpQuatPiPauli()) => _merge_same_axis(op1, op2, ExpHalfPiPauli)
        (ExpHalfPiPauli(),ExpHalfPiPauli()) => begin
            (p1,p2) = complete_paulis(op1,op2)
            if p1.xz == p2.xz
                # Two pi/2 rotations about the same axis compose to +/-identity
                # (a global phase), regardless of the axis signs
                return ()
            else return nothing
            end
        end
        _ => nothing
    end
end
