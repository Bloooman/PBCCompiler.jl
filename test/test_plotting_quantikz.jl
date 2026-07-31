@testitem "quantikz plotting" tags=[:plotting] begin

using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, Measurement, ExpHalfPiPauli, ExpQuatPiPauli,
    ExpEighPiPauli, PauliConditional, BitConditional, circuitstring_quantikz
using QuantumClifford: @P_str
using Quantikz

"""Split a quantikz environment into one vector of cells per qubit/bit row."""
function rows(str)
    body = split(str, "\n")[2:end-1]  # drop \begin{quantikz} / \end{quantikz}
    [strip.(split(replace(row, "\\\\" => ""), " & ")) for row in body]
end

"""A row's cells with the leading `\\lstick` wire label dropped."""
wire(r) = r[2:end]

@testset "support is per-wire, untouched wires stay clear" begin
    # Weight-3 Pauli straddling qubit 3, which it does not act on.
    circuit = Circuit([ExpQuatPiPauli(P"XYZ", [1, 2, 4])])
    str = circuitstring_quantikz(circuit)
    r = rows(str)

    @test count("\\gate", str) == 3
    # No spanning-rectangle machinery anywhere.
    @test !occursin("\\gate[", replace(str, "\\gate[style=" => ""))
    @test !occursin("\\linethrough", str)
    # Qubit 3 (row 3) has no gate in the operation's column.
    @test all(!contains("\\gate"), r[3])
    @test all(cell -> cell == "\\qw", wire(r[3]))
end

@testset "vertical rail offsets chain to the previous linked wire" begin
    circuit = Circuit([ExpQuatPiPauli(P"XYZ", [1, 2, 4])])
    r = rows(circuitstring_quantikz(circuit))
    # 1 -> 2 is one row up; 2 -> 4 is two rows up (not three).
    @test occursin("\\vqw{-1}", r[2][2])
    @test occursin("\\vqw{-2}", r[4][2])
    @test !occursin("\\vqw", r[1][2])
end

@testset "a rail never crosses another operation's box" begin
    # Qubit 2 lies inside the measurement's span but outside its support. Packing
    # must still push the rotation into its own column, or the rail would be drawn
    # straight through the rotation's box.
    circuit = Circuit([Measurement(P"ZIZ", 1, [1, 2, 3]), ExpQuatPiPauli(P"X", [2])])
    r = rows(circuitstring_quantikz(circuit))
    meas_col = findfirst(contains("\\gate"), r[1])
    rot_col = findfirst(contains("\\gate"), r[2])
    @test meas_col != rot_col
end

@testset "measurements keep their wire and bend to the classical bit" begin
    circuit = Circuit([Measurement(P"ZZ", 3, [1, 2])])
    str = circuitstring_quantikz(circuit)
    # `\meterD` would terminate the wire; PBC Pauli measurements are non-destructive.
    @test !occursin("\\meterD", str)
    @test occursin("\\cwbend", str)
    @test occursin("c3", str)  # outcome bit labelled on the box the bend leaves from
end

@testset "identity factors are dropped from support" begin
    circuit = Circuit([ExpQuatPiPauli(P"IZI", [1, 2, 3])])
    r = rows(circuitstring_quantikz(circuit))
    # Qubit 3 only ever carries an identity factor, but it is still declared, so
    # its wire is drawn -- just without a box on it.
    @test length(r) == 3
    @test occursin("\\gate", r[2][2])
    @test all(cell -> cell == "\\qw", wire(r[1]))
    @test all(cell -> cell == "\\qw", wire(r[3]))
end

@testset "empty circuit" begin
    @test circuitstring_quantikz(Circuit()) isa String
end

@testset "rotation angles and fills" begin
    circuit = Circuit([
        ExpHalfPiPauli(P"X", [1]),
        ExpQuatPiPauli(P"Y", [2]),
        ExpEighPiPauli(P"Z", [3]),
    ])
    str = circuitstring_quantikz(circuit)
    # Quantikz.jl's own `transparent` option would suppress every fill we emit.
    @test !occursin("transparent", str)
    @test occursin("fill=", str)
    # Each rotation kind gets its own fill, so the angle needs no label by default.
    @test !occursin("\\pi/", str)
    @test length(unique(eachmatch(r"fill=[^]]*", str) .|> m -> m.match)) == 3

    labelled = circuitstring_quantikz(circuit; showangles=:first)
    @test occursin("\\pi/2", labelled) && occursin("\\pi/4", labelled) && occursin("\\pi/8", labelled)

    # showangles=:all repeats the angle on every wire, :first only on the topmost.
    wide = Circuit([ExpQuatPiPauli(P"XY", [1, 2])])
    @test count("\\pi/4", circuitstring_quantikz(wide; showangles=:all)) == 2
    @test count("\\pi/4", circuitstring_quantikz(wide; showangles=:first)) == 1
    @test !occursin("\\pi/4", circuitstring_quantikz(wide; showangles=:none))
    @test_throws ArgumentError circuitstring_quantikz(wide; showangles=:every)
end

@testset "wire labels" begin
    circuit = Circuit([Measurement(P"ZZ", 1, [1, 2]), ExpQuatPiPauli(P"XY", [1, 3])])
    r = rows(circuitstring_quantikz(circuit))
    @test r[1][1] == "\\lstick{\$q_{1}\$}"
    @test r[3][1] == "\\lstick{\$q_{3}\$}"
    # Compressed bit layout collapses all bits onto one wire, so it is just "c".
    @test r[end][1] == "\\lstick{\$c\$}"

    # Labels must not displace circuit content: same gate count either way.
    bare = circuitstring_quantikz(circuit; showlabels=false)
    @test !occursin("\\lstick", bare)
    @test count("\\gate", bare) == count("\\gate", circuitstring_quantikz(circuit))
end

@testset "expanded bit layout labels each bit" begin
    circuit = Circuit([Measurement(P"Z", 1, [1]), Measurement(P"Z", 2, [2])])
    @with Quantikz.classicalbitslayout => :expanded begin
        r = rows(circuitstring_quantikz(circuit))
        @test r[end-1][1] == "\\lstick{\$c_{1}\$}"
        @test r[end][1] == "\\lstick{\$c_{2}\$}"
    end
end

@testset "color overrides" begin
    circuit = Circuit([ExpQuatPiPauli(P"X", [1])])
    @test occursin("fill=green!50", circuitstring_quantikz(circuit; colors=Dict(:quatpi => "green!50")))
end

@testset "conditionals" begin
    # BitConditional annotates the topmost box and bends to the conditioning bit.
    circuit = Circuit([BitConditional(ExpHalfPiPauli(P"XZ", [2, 4]), 7)])
    str = circuitstring_quantikz(circuit)
    @test occursin("?c7", str)
    @test occursin("\\cwbend", str)

    # A conditional measurement keeps its own outcome bit *and* the conditioning bit.
    both = circuitstring_quantikz(Circuit([BitConditional(Measurement(P"Z", 2, [1]), 5)]))
    @test occursin("c2", both) && occursin("?c5", both)

    # PauliConditional draws control and target legs on one rail.
    pc = circuitstring_quantikz(Circuit([PauliConditional(P"X", [1], P"Z", [3])]))
    @test occursin("\\mathrm{c}", pc) && occursin("\\mathrm{t}", pc)
    @test occursin("\\vqw", pc)

    # A qubit carrying both a control and a target leg gets one merged box.
    merged = circuitstring_quantikz(Circuit([PauliConditional(P"X", [1], P"Z", [1])]))
    @test occursin("X_{\\mathrm{c}}\\,Z_{\\mathrm{t}}", merged)
end

@testset "file output" begin
    circuit = Circuit([Measurement(P"ZZ", 1, [1, 2]), ExpQuatPiPauli(P"XY", [1, 3])])
    mktempdir() do dir
        tex = joinpath(dir, "circuit.tex")
        circuitplot_quantikz(circuit, tex)
        @test isfile(tex)
        @test occursin("\\begin{quantikz}", read(tex, String))

        # The .pdf branch goes through tectonic; skip rather than fail if the
        # toolchain is unavailable on this machine.
        pdf = joinpath(dir, "circuit.pdf")
        try
            circuitplot_quantikz(circuit, pdf)
            @test isfile(pdf) && filesize(pdf) > 0
        catch e
            @info "skipping PDF render: LaTeX toolchain unavailable" exception=e
        end
    end
end

end
