module PBCCompiler

# `run` and `to_result` are deliberately not exported: `run` would shadow
# `Base.run` for downstream users; call them qualified as `PBCCompiler.run`.
export Circuit, CircuitOp, MeasurementResult,
    AbstractRuntime, SimRuntime, StabilizerRuntime, DummyRuntime, TraversalRuntime,
    CompilerState, CompilationResult,
    preprocess_circuit, traversal, affectedqubits,
    parse_input, save_result, load_result,
    random_test_circuit,
    get_distribution, get_graph, get_hypergraph, weight_std_graph,
    circuitplot, circuitplot!, circuitplot_axis,
    circuitplot_quantikz, circuitstring_quantikz

include("type.jl")
include("io.jl")
include("traversal.jl")
include("affectedqubits.jl")
include("plotting.jl")
include("pair_transformation.jl")
include("preprocess.jl")
include("random_circuit.jl")
include("joint_measurement_check.jl")
include("logic.jl")
include("statistics.jl")

end # module PBCCompiler
