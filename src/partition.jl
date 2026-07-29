"""
I need a function to expand a node and return expanded nodes(addition to frontier)
"""
function node_expand(node::CompilerState{TraversalRuntime})
    len=length(node.classical_register)
    if node.instruction_pointer>len
        return nothing
    else
        p1_state = copy(node)
        @reset p1_state.runtime = TraversalRuntime(1)
        m1_state = copy(node)
        @reset m1_state.runtime = TraversalRuntime(0)
        p1_state = do_quantum_step(p1_state)
        m1_state = do_quantum_step(m1_state)
        resolve_conditionals(p1_state)
        resolve_conditionals(m1_state)
        return (p1_state, m1_state)
    end
end

##
"""BFS function"""
function breath_first_search(state::CompilerState{TraversalRuntime})
    frontier = CompilerState[]
    len = 2^length(state.classical_register)
    push!(frontier,state)
    while length(frontier)<len
        node = popfirst!(frontier)
        expanded_nodes = node_expand(node)
        if isnothing(expanded_nodes)
            break
        else
            append!(frontier,node_expand(node))
        end
    end
    return frontier
end
