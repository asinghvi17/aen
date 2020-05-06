function construct_system(systems::Vector{ODESystem}, graph::LightGraphs.AbstractGraph, connector::Function)

    couplings = Equation[]

    edgs = edges(graph)
    graph_weights = SparseArrays.SparseMatrixCSC(weights(graph))
    for vertex in vertices(graph)
        sys_to = systems[vertex]
        weights = []
        systems_from = ODESystem[]
        for neighbor in inneighbors(graph, vertex)
            push!(systems_from, systems[neighbor])
            push!(weights, graph_weights[neighbor, vertex])
        end
        append!(
            couplings,
            couple(
                systems[vertex],
                systems_from,
                weights
            )
        )
    end

    return ODESystem(
        couplings,
        t,
        [],
        [];
        systems = systems
    )
end

function construct_system(eqs::Vector{Equation}, graph::LightGraphs.AbstractGraph, connector::Function)
    return construct_system(
        [
            ODESystem(eqs)
            for _ in vertices(graph)
        ],
        graph,
        connector
    )
end
