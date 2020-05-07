using LightGraphs, SimpleWeightedGraphs, SparseArrays, ModelingToolkit, DifferentialEquations, Plots
using LightGraphs: src, dst, weights
using SimpleWeightedGraphs: SimpleWeightedDiGraph, SimpleWeightedEdge

# define the connectome graph

cpg = SimpleWeightedDiGraph(12)

Dinhibs = fill(0.05, 6)
Dinhibs[1] = 0.1 # head oscillator
Dgaps = fill(-0.05, 5)

# Add mutual inhibition
for (ventral, dorsal, Dinhib) in zip(1:6, 7:12, Dinhibs)
    e1 = SimpleWeightedEdge(ventral, dorsal, Dinhib)
    e2 = SimpleWeightedEdge(dorsal, ventral, Dinhib)
    LightGraphs.add_edge!(cpg, e1)
    LightGraphs.add_edge!(cpg, e2)
end

# We handle gaps separately, since
# the head oscillator does not have
# a forcing neuron pair upstream
for (ventral, dorsal, Dgap) in zip(2:6, 8:12, Dgaps)
    e1 = SimpleWeightedEdge(ventral - 1, ventral, Dgap)
    e2 = SimpleWeightedEdge(dorsal - 1, dorsal, Dgap)
    LightGraphs.add_edge!(cpg, e1)
    LightGraphs.add_edge!(cpg, e2)
end

# define the unit equation

@parameters t g e b

@variables v(t) w(t) F(t)

@derivatives D'~t

f(v) = v - v^3/3

fhn = [
    D(v) ~ f(v) - w + F,
    D(w) ~ e * (v - g * w + b)
]


function couple(system_to, systems_from, weights)
    return [0 ~ system_to.F +  sum(weights .* getproperty.(systems_from, :v))]
end

# generate a single equation spec from the graph

# first, create all relevant systems
systems = [
    ODESystem(
        fhn, t, [v,w,F], [g,e,b];
        name = Symbol(
                "n" * string(i) # ModelingToolkit.map_subscripts(string(i))
        )
    )
    for i in vertices(cpg)
]

# populate couplings from edge

# couplings = Equation[]
#
# edgs = edges(cpg)
# graph_weights = SparseArrays.SparseMatrixCSC(weights(cpg))
# for vertex in vertices(cpg)
#     sys_to = systems[vertex]
#     weights = []
#     systems_from = ODESystem[]
#     for neighbor in inneighbors(cpg, vertex)
#         push!(systems_from, systems[neighbor])
#         push!(weights, graph_weights[neighbor, vertex])
#     end
#     append!(
#         couplings,
#         couple(
#             systems[vertex],
#             systems_from,
#             weights
#         )
#     )
# end
#
#
# connected = ODESystem(
#     couplings,
#     t,
#     [],
#     [];
#     systems = systems
# )

include("graph_system_constructor.jl")
connected = construct_system(systems, cpg, couple)

# Construct the initial condition, special-casing the first
# 2 neurons which comprise the head oscillator

u0 = [
    systems[1].v => -1.0,
    systems[1].w => -0.51,
    systems[1].F => -0.01,
    systems[2].v =>  1.0,
    systems[2].w => -0.49,
    systems[2].F => 0.01,
]


for sys in systems[3:end]
    append!(u0,
        [
        sys.v =>  1.0,
        sys.w => -0.49,
        sys.F => 0.01,
        ]
    )
end

# Construct the parameters, which are the same across
# all systems in this case

p0 = Pair{Operation, Float64}[]

for sys in systems
    append!(
        p0,
        [
            sys.g => 0.8,
            sys.b => 0.46,
            sys.e => 0.04,
        ]
    )
end

# Construct the ODEProblem

prob = ODEProblem(connected, u0, (0.0, 2000.0), p0)
sol = solve(prob, Rodas5())

p = palette(:diverging_linear_bjr_30_55_c53_n256, 14)[[(1:6)..., (9:14)...]]


Plots.plot(sol; vars = collect(1:3:36), tspan = (1600, 2000), palette = p, legend = :outerright)

##

# graphplot(cpg; names = 1:12)
