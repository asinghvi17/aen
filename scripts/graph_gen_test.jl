using LightGraphs, SimpleWeightedGraphs, ModelingToolkit, DifferentialEquations, Plots
using LightGraphs: src, dst, weights
using SimpleWeightedGraphs: SimpleWeightedDiGraph, SimpleWeightedEdge

# define the connectome graph

cpg = SimpleWeightedDiGraph(12)

Dinhibs = fill(0.1, 6)
Dgaps = fill(-0.01, 5)

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


function couple(sys_from, sys_to, weight)
    return [0 ~ sys_from.F - weight * sys_to.v]
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

couplings = Equation[]

for edge in edges(cpg)
    append!(
        couplings,
        couple(
            systems[src(edge)],
            systems[dst(edge)],
            weight(edge)
        )
    )
end

connected = ODESystem(
    couplings,
    t,
    [],
    [];
    systems = systems
)

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
