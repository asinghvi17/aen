using LightGraphs, MetaGraphs

using DifferentialEquations, ModelingToolkit

t = Variable(:t; known = true)()
σ = Variable(:σ; known = true)()
ρ = Variable(:ρ; known = true)()
β = Variable(:β; known = true)()
x = Variable(:x)(t)
y = Variable(:y)(t)
z = Variable(:z)(t)
D = Differential(t)

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

de = ODESystem(eqs)

f = ODEFunction(de, [x, y, z], [σ,ρ,β]) # pass vars and then params....

prob = ODEProblem(f, [1.0,0.0,0.0], (0.0, 100.0), [10.0,28.0,8/3])

sol = solve(prob)

using Plots

plot(sol, vars = (1, 2, 3))

struct ODESpec

    lhs

    rhs

    coupling

end

struct ODENode

    lhs

    rhs

    params

    coupling

end

ODENode(ODESpec, params) = ODENode(ODESpec.lhs, ODESpec.rhs, params, ODESpec.coupling)

function couple!(o1::ODESys, o2::ODESys, D::Number)
    return
end
