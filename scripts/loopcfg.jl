using DifferentialEquations, Plots

@inline cubic(v) = v - v^3/3
@inline pwl(v) = min(max(-2 - v, v), 2 - v)

# plot(-2:0.001:2, cubic)
# plot!(-2:0.001:2, pwl)

"""
    f(v)

The nonlinear voltage function, for FitzHugh-Nagumo
"""
@inline function f(v)
    return pwl(v)
end

"""
    CPG!(du, u, p, t)

In-place function to calculate the derivatives of the simplified CPG.
"""
function CPG!(du, u, p, t)

    # unpack the state array

    vxs = @view(u[1:4:end]) # 4-period of states.
    dxs = @view(u[2:4:end])
    vys = @view(u[3:4:end])
    dys = @view(u[4:4:end])

    e0, e1, g, b0, b1, Dhead, Drest, Dgap, J = p # unpack the parameter array

    # Head oscillator
    # neurons 0 (ventral and dorsal) coupled by 1-way diffusion,
    # negative diffusion constant `Dhead` simulated inhibitory synapse.
    du[1] = f(vxs[1]) - vys[1] + Dhead * dxs[1] + J
    du[2] = e0 * (vxs[1] - g * vys[1] + b0)
    du[3] = f(dxs[1]) - dys[1] + Dhead * vxs[1] + J
    du[4] = e0 * (dxs[1] - g * dys[1] + b0)


    # @infiltrate

    # body oscillators
    for i in 2:length(vxs)
        ind = (i-1)*4 + 1 # account for the periodicity of the array.
        du[ind]   = f(vxs[i]) - vys[i] + Drest * dxs[i] + Dgap * (vxs[i-1] - vxs[i]) + J
        du[ind+1] = e1 * (vxs[i] - g * vys[i] + b1)
        du[ind+2] = f(dxs[i]) - dys[i] + Drest * vxs[i] + Dgap * (dxs[i-1] - dxs[i]) + J
        du[ind+3] = e1 * (dxs[i] - g * dys[i] + b1)
    end

    return du


end


u0 = [
    1.0, -0.49, 1.0, -0.51,
    1.0, -0.5, 1.0, -0.5,
    1.0, -0.5, 1.0, -0.5,
    1.0, -0.49, 1.0, -0.51,
    1.0, -0.5, 1.0, -0.5,
    1.0, -0.5, 1.0, -0.5
]

## parameters
e0 = 0.08
e1 = 0.08
g = 0.8
b0 = 0.46
b1 = 0.47
Dhead = -0.2 # diffusion constant
Drest= -0.02
Dgap=0.03
J = 0.5

p = [e0, e1, g, b0, b1, Dhead, Drest, Dgap, J]


Juno.@enter CPG!(zeros(length(u0)), u0, p, 0.0)
## times
tspan = (0.0, 500.0)

prob = ODEProblem(CPG!, u0, tspan, p)

sol = solve(prob)

plot(sol; vars = 1:2:4, tspan = (0, 250), legend = :none)
savefig("lol.png")
