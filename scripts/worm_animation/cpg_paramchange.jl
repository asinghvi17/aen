tspan = (0.0, 3000.0)


function condition(u, t, integrator)
    return (t % 500) - 250
end

function affect!(integrator)
    integrator.p = construct_params(
        1, 1.0, 1, 1.0;
        e = integrator.p[1][1] + .01
    )
end

cb = ContinuousCallback(condition, affect!)

prob = ODEProblem(CPG!, u0, tspan, construct_params(1, 1.0, 1, 1.0))

sol = solve(prob, Tsit5(); reltol = 1e-5, abstol = 1e-5, callback = cb)

p1 = Plots.plot(sol; vars = 3:4:NUM_NEURONS, title = "Original solution")
p2 = Plots.plot(sol_pc; vars = 3:4:NUM_NEURONS, title = "Solution with parameter change")
Plots.plot(p1, p2; layout = (2, 1), size = (1200, 1000))
