
# We define the initial conditions

# These are for the head oscillator,
u0 = [
    1.0,    # ventral V₀
    -0.667, # ventral W₀
    1.0,    #  dorsal V₁
    -0.7    #  dorsal W₁
]

# and these define the descending pathway.

for neuron in 3:NUM_NEURONS
    append!(
        u0,
        [1.0, -0.5]
    )
end

# ## Initial solution

tspan = (0.0, 3000.0)

p0 = construct_params(
    1, 1.0, 1, 1.0;
    # Dgap = 0.1,
    # J = 0.01
)

p0[4][1:2:end] .= -.2

prob = ODEProblem(CPG!, u0, tspan, p0)
sol = solve(prob, Tsit5(); reltol = 1e-5, abstol = 1e-5)

Plots.plot(
    Plots.plot(sol; vars = 1:4:NUM_NEURONS, title = "Ventral neurons"),
    Plots.plot(sol; vars = 3:4:NUM_NEURONS, title = "Dorsal neurons");
    layout = (2,1)
)
