
# ## Parameters and initial conditions

# The parameter structure is the following:
# - `e`, `b`, `g` and `J` specified per neuron
# - `Dgap` and `Dinhib` specified per neuron in the descending pathway - they control
#   the connections coming in to that neuron
# - `Kgap` and `Kinhib` constants which control asymmetry, specified per neuron, multiply with previous
# - The head oscillator does not use `Dgap` or `Kgap`; other than that, it's a normal neuron.

NUM_NEURONS = 24

function construct_params(
        dinh_ind,
        dinh_val,
        dgap_ind,
        dgap_val;
        e = 0.04,
        b = 0.47,
        g = 0.8,
        J = 0.0,
        Dgap = 0.05,
        Dinhib = -0.02,
        Dhead = -0.1
    )
    es      = fill(e,   NUM_NEURONS)
    bs      = fill(b,   NUM_NEURONS)
    gs      = fill(g,    NUM_NEURONS)
    Js      = fill(J,    NUM_NEURONS)
    Dgaps   = fill(Dgap,   NUM_NEURONS)
    Dinhibs = fill(Dinhib,  NUM_NEURONS)
    Kgaps   = fill(1.0,    NUM_NEURONS)
    Kinhibs = fill(1.0,    NUM_NEURONS)

    # These have to be special-cased for the head oscillator.
    Dinhibs[1:2] .= Dhead
    # Js[3:end] .= 0.1

    # Here, we apply the overrides in the function description.
    Kinhibs[dinh_ind] = dinh_val
    Kgaps[dgap_ind] = dgap_val
    # This puts all the parameters together into a Vector of Vectors, which is useful to unpack.
    ps = [es, bs, gs, Js, Dgaps, Dinhibs, Kgaps, Kinhibs]
    return ps
end

# We define the initial conditions

# These are for the head oscillator,
u0 = [
    1.0,    # ventral V₀
    -0.667, # ventral W₀
    -1.0,   #  dorsal V₁
    -0.7    #  dorsal W₁
]

# and these define the descending pathway.

for neuron in 3:NUM_NEURONS
    append!(
        u0,
        [1.0, -0.5]
    )
end


# ## Core function

f(v) = #=v - v^3 / 3=# min(max(-2 - v, v), 2 - v)

function CPG!(du, u, p, t)

    num_neurons = length(u) ÷ 2
    num_pairs   = num_neurons ÷ 2

    # `@inbounds` disables bounds checking, making
    # array accesses much faster.
    @inbounds begin

    # First, we unpack the parameters.
    es, bs, gs, Js, Dgaps, Dinhibs, Kgaps, Kinhibs = p

    # Now, we get to the calculations.

    # First, we handle the head oscillator.
    # 1 is the ventral neuron,

    du[1] = f(u[1]) - u[2] + Kinhibs[1] * Dinhibs[1] * u[3] + Js[1]

    du[2] = es[1] * (u[1] - gs[1] * u[2] + bs[1])

    # and 2 is the dorsal neuron.

    du[3] = f(u[3]) - u[4] + Kinhibs[2] * Dinhibs[2] * u[1] + Js[2]

    du[4] = es[2] * (u[3] - gs[2] * u[4] + bs[2])

    # Now, we handle the descending pathway.
    # The symmetry of the C. elegans CPG means that
    # the loop is surprisingly simple, and we can
    # enable some performance optimizations.
    for i in 2:(num_neurons - 1)
        nn = i * 2 + 1
        v = u[nn]         # voltage channel of the neuron
        w = u[nn + 1]     # membrane potential of the neuron
        pre_v = u[nn - 4] # voltage of the neuron to the left (which drives it through a gap jcn)
        opp_v = u[nn - 2] # voltage of the neuron opposite (mutual inhibitory coupling)

        du[nn] = f(v) - w + Kinhibs[i] * Dinhibs[i] * opp_v + Kgaps[i] * Dgaps[i] * max(pre_v - v, 0) + Js[i]

        du[nn + 1] = es[i] * (v - gs[i] * w + bs[i])

    end
    end
    return du
end

# ## Initial solution

tspan = (0.0, 3000.0)

prob = ODEProblem(CPG!, u0, tspan, construct_params(1, 1.0, 1, 1.0; J = .238))
sol = solve(prob, Tsit5(); reltol = 1e-5, abstol = 1e-5)

Plots.plot(sol; vars = collect(3:4:NUM_NEURONS))
Plots.plot(sol; vars = collect(1:4:NUM_NEURONS))
