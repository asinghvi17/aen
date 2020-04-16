
using DifferentialEquations,
      ParameterizedFunctions

using Makie, Observables

# create ode
# The parameter structure is the following:
# - `e`, `b`, `g` and `J` specified per neuron
# - `Dgap` and `Dinhib` specified per neuron in the descending pathway - they control
#   the connections coming in to that neuron
# - `Kgap` and `Kinhib` constants which control asymmetry, specified per neuron, multiply with previous
# - The head oscillator does not use `Dgap` or `Kgap`; other than that, it's a normal neuron.

es      = fill(0.04,   12)
bs      = fill(0.0,    12)
gs      = fill(0.8,    12)
Js      = fill(0.5,    12)
Dgaps   = fill(0.03,   12)
Dinhibs = fill(-0.02, 12)
Kgaps   = fill(1.0,    12)
Kinhibs = fill(1.0,    12)

# These have to be special-cased for the head oscillator.
Dinhibs[1:2] .= -0.047
bs[1:2] .= 0
ps = [es, bs, gs, Js, Dgaps, Dinhibs, Kgaps, Kinhibs]

# We define the initial conditions;

# These are for the head oscillator,
vent_v0 = -2.0
vent_w0 = -0.667
dors_v0 = -2.0
dors_w0 = -0.7

# and these define the descending pathway.
vent_v1 = -2.0
vent_w1 = -0.5
dors_v1 = -2.0
dors_w1 = -0.5
vent_v2 = -2.0
vent_w2 = -0.5
dors_v2 = -2.0
dors_w2 = -0.5
vent_v3 = -2.0
vent_w3 = -0.49
dors_v3 = -2.0
dors_w3 = -0.51
vent_v4 = -2.0
vent_w4 = -0.5
dors_v4 = -2.0
dors_w4 = -0.5
vent_v5 = -2.0
vent_w5 = -0.5
dors_v5 = -2.0
dors_w5 = -0.5

u0 = [vent_v0,vent_w0,dors_v0,dors_w0,vent_v1,vent_w1,dors_v1,dors_w1,vent_v2, vent_w2, dors_v2, dors_w2, vent_v3, vent_w3, dors_v3, dors_w3, vent_v4, vent_w4, dors_v4, dors_w4, vent_v5, vent_w5, dors_v5, dors_w5 ]


f(v) = min(max(-2 - v, v), 2 - v)

function CPG!(du, u, p, t)

    # du = zeros(length(u))

    # `@inbounds` disables bounds checking, making
    # array accesses much faster.

    # First, we unpack the parameters and the
    # states.  The first chunk handles parameters,
    # and the second and third chunks handle
    # derivatives and states.
    # The `@view` macro makes a "view" of an array
    # without copying it, but it can have nonlinear
    # indexing; this is very convenient.
    @inbounds begin
    es, bs, gs, Js, Dgaps, Dinhibs, Kgaps, Kinhibs = p

    # Now, we get to the calculations.

    # First, we handle the head oscillator.
    # first, the ventral side

    du[1] = f(u[1]) - u[2] + Dinhibs[1] * max(u[3] - u[1], 0)

    du[2] = es[1] * (u[1] - gs[1] * u[2] + bs[1])

    du[3] = f(u[3]) - u[4] + Dinhibs[2] * max(u[1] - u[2], 0)

    du[4] = es[2] * (u[3] - gs[2] * u[4] + bs[2])
    # du[1] = f(vs[1])
    #     - ws[1]
    #     + Kinhibs[1] * Dinhibs[1] * max(vs[2] - vs[1], 0)
    #     + Js[1]
    #
    # du[2] = es[1] * (vs[1] - gs[1] * ws[1] + bs[1])
    #
    # du[3] = f(vs[2])
    #     - ws[2]
    #     + Kinhibs[2] * Dinhibs[2] * max(vs[1] - vs[2], 0)
    #     + Js[2]
    #
    # du[4] = es[2] * (vs[2] - gs[2] * ws[2] + bs[2])

    # Now, we handle the descending pathway.
    # The symmetry of the C. elegans CPG means that
    # the loop is surprisingly simple, and we can
    # enable several performance optimizations.
    for i in 3:12
        nn = i * 2
        v = u[nn]
        w = u[nn + 1]
        pre_v = u[nn - 4]
        opp_v = u[nn - 2]
        du[nn] = f(v)
            - w
            + Kinhibs[i] * Dinhibs[i] * max(opp_v - v, 0)
            + Kgaps[i] * Dgaps[i] * max(pre_v - v, 0)
            + Js[i]

        du[nn + 1] = es[i] * (v - gs[i] * w + bs[i])

    end
    end
    return du
end

# function FHND_(u, p, t)
#     es, bs, gs, Js, Dgaps, Dinhibs, Kgaps, Kinhibs = p
#
#     return [
#         f(u[1]) - u[2] + Dinhibs[1] * max(u[3] - u[1], 0),
#
#         es[1] * (u[1] - gs[1] * u[2] + bs[1]),
#
#         f(u[3]) - u[4] + Dinhibs[2] * max(u[1] - u[2], 0),
#
#         es[2] * (u[3] - gs[2] * u[4] + bs[2]),
#     ]
# end

prob = ODEProblem(CPG!, u0, (0.0, 2000.0), ps)

sol = solve(prob)

using MakieLayout

oscene, layout = layoutscene()

ax = layout[1, 1] = LAxis(oscene)

Dkr = LinRange(0, 1, 100)

# sl = layout[2, 1] = LSlider(oscene; range = Dkr, startvalue = .5)
Dk = Observable(1.0)

oscene



x0 = -2
y0 =-0.6667
x1 = -2.0
y1 = -0.7

u0 = [x0, y0, x1, y1]

D  = -0.047
g  = 0.8
e0 = 0.04
b0 = 0
e1 = 0.04
b1 = 0.2 # was 0.25

txt = layout[0, 1] = LText(oscene, @lift("Dₙₑₜ = " * string(round(D * $Dk; sigdigits=3))); width = Auto(false), textsize = 30)

p = @lift [D, $Dk, g, e0, b0, e1, b1]

ts = LinRange(1000, 1750, 10000)

prob = ODEProblem(FHND, u0, (000.0, 2000.0), to_value(p), AutoTsit5(Vern9()))

sol = @lift solve(remake(prob; p = $p))(ts)

plot_v0 = @lift $sol[1, :]
plot_v1 = @lift $sol[3, :]

lines!(ax, ts, plot_v0; color=:red)
lines!(ax, ts, plot_v1; color=:blue)

record(oscene, "coupling.mp4", LinRange(0, 1, 240)) do i
    Dk[] = i
end
