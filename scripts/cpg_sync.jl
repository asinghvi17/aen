# # An interactive CPG dashboard

using DifferentialEquations, ImageFiltering, Interpolations

using Makie, Observables

# ## Parameters and initial conditions

# The parameter structure is the following:
# - `e`, `b`, `g` and `J` specified per neuron
# - `Dgap` and `Dinhib` specified per neuron in the descending pathway - they control
#   the connections coming in to that neuron
# - `Kgap` and `Kinhib` constants which control asymmetry, specified per neuron, multiply with previous
# - The head oscillator does not use `Dgap` or `Kgap`; other than that, it's a normal neuron.

function construct_params(dinh_ind, dinh_val, dgap_ind, dgap_val)
    es      = fill(0.04,   12)
    bs      = fill(0.47,   12)
    gs      = fill(0.8,    12)
    Js      = fill(0.0,    12)
    Dgaps   = fill(0.05,   12)
    Dinhibs = fill(-0.02,  12)
    Kgaps   = fill(1.0,    12)
    Kinhibs = fill(1.0,    12)

    # These have to be special-cased for the head oscillator.
    Dinhibs[1:2] .= -0.1
    # Js[3:end] .= 0.1

    # Here, we apply the overrides in the function description.
    Kinhibs[dinh_ind] = dinh_val
    Kgaps[dgap_ind] = dgap_val
    # This puts all the parameters together into a Vector of Vectors, which is useful to unpack.
    ps = [es, bs, gs, Js, Dgaps, Dinhibs, Kgaps, Kinhibs]
    return ps
end

ps = construct_params(1, 1.0, 1, 1.0)
# We define the initial conditions

# These are for the head oscillator,
vent_v0 = -1.0
vent_w0 = -0.667
dors_v0 = -1.0
dors_w0 = -0.7

# and these define the descending pathway.
vent_v1 = -1.0
vent_w1 = -0.5
dors_v1 = -1.0
dors_w1 = -0.5
vent_v2 = -1.0
vent_w2 = -0.5
dors_v2 = -1.0
dors_w2 = -0.5
vent_v3 = -1.0
vent_w3 = -0.49
dors_v3 = -1.0
dors_w3 = -0.51
vent_v4 = -1.0
vent_w4 = -0.5
dors_v4 = -1.0
dors_w4 = -0.5
vent_v5 = -1.0
vent_w5 = -0.5
dors_v5 = -1.0
dors_w5 = -0.5

u0 = [vent_v0,vent_w0,dors_v0,dors_w0,vent_v1,vent_w1,dors_v1,dors_w1,vent_v2, vent_w2, dors_v2, dors_w2, vent_v3, vent_w3, dors_v3, dors_w3, vent_v4, vent_w4, dors_v4, dors_w4, vent_v5, vent_w5, dors_v5, dors_w5 ]

tspan = (0.0, 3000.0)
# ## Core function

f(v) = v - v^3 / 3 # min(max(-2 - v, v), 2 - v)

function CPG!(du, u, p, t)

    # du = zeros(length(u))

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
    for i in 2:11
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

prob = ODEProblem(CPG!, u0, tspan, ps)
sol = solve(prob; reltol = 1e-5, abstol = 1e-5)
# We've solved the equation; now we need plots.
import Plots
# This plots all the voltages:
Plots.plot(sol; vars = 1:2:24, tspan = (1000, 1600), palette = p)
Plots.plot(sol; vars = 1:2:4)
# this plots only the ventral voltages
Plots.plot(sol; vars = 1:4:24, tspan = 1000:2000)
# and this plots only the dorsal ones.
Plots.plot(sol; vars = 3:4:24, tspan = 1000:2000)

neuron_range = 1:12

function eff_signal(
        sol;
        σ = 40,
        kernel = reflect(ImageFiltering.Kernel.gaussian((σ,), (σ * 4 * 2 + 1,))),
        trange = (0, 2000),
        n_samples = 20_000
    )

    num_neurons = size(sol, 1) ÷ 2 # 2 eqs/neuron
    num_pairs   = num_neurons  ÷ 2 # 2 neurons per pair

    ts = sol(LinRange(trange..., n_samples))

    effective_signal = zeros(n_samples, num_pairs)

    for nn in 1:num_pairs
        neuron_number = nn - 1 # offset for indexing
        effective_signal[:, nn]  = (ts[4*neuron_number + 1, :]  .+ abs.(ts[4*neuron_number + 1, :])) ./ 2
        effective_signal[:, nn] -=(ts[4*neuron_number+2 + 1, :] .+ abs.(ts[4*neuron_number+2 + 1, :])) ./ 2
        effective_signal[:, nn]  = imfilter(effective_signal[:, nn], kernel)
    end

    return effective_signal
end

signal = eff_signal(sol; σ = 40)

function rotmat(θ)
    return AbstractPlotting.Mat2(cos(θ), sin(θ), -sin(θ), cos(θ))
end

function signal_to_worm(
        Δθ::AbstractVector;
        tightness_of_bend = 0.5,
    )

    xpos = 0
    ypos = 0

    xcoords = Float64[xpos]
    ycoords = Float64[ypos]

    θ = 0

    number_of_joints = length(Δθ)

    for j in 1:number_of_joints
        δθ = Δθ[j]
        θ += tightness_of_bend * δθ

        xpos += cos(θ)
        ypos += sin(θ)
        append!(xcoords, xpos)
        append!(ycoords, ypos)
    end
        # print(i,deltatheta, xpos, ypos)
        # @show θ

    # spline fit
    nodes=range(0, length = number_of_joints+1)

    csx = CubicSplineInterpolation(nodes,xcoords)
    csy = CubicSplineInterpolation(nodes,ycoords)

    worm_coords = range(0, stop = number_of_joints, step = 0.01)
    # I am incrementing the position in steps of
    # 0.01*(distance between nodes)

    smoothwormx=csx.(worm_coords)
    smoothwormy=csy.(worm_coords)

    points = Point2f0.(smoothwormx, smoothwormy)


    rotation = rotmat(-tightness_of_bend * Δθ[1])

    return Ref(rotation) .* points
end

using Makie
framerate = 1/40

linobs = Node(Point2f0[(0, 0)])
sc = lines(linobs; limits = Rect(-8, -8, 16, 16), resolution = (1000, 1000), linewidth = 10)
for i in 8000:8800
    linobs[] = signal_to_worm(signal[i, :])
    sleep(framerate)
end

# ## Dashboard
using Makie: lines, scatter
using MakieLayout

# Construct Scene and grid layout
scene, layout = layoutscene()

# Construct axes to plot into, and set them up
axs = layout[1:2, 1:3] = [LAxis(scene) for i in 1:6]
hidexdecorations!.(axs)
hideydecorations!.(axs)
linkxaxes!(axs...)
linkyaxes!(axs...)

# Construct sliders and text indicators
dinh_num = layout[4, 1:3] = LSlider(scene; range = 1:12, startval = 2)
dinh_val = layout[5, 1:3] = LSlider(scene; range = LinRange(0, 1, 100), startval = 1)
layout[3, 1:3] = LText(scene, lift((ind, val) -> "Dinhib = $val at $ind", dinh_num.value, dinh_val.value); textsize = 30, width = Auto(false))
dgap_num = layout[7, 1:3] = LSlider(scene; range = 2:11, startval = 2)
dgap_val = layout[8, 1:3] = LSlider(scene; range = LinRange(0, 1, 100), startval = 1)
layout[6, 1:3] = LText(scene, lift((ind, val) -> "Dgap = $val at $ind", dgap_num.value, dgap_val.value); textsize = 30, width = Auto(false))

# Construct time slider
t_sl = layout[1:2, 4] = LSlider(scene; range = LinRange(tspan..., 100), startval = 1000, horizontal = false,
    width = Auto(true), height = nothing)


# Set up infrastructure for diffeq solve
params = lift(construct_params, dinh_num.value, dinh_val.value, dgap_num.value, dgap_val.value)

tspan = (0.0, 2000.0)
prob = ODEProblem(CPG!, u0, tspan, ps)
oprob = @lift(remake(prob; p = $params))

tsv = t_sl.value
t = @lift LinRange(clamp($tsv - 350, tspan...), clamp($tsv + 350, tspan...), 1000)

sol = @lift(solve($oprob)($t))

vs = [lift(s -> getindex(s, i, :), sol) for i in 1:2:24]

# Panel 1 is the head oscillator
lines!(axs[1], vs[1]; color = :blue)
lines!(axs[1], vs[2]; color = :red)

lines!(axs[2], vs[3]; color = :blue)
lines!(axs[2], vs[4]; color = :red)

lines!(axs[3], vs[5]; color = :blue)
lines!(axs[3], vs[6]; color = :red)

lines!(axs[4], vs[7]; color = :blue)
lines!(axs[4], vs[8]; color = :red)

lines!(axs[5], vs[9]; color = :blue)
lines!(axs[5], vs[10]; color = :red)

lines!(axs[6], vs[11]; color = :blue)
lines!(axs[6], vs[12]; color = :red)
# if i is odd, then it is a single;
# else, it is a double.
# however, we need to be careful of which lines we color.



# cpg = layout[1, 2] = LAxis(scene)
# hidexdecorations!(cpg)
# hideydecorations!(cpg)
# cpg.xgridvisible[] = false
# cpg.ygridvisible[] = false
#
# nodes = Point2f0[
#     Point2f0(x, y) for y in 1:2, x in 1:6
# ]
#
#
# doubles = [nodes[1, i] => nodes[2, i] for i in 1:6]
#
# up = map(x -> first(x) + Point2f0(-0.1, 0) => last(x) + Point2f0(-0.1, 0), doubles)
# down = map(x -> first(x) + Point2f0(0.1, 0) => last(x) + Point2f0(0.1, 0), doubles)
#
# singles = vec([nodes[j, i-1] => nodes[j, i] for j in 1:2, i in 2:6])
#
# upsegs = linesegments!(cpg, up; color = fill(:black, length(up)))
# dnsegs = linesegments!(cpg, down; color = fill(:black, length(down)))
#
# c = fill(:black, length(singles))
# c[4] = :red
# sgsegs = linesegments(singles; color = c)
#
#
# lift(dgap_num.value) do n
#     c = fill(:black, length(singles))
#     c[n] = :red
#     sgsegs.color[] = c
# end
#
# lift(dinh_num.value) do n
#     seg, newn = if isodd(n) # n is going to an upper neuron
#         (upsegs, (n-1)÷2)
#     else
#         (dnsegs, n÷2)
#     end
#     c = fill(:black, length(seg))
#     c[newn] = :red
#     seg.color[] = c
# end
scene

linobs = Node(Point2f0[(0, 0)])

sc = lines(linobs; limits = Rect(-4, -6, 10, 12), resolution = (750, 900), linewidth = 10, scale_plot = false)

record(sc, "worm.mp4", 8_000:12_000; framerate = 60, sleep = false) do i
    linobs[] = signal_to_worm(signal[i, :])
end
