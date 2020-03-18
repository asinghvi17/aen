using DifferentialEquations
using RecurrenceAnalysis, ChaosTools, DynamicalSystems
using ImageFiltering, DataInterpolations
using DSP, DataFrames, Statistics
using PyCall, Conda
import Plots, PlotUtils
using Makie

# Pycall setup:
# Conda.add("scipy")

py"""
import scipy.ndimage.filters as filt
def pyfilter(x):
    return filt.gaussian_filter1d(x,40)
"""

pyfilter(x) = py"pyfilter"(x)

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

function CPG(u, p, t)

    vx0, vy0, dx0, dy0, vx1, vy1, dx1, dy1, vx2, vy2, dx2, dy2, vx3, vy3, dx3, dy3, vx4, vy4, dx4, dy4, vx5, vy5, dx5, dy5 = u

    e0, e1, g, b0, b1, Dhead, Drest, Dgap, J = p # unpack the parameter array

    dvdt=SVector(
          f(vx0)- vy0+ Dhead*dx0                  + J, e0*(vx0-g*vy0+b0),   # neurons 0 (ventral and dorsal)
          f(dx0)- dy0+ Dhead*vx0                  + J, e0*(dx0-g*dy0+b0),   # coupled by 1-way diffusion
                                                                                    # negative diffusion constant Dhead
                                                                                    # simulated inhibitory synapse

          f(vx1)- vy1+ Drest*dx1 + Dgap*(vx0-vx1) + J, e1*(vx1-g*vy1+b1),   # neurons 1 (ventral and dorsal)
          f(dx1)- dy1+ Drest*vx1 + Dgap*(dx0-dx1) + J, e1*(dx1-g*dy1+b1),   # driven by neurons 0 through one-way
                                                                                    # diffusion (selective gap junction)
                                                                                    # Dgap
                                                                                    # Drest simulates inhibitary synapse
                                                                                    # coupling ventral and dorsam neurons 1

          f(vx2)- vy2+ Drest*dx2 + Dgap*(vx1-vx2) + J, e1*(vx2-g*vy2+b1),   # neurons 2
          f(dx2)- dy2+ Drest*vx2 + Dgap*(dx1-dx2) + J, e1*(dx2-g*dy2+b1),

          f(vx3)- vy3+ Drest*dx3 + Dgap*(vx2-vx3) + J, e1*(vx3-g*vy3+b1),   # neurons 3
          f(dx3)- dy3+ Drest*vx3 + Dgap*(dx2-dx3) + J, e1*(dx3-g*dy3+b1),

          f(vx4)- vy4+ Drest*dx4 + Dgap*(vx3-vx4) + J, e1*(vx4-g*vy4+b1),   # neurons 4
          f(dx4)- dy4+ Drest*vx4 + Dgap*(dx3-dx4) + J, e1*(dx4-g*dy4+b1),

          f(vx5)- vy5+ Drest*dx5 + Dgap*(vx4-vx5) + J, e1*(vx5-g*vy5+b1),   # neurons 5
          f(dx5)- dy5+ Drest*vx5 + Dgap*(dx4-dx5) + J, e1*(dx5-g*dy5+b1))


    # diffusion is represented by the term D*(x1-⁠x0)
    return dvdt
end
# constants

## initial conditions
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

## times
tspan = (0.0, 2000.0)

# Define the problem

prob = ODEProblem(CPG, u0, tspan, p)

sol = solve(prob)

Plots.plot(sol; vars = 1:2:8, tspan = (500, 750), legend = :none)

phase_lines = sol[1:2, :]

phasepoints = Point2f0.(phase_lines[1, :], phase_lines[2, :])

tmin = Node(1)
tmax = Node(length(phasepoints))
trange = lift(((x, y) -> x:y), tmin, tmax)

scat = scatter(trange, lift(tr -> phasepoints[tr], trange))
# Ensemble simulation

trs = 1000

function mutate_prob!(prob, i, repeat)
    prob.p[end] = i/trs
    return prob
end

ens = EnsembleProblem(prob; prob_func = mutate_prob!)

esol = solve(ens, AutoTsit5(Rosenbrock23()), EnsembleThreads(); trajectories = trs)

plotint = 500..750

titlestr = Node("j")

nv = [esol[500](LinRange(plotint.left, plotint.right, 1000))[i, :] for i in 1:2:8]

n1, n2, n3, n4 = Node.(nv)

lines_sc = lines(plotint, n1; color = :blue, linewidth = 5)
lines!(lines_sc, plotint, n2; color = :orange, linewidth = 5)
lines!(lines_sc, plotint, n3; color = :green, linewidth = 5)
lines!(lines_sc, plotint, n4; color = :red, linewidth = 5)

pixelchild = campixel(lines_sc) # create a pixel-level child scene
legend_sc = legend!(
                pixelchild,
                lines_sc.plots[end-3:end],
                ["n1", "n2", "n3", "n4"];
                position = (0.85, 0.8), # in normalized figure coordinates
                textsize = 20,
                gap = 25,
                backgroundcolor = colorant"light gray"
            ) # plot to the pixel child, so we can have an inline legend
lines_sc

sc = title(lines_sc, titlestr; textsize = 50)

record(sc, "out.mp4", 1:trs; framerate = 30) do i
    titlestr[] = "j = " * string(round(i/trs; digits = 3))
    nv = [esol[i](LinRange(plotint.left, plotint.right, 1000))[i, :] for i in 1:2:8]
    n1[] = nv[1]
    n2[] = nv[2]
    n3[] = nv[3]
    n4[] = nv[4]
    update_limits!(sc)
    @debug "Recording" progress=i/trs
end

ds = ContinuousDynamicalSystem(CPG, u0, p; t0 = 0.0)

tr = trajectory(ds, 500)

# trimmed_ds = tr[25000:50000] |> Dataset

# rm = RecurrenceMatrix(trimmed_ds, 0.01)

plot(
    plot(sol; vars = 1:2:6, title = "neurons 0 and 1", tspan = (500, 750), labels = ["v0", "d0", "v1", "v2"]),
    plot(sol; vars = 1:4:20, title = "ventral voltages", tspan = (500, 750)),
    plot(sol; vars = 2:4:22, title = "dorsal voltages", tspan = (500, 750)),
    plot(sol; vars = [12, 16, 20, 14, 18, 22], title = "last 3 voltages", tspan = (500, 750));
    legend = :none
)

gauss_40 = ImageFiltering.KernelFactors.IIRGaussian((40,))

effectivesignal = zeros(size(tr))

Plots.plot(tr[:, collect(1:5) .* 4])

for neuron_number in 1:5
    effectivesignal[:, neuron_number] = (
        tr[:, 4*neuron_number]
        .+ abs.(tr[:, 4*neuron_number])
        ) ./ 2 .- (tr[:,4*neuron_number+2] .+ abs.(tr[:,4*neuron_number+2]))/2

    effectivesignal[:, neuron_number] = pyfilter(effectivesignal[:, neuron_number])
end

Plots.plot(effectivesignal[round(Int, end*8/10):end, 1:5])
plot(effectivesignal[round(Int, end*8/10):end, 1])

number_of_joints  = 6   # check #  # actually muscles #
tightness_of_bend = 0.5

# xlim=(-8, 8), ylim=(-8, 8)
# begin outer loop over frames in video
function worm_motion(framenumber)
    # s = 0         to be added
    xpos = 0.0
    ypos = 0.0

    xcoords = zeros(number_of_joints)
    ycoords = zeros(number_of_joints)
    # scoords=[s]    to be added
    # deltatheta # this will be vector of neuronal inputs
    phase=0.0
    theta=0.0

    for j in range(1, stop = number_of_joints)
        deltatheta = effectivesignal[8000+framenumber, j]
        theta += tightness_of_bend * deltatheta
        xpos = xpos + cos(theta)
        ypos = ypos + sin(theta)
        xcoords[j] = xpos
        ycoords[j] = ypos
    end
        # print(i,deltatheta, xpos, ypos)


    # spline fit
    nodes = range(1, stop = number_of_joints) .|> Float64
    csx = CubicSpline(nodes, xcoords)
    csy = CubicSpline(nodes, ycoords)
    worm_coords = range(0, number_of_joints+1.01, step=0.01) # I am incrementing the position in steps of
                                                             # 0.01*(distance between nodes)
    smoothwormx = csx.(worm_coords)
    smoothwormy = csy.(worm_coords)

    return (smoothwormx, smoothwormy)
end

@gif for i in 1:5:1000
    plot(worm_motion(i)...; xlims = (-2, 14), ylims = (-10, 10), lw=10, title = "t = $i")
    @debug "giffing" progress=i/length(1:5:1000)
end

# using PyCall
#
# duPy = py"""
# import math
#
# #constants
# dx0= 1.0
# dy0=-0.51
# vx0= 1
# vy0= -0.49
# dx1= 1.0
# dy1=-0.5
# vx1= 1
# vy1= -0.5
# dx2= 1
# dy2=-0.5
# vx2= 1
# vy2= -0.5
# dx3= 1.0
# dy3=-0.51
# vx3= 1
# vy3= -0.49
# dx4= 1.0
# dy4=-0.5
# vx4= 1
# vy4= -0.5
# dx5= 1
# dy5=-0.5
# vx5= 1
# vy5= -0.5
# # amplitude
# a0= 0
# # frequency
# omega= 0.1
#
# #vector of x and y values
# v=[vx0,vy0,dx0,dy0,vx1,vy1,dx1,dy1,vx2, vy2, dx2, dy2, vx3, vy3, dx3, dy3, vx4, vy4, dx4, dy4, vx5, vy5, dx5, dy5 ]
#
# #tests for times
# #t = np.linspace(0, 20000, 100001)
# # t = np.linspace(0, 2000, 20001)
#
# def FHND(v, t):
#     # calculates the derivatives in FHN with diffusion
#     vx0, vy0, dx0, dy0, vx1, vy1, dx1, dy1, vx2, vy2, dx2, dy2, vx3, vy3, dx3, dy3, vx4, vy4, dx4, dy4, vx5, vy5, dx5, dy5 = v
#     # v denotes vector (x0, y0, x1, y1, x2, y2, x3, y3) passed to the function FHND
#     # the statement x0, y0, x1, y1, x2, y2, x3, y3 = v assigns components in the vector v to the
#     # corresponding components.
#     # Thus, do not use v = x0, y0, x1, y1, x2, y2, x3, y3
#     # this allows for four neurons
#     # for now, just use one-⁠way coupling from neuron 0 to neuron 1
#     # allow for different epsilons
#     e0 = 0.08
#     e1 = 0.08
#     g = 0.8
#     b0 = 0.46
#     b1 = 0.47
#     # diffusion constant
#     Dhead = -0.2
#     Drest= -0.02
#     Dgap=0.03
#     J = 0
#
#
#     #
#     dvdt=[vx0-(vx0**3/3)- vy0+ Dhead*dx0                  + J, e0*(vx0-g*vy0+b0),   # neurons 0 (ventral and dorsal)
#           dx0-(dx0**3/3)- dy0+ Dhead*vx0                  + J, e0*(dx0-g*dy0+b0),   # coupled by 1-way diffusion
#                                                                                     # negative diffusion constant Dhead
#                                                                                     # simulated inhibitory synapse
#
#           vx1-(vx1**3/3)- vy1+ Drest*dx1 + Dgap*(vx0-vx1) + J, e1*(vx1-g*vy1+b1),   # neurons 1 (ventral and dorsal)
#           dx1-(dx1**3/3)- dy1+ Drest*vx1 + Dgap*(dx0-dx1) + J, e1*(dx1-g*dy1+b1),   # driven by neurons 0 through one-way
#                                                                                     # diffusion (selective gap junction)
#                                                                                     # Dgap
#                                                                                     # Drest simulates inhibitary synapse
#                                                                                     # coupling ventral and dorsam neurons 1
#
#           vx2-(vx2**3/3)- vy2+ Drest*dx2 + Dgap*(vx1-vx2) + J, e1*(vx2-g*vy2+b1),   # neurons 2
#           dx2-(dx2**3/3)- dy2+ Drest*vx2 + Dgap*(dx1-dx2) + J, e1*(dx2-g*dy2+b1),
#
#           vx3-(vx3**3/3)- vy3+ Drest*dx3 + Dgap*(vx2-vx3) + J, e1*(vx3-g*vy3+b1),   # neurons 3
#           dx3-(dx3**3/3)- dy3+ Drest*vx3 + Dgap*(dx2-dx3) + J, e1*(dx3-g*dy3+b1),
#
#           vx4-(vx4**3/3)- vy4+ Drest*dx4 + Dgap*(vx3-vx4) + J, e1*(vx4-g*vy4+b1),   # neurons 4
#           dx4-(dx4**3/3)- dy4+ Drest*vx4 + Dgap*(dx3-dx4) + J, e1*(dx4-g*dy4+b1),
#
#           vx5-(vx5**3/3)- vy5+ Drest*dx5 + Dgap*(vx4-vx5) + J, e1*(vx5-g*vy5+b1),   # neurons 5
#           dx5-(dx5**3/3)- dy5+ Drest*vx5 + Dgap*(dx4-dx5) + J, e1*(dx5-g*dy5+b1)]
#
#
#     # diffusion is represented by the term D*(x1-⁠x0)
#     return dvdt
# """
#
# py"FHND"(v, 0.0) |> plot!


# function CPG!(du, u, p, t)
#
#     # unpack the state array
#
#     vxs = @view(u[1:4:end]) # 4-period of states.
#     dxs = @view(u[2:4:end])
#     vys = @view(u[3:2:end])
#     dys = @view(u[4:2:end])
#
#     e0, e1, g, b0, b1, Dhead, Drest, Dgap, J = p # unpack the parameter array
#
#     # Head oscillator
#     # neurons 0 (ventral and dorsal) coupled by 1-way diffusion,
#     # negative diffusion constant `Dhead` simulated inhibitory synapse.
#     du[1] = f(vxs[1]) - vys[1] + Dhead * dxs[1] + J
#     du[2] = e0 * (vxs[1] - g * vys[1] + b0)
#     du[3] = f(dxs[1]) - dys[1] + Dhead * vxs[1] + J
#     du[4] = e0 * (dxs[1] - g * dys[1] + b0)
#
#
#     # body oscillators
#     for i in 2:length(vxs)-1
#         ind = i*4 # account for the periodicity of the array.
#         du[ind]   = f(vxs[i]) - vys[i] + Drest * dxs[i] + Dgap * (vxs[i-1] - vxs[i]) + J
#         du[ind+1] = e1 * (vxs[i] - g * vys[i] + b1)
#         du[ind+2] = f(dxs[i]) - dys[i] + Drest * vxs[i] + Dgap * (dxs[i-1] - dxs[i]) + J
#         du[ind+3] = e1 * (dxs[i] - g * dys[i] + b1)
#     end
#
# end

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

    # body oscillators
    for i in 2:length(vxs)
        ind = (i-1)*4+1 # account for the periodicity of the array.
        du[ind]   = f(vxs[i]) - vys[i] + Drest * dxs[i] + Dgap * (vxs[i-1] - vxs[i]) + J
        du[ind+1] = e1 * (vxs[i] - g * vys[i] + b1)
        du[ind+2] = f(dxs[i]) - dys[i] + Drest * vxs[i] + Dgap * (dxs[i-1] - dxs[i]) + J
        du[ind+3] = e1 * (dxs[i] - g * dys[i] + b1)
    end

    # dx1-(dx1^3/3)- dy1+ Drest*vx1 + Dgap*(dx0-dx1) + J, e1*(dx1-g*dy1+b1)
    any(isnan.(du)) && @infiltrate


end

using DifferentialEquations|

prob = ODEProblem(CPG!, u0, tspan, p)

using Plots

sol = solve(prob)

plot(sol; vars = 1:2:4, tspan = (500, 750), legend = :none)
savefig("lol.png")
