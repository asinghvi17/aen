cd("/Users/anshul/cw/Spring 2020/aen/notebooks")

import Pkg
Pkg.activate(".")
Pkg.pkg"add Plots, DifferentialEquations, StaticArrays, PyCall, ImageFiltering, DataInterpolations"

using Plots, DifferentialEquations, StaticArrays, PyCall, ImageFiltering, DataInterpolations

function f(v)
    return min(max(-2 - v, v), 2 - v) # linear activation function
end

function FHN(u, p, t)
    Dhead, Drest, Dgap, vent_diff, e0, b0, e1, b1, g, J = p # unpack parameters
    vx0, vy0, dx0, dy0, vx1, vy1, dx1, dy1, vx2, vy2, dx2, dy2, vx3, vy3, dx3, dy3, vx4, vy4, dx4, dy4, vx5, vy5, dx5, dy5 = u # unpack states

    return SVector(
        f(vx0)- vy0+ Dhead*vent_diff*dx0 + J, e0*(vx0-g*vy0+b0),
        # neurons 0 (ventral and dorsal)
        f(dx0)- dy0+ Dhead*vx0 + J, e0*(dx0-g*dy0+b0),
        # coupled by 1-way diffusion
        # negative diffusion constant Dhead
        # simulated inhibitory synapse

        # neurons 1 (ventral and dorsal)
        # driven by neurons 0 through one-way
        # diffusion (selective gap junction) Dgap.
        # Drest simulates inhibitory synapse
        # coupling ventral and dorsam neurons 1
        f(vx1)- vy1+ Drest*dx1 + Dgap*(vx0-vx1) + J, e1*(vx1-g*vy1+b1),
        f(dx1)- dy1+ Drest*vx1 + Dgap*(dx0-dx1) + J, e1*(dx1-g*dy1+b1),


        f(vx2)- vy2+ Drest*dx2 + Dgap*(vx1-vx2) + J, e1*(vx2-g*vy2+b1),   # neurons 2
        f(dx2)- dy2+ Drest*vx2 + Dgap*(dx1-dx2) + J, e1*(dx2-g*dy2+b1),

        f(vx3)- vy3+ Drest*dx3 + Dgap*(vx2-vx3) + J, e1*(vx3-g*vy3+b1),   # neurons 3
        f(dx3)- dy3+ Drest*vx3 + Dgap*(dx2-dx3) + J, e1*(dx3-g*dy3+b1),

        f(vx4)- vy4+ Drest*dx4 + Dgap*(vx3-vx4) + J, e1*(vx4-g*vy4+b1),   # neurons 4
        f(dx4)- dy4+ Drest*vx4 + Dgap*(dx3-dx4) + J, e1*(dx4-g*dy4+b1),

        f(vx5)- vy5+ Drest*dx5 + Dgap*(vx4-vx5) + J, e1*(vx5-g*vy5+b1),   # neurons 5
        f(dx5)- dy5+ Drest*vx5 + Dgap*(dx4-dx5) + J, e1*(dx5-g*dy5+b1)
    )
end













# setup

# initial conditions
u0 = Float64[
    1,     # vx0
    -0.49, # vy0
    1.0,   # dx0
    -0.51, # dy0
    1,     # vx1
    -0.5,  # vy1
    1.0,   # dx1
    -0.5,  # dy1
    1,     # vx2
    -0.5,  # vy2
    1,     # dx2
    -0.5,  # dy2
    1,     # vx3
    -0.49, # vy3
    1.0,   # dx3
    -0.51, # dy3
    1,     # vx4
    -0.5,  # vy4
    1.0,   # dx4
    -0.5,  # dy4
    1,     # vx5
    -0.5,  # vy5
    1,     # dx5
    -0.5   # dy5
]

p = Float64[
    -0.2,  # Dhead,
    -0.02, # Drest,
    0.03,  # Dgap,
    0.0,   # vent_diff,
    0.08,  # e0,
    0.46,  # b0,
    0.08,  # e1,
    0.47,  # b1,
    0.8,   # g,
    0.5,   # J
]

prob = ODEProblem(FHN, u0, (0.0, 2000.0), p; dtmax = 0.1)

sol = solve(prob)

pyfhn(v, t) = FHN(v, p, t)

scipy = PyCall.pyimport(:scipy)

function eff_signal(
        sol;
        filter = ImageFiltering.reflect(ImageFiltering.Kernel.gaussian((40,))) # convolution, not correlation
    )
    signal = zeros(size(sol, 2), size(sol, 1))
    for neuron_number in 1:6
        signal[:, neuron_number] = (sol[4*neuron_number, :] + abs.(sol[4*neuron_number, :]))/2 - (sol[4*neuron_number+2, :] + abs.(sol[4*neuron_number+2, :]))/2
        signal[:, neuron_number] = imfilter(signal[:, neuron_number], filter)
        # Gaussian filter simulates a combination of muscle reponse to neural
        # stimulus; also effects of elastic properties of worm and interaction with
        # fluid
    end
    return signal
end

effective_signal = eff_signal(sol)
effective_signal |> size

number_of_bends = 6
# Increment the position in steps of
# 0.01*(distance between nodes)
worm_coords=range(0, stop = number_of_bends+1.01, step = 0.01)
sthetas = []


for framenumber in 0:5:400
    s = 0     # arclength
    theta = 0
    svals = Float64[s]
    thetavals= Float64[theta]
    for j in 1:number_of_bends
        s += 1
        theta+= effective_signal[8000+framenumber,j]
        append!(svals, s)
        append!(thetavals, theta)
    end
    @debug(svals)
    @debug(thetavals)
    # spline fit
    cstheta=DataInterpolations.CubicSpline(thetavals, svals)


    push!(sthetas, cstheta.(worm_coords))
end

plot(worm_coords, sthetas[1]; aspect_ratio = 1)

g = @gif for stheta in sthetas
    plot(worm_coords, stheta; xlims = (0, 8), ylims = (-8, 8), aspect_ratio = 1)
end
