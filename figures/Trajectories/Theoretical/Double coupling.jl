using DifferentialEquations,
      DelimitedFiles,
      Plots

#constants
x0= -2
y0=-0.6667
x1= -1.0
y1= 0.7
x2=2.01
y2=-0.5
x3=-0.5
y3=0.5

J = 0.1

#vector of x and y values
u0=[x0,y0,x1,y1]

#tests for times
#t = np.linspace(0, 20000, 100001)
ts = LinRange(0, 100, 5001)

function FHND!(du, u, p, t)
    # calculates the derivatives in FHN with diffusion
    x0, y0, x1, y1 = u
    # v denotes vector (x0, y0, x1, y1, x2, y2, x3, y3) passed to the function FHND
    # the statement x0, y0, x1, y1, x2, y2, x3, y3 = v assigns components in the vector v to the
    # corresponding components.
    # Thus, do not use v = x0, y0, x1, y1, x2, y2, x3, y3
    # this allows for four neurons
    # for now, just use one-⁠way coupling from neuron 0 to neuron 1
    # allow for different epsilons
    e0 = 0.04
    e1 = 0.04
    e2 = 0.04
    e3 = 0.04
    g  = 0.8
    b0 = 0.42
    b1 = 0.46
     # diffusion constant
    D  = 0.03

    # here we will have coupling from neuron 0 to neuron 1
    du[1] = x0-(x0^3 / 3)- y0 + J # + D * max((x0 - x1),0)
    du[2] = e0 * (x0 - g * y0 + b0)
    du[3] = x1-(x1^3 / 3)- y1 + D * max((x1 - x0),0) + J
    du[4] = e1 * (x1 - g * y1 + b1)
    end


tr = (0.0, 2500.0)

ts = LinRange(200.0, 350.0, 500)

prob = ODEProblem(FHND!, u0, tr)

sol = solve(prob, dense=true, save_idxs = [1, 3])
#
# using Plots
#
## different approach




x0= -2
y0=-0.6667
x1= -2.0
y1= -0.7

u0 = [x0, y0, x1, y1]

function FHND(v, p, t)
  # calculates the derivatives in FHN with diffusion
  x0, y0, x1, y1 = v
  e0 = 0.04
  e1 = 0.04
  e2 = 0.04
  e3 = 0.04
  g = 0.8
  b0 = 0
  b1 = 0.2 # was 0.25
   # diffusion constant
  D = -0.047

  # here we will have coupling from neuron 0 to neuron 1
  dvdt=[
        min(max(-2 - x0, x0), 2 - x0) - y0 + D*max(x1 - x0, 0),
        e0 * (x0 - g * y0 + b0),
        min(max(-2 - x1, x1), 2 - x1) - y1 + D*max(x0 - x1, 0) + 0.1,
        e1 * (x1 - g * y1 + b1)
    ]
  # dvdt=[x0-(x0**3/3)- y0 , e0*(x0-g*y0+b0),
  # x1-(x1**3/3)- y1+ D*max((x1-x0),0) , e1*(x1-g*y1+b1)]


  # diffusion is represented by the term D*(x1-⁠x0)
  return dvdt
end

ts = LinRange(1000, 1250, 10000)
sol = solve(ODEProblem(FHND, u0, (0f0, 2000f0)), save_idxs = [1, 3])

plot(
    sol,
    legend=false,
    tspan = (1000.0, 1500.0)
    )


tr = LinRange(1000, 1250, 1000)


sol_a = hcat(tr, transpose(hcat(sol.(tr)...)))

plot(sol_a[:, 1], sol_a[:, 3]);
plot!(sol_a[:, 1], sol_a[:, 2])

sol_arr = sol_a |> x -> vcat(["t" "v0" "v1"], x)



writedlm("./double.dat", sol_arr)
#
#
#
# sol_a = hcat(ts, transpose(hcat(sol.(ts)...)))
#
# plot(sol_a[:, 1], sol_a[:, 2]);
# plot!(sol_a[:, 1], sol_a[:, 3])
#
# sol_arr = sol_a |> x -> vcat(["t" "v1" "v0"], x)
#
#
#
# writedlm("./double.dat", sol_arr)

using FFTW

fft(sol_a[:, 3])

scatter(fft(sol_a[:, 3]), marker_z = abs.(fft(sol_a[:, 3])))
