using DifferentialEquations,
      DelimitedFiles,
      Plots

x0= -2
y0=-0.6667
x1= -1.0
y1= 0.7

u0 = [x0, y0, x1, y1]

function FHND(v, p, t)
  # calculates the derivatives in FHN with diffusion
  x0, y0, x1, y1 = v
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
  g = 0.8
  b0 = 0
  b1 = 0.2 # was 0.25
   # diffusion constant
  D = 0.003

  # here we will have coupling from neuron 0 to neuron 1
  dvdt=[min(max(-2-x0,x0), 2-x0)- y0 , e0*(x0-g*y0+b0),
        min(max(-2-x1,x1), 2-x1)- y1+ D*max((x0-x1),0) + 0.2, e1*(x1-g*y1+b1)]
  # dvdt=[x0-(x0**3/3)- y0 , e0*(x0-g*y0+b0),
  # x1-(x1**3/3)- y1+ D*max((x1-x0),0) , e1*(x1-g*y1+b1)]


  # diffusion is represented by the term D*(x1-⁠x0)
  return dvdt
end

ts = LinRange(1000, 1250, 10000)
sol = solve(ODEProblem(FHND, u0, (0f0, 2000f0)), save_idxs = [1, 3])

plot(sol, legend=false, tspan = (1000.0, 1180.0))


tr = LinRange(1000, 1250, 1000)


sol_a = hcat(tr, transpose(hcat(sol.(tr)...)))

plot(sol_a[:, 1], sol_a[:, 3]);
plot!(sol_a[:, 1], sol_a[:, 2])

sol_arr = sol_a |> x -> vcat(["t" "v0" "v1"], x)



writedlm("./single.dat", sol_arr)
