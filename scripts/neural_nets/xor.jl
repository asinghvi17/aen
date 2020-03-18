using Plots
using Flux, DiffEqFlux
using RDatasets

data = dataset("datasets", "iris")

xs = rand(0:1 , 20)
ys = rand(0:1 , 20)
zs = xor.(xs, ys)

loss(x, y) = Float32(xor(round(Int, x), round(Int, y)))

model = Chain(Dense(2, 2), Dense(2, 1))

optimizer = Descent()

params = Flux.params(model)

Flux.train!(model, params, zip(Float32.(xs), Float32.(ys)), optimizer)
