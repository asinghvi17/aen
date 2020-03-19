using DrWatson

@quickactivate

DrWatson.greet()

DrWatson.plotsdir() = joinpath(dirname(@__DIR__), "figures")
