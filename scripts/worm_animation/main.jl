using DifferentialEquations, ModelingToolkit
using ImageFiltering, Interpolations
using Makie, Observables
using Makie.AbstractPlotting.MakieLayout
using CoordinateTransformations, Rotations
using FFTW, Colors
using CairoMakie
import Plots
Makie.GLMakie.activate!()

function warmup_plotting()
    println("Plots")
    @time Plots.plot(rand(10))
    println("Makie")
    @time scene, layout = layoutscene()
    @time la = layout[1, 1] = LAxis(scene)
    @time heatmap!(la, rand(10, 10))
end
warmup_plotting()

include("cpg.jl")
include("utils.jl")
include("animation.jl")
