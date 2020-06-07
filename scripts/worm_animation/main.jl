using DifferentialEquations, ModelingToolkit
using ImageFiltering, Interpolations
using Makie, MakieLayout, Observables
using CoordinateTransformations, Rotations
using FFTW, Colors
using CairoMakie
import Plots
Makie.GLMakie.activate!()
include("cpg.jl")
include("utils.jl")
include("animation.jl")
