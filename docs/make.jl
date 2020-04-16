using Pkg

Pkg.activate(dirname(@__DIR__)) # activate the main project
using Literate
# using DifferentialEquations, ModelingToolkit

# Set up the index / home page
cp("$(dirname(@__DIR__))/README.md", "$(@__DIR__)/src/index.md"; force = true)

# Generate the Literate source files

scriptdir = joinpath(dirname(@__DIR__), "scripts")

gendir = joinpath(@__DIR__, "src", "generated")

exe_files = [
    "neural_sync.jl",
    "cpg_sync.jl",
]

noexe_files = [
    "neural_sync_mtk.jl",
    "cpg_sync_mtk.jl"
]

Literate.markdown.(
    joinpath.(scriptdir, exe_files),
    gendir;
    documenter = false,
    execute = false
)

makedocs(...)
