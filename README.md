# Analog Electronic Neurons

This is a project concerning the implementation of biomimetic neurons (of the FitzHugh-Nagumo model) onto analog circuits, with the ultimate goal being to simulate a worm (_C. elegans_).

# Installing dependencies

You will need to have the `Julia` language installed, on version 1.4 or higher.  You can download it at https://julialang.org/downloads/; once you have done that, clone the repo, and run this command in the root directory of the project:

`julia --project=. -e 'using Pkg; Pkg.instantiate()'`

This will install all the Julia package dependencies.

# Project Setup
```
aen          <- Project's main folder.
│
├── _research        <- WIP scripts, code, notes, comments,
│   |                   to-dos and anything in an alpha state.
│   └── tmp          <- Temporary data folder.
│
├── data             <- **Immutable and add-only!**
│   ├── sims         <- Data resulting directly from simulations.
│   ├── exp_pro      <- Data from processing experiments.
│   └── exp_raw      <- Raw experimental data.
│
├── plots            <- Self-explanatory.
├── notebooks        <- Jupyter, Weave or any other mixed media notebooks.
│
├── papers           <- Scientific papers resulting from the project.
│
├── scripts          <- Various scripts, e.g. simulations, plotting, analysis,
│   │                   The scripts use the `src` folder for their base code.
│   └── intro.jl     <- Simple file that uses DrWatson and uses its greeting.
│
├── src              <- Source code for use in this project. Contains functions,
│                       structures and modules that are used throughout
│                       the project and in multiple scripts.
│
├── test             <- Tests for code in src.
│   └── runtests.jl  <- Script to run all tests.
│
├── README.md        <- You are here.
├── .gitignore       <- ignores latex-compilation related files and filesystem artifacts.
│
├── Manifest.toml    <- Contains full list of exact Julia package versions used currently.
└── Project.toml     <- Main Julia project file, allows activation and installation.
                        Includes DrWatson by default.
```

# Making changes

If you want to make a change to a file, please branch out from `master` and submit a pull request.  Since this is a shared repo, the recommended naming scheme or a branch is `<initials>/<short-description>`.  An example is `as/depwarns`.
