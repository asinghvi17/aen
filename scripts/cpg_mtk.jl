using ModelingToolkit, DifferentialEquations
using Makie, MakieLayout

# ## Initial definitions

# First, we define the parameters of the atomic
# FHN system.  Time is treated as a "parameter"
# for ModelingToolkit, and g, e and b are the standard
# FHN params.
@parameters t g e b

# Now, we define variables.  $v(t)$ and $w(t)$ are
# again the standard FHN variables, and $F(t)$ is a
# general "flux" term which is used as a sink for
# all of the diffusion terms.
@variables v(t) w(t) F(t)

# This block is pretty simple - it defines a
# derivative operator $D$ as the time derivative.
@derivatives D'~t

# ## Atomic system definition

# We now define the atomic FHN equations.
# These have the standard structure, and on their
# own would simulate the FHN as usual.

fhn = [
    D(v) ~ min(max(-2 - v, v), 2 - v) - w + F,
    D(w) ~ e * (v - g * w + b)
]

# Here, we create a vector of systems.  This is just
# a collection of systems which doesn't encode any crucial meaning - yet.

# First
Npair = 6

systems = [
    ODESystem(
        fhn,
        t,
        [v,w,F],
        [g,e,b];
        name = Symbol(:n, i)
    )
    for i in 1:2*Npair
]

ventral = systems[1:2:end]
dorsal  = systems[2:2:end]

# ## Connections

@parameters Dhead J Hk[1:2]

gaps, inhibs = @parameters Dgap[1:(Npair-1)] Dinhib[1:(Npair-1)]

head_connections = [
    0 ~ ventral[1].F + Dhead * Hk[1] * max(ventral[2].v - ventral[1].v)
    0 ~ ventral[2].F + Dhead * Hk[2] * max()
]

# First, we handle the head oscillator.  This only has inhibitory coupling between the ventral and dorsal neurons, so it's not that big of a deal.

# ## Bringing it all together

# ## Solving the equation
