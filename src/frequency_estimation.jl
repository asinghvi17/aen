using Hilbert

"""
    arg(z::Complex)

Computes the "angle" of the given complex number
as a vector away from the real axis, counterclockwise.
"""
arg(z::Complex) = atan(z.im, z.re)

"""
    find_period(data; atol = 0.005, ptol = 5*atol)

Uses the Hilbert transform to determine
the periodicity of the dataset.
This function assumes that the data is known to be
periodic.

Returns the index of the first entry that is 2π away from the initial entry.

Arguments
***

- atol - The error tolerance for the phase difference between indices.
- ptol - The error tolerance for what constitutes "one rotation" of the tangent vector.

!!!note
Currently implemented only for one dimension.
"""
function find_period(data; atol = 0.005, ptol = 5*atol)

    # the function expects a 2d array but we're giving it a 1d timeseries
    H = Hilbert.hilbert(data[:, :])

    ϕ = arg.(H)

    ϕ₀ = ϕ[1]

    Δϕ = 0

    for i in eachindex(ϕ)[2:end]

        Δϕ += abs(ϕ[i] - ϕ[i-1])

        # if the tangent vector has completed a full rotation
        if isapprox(Δϕ, 4π; atol = ptol) # why 4π?

            println(i) # debuggingg info

            println(Δϕ/π) # debugging info

            isapprox(ϕ[i], ϕ₀; atol = atol) && return i # second check - phase angle should be accurate...

        end

    end

    return -1

end

find_period(sin.(0:0.01:8π))
(0:0.01:8π)[629]/π
data = sin.(0:π/10:3π)

H = Hilbert.hilbert(data[:, :]) # the function expects a 2d array but we're giving it a 1d timeseries

ϕ₀ = arg(H[1]) # the initial angle

arg(H[21]) ≈ ϕ₀

using Plots
plotly()
plot(0:π/10:3π, arg.(H[:]))

Δϕ
collect(0:π/10:3π)[21]

for i in eachindex(H)[2:end]
    # if the tangent vector has completed a full rotation
    arg(H[i]) ≈ ϕ₀ && print(i);
end

using DifferentialEquations, ParameterizedFunctions

# create ode
FHND = @ode_def begin

    dv₀ = min(max(-2 - v₀, v₀), 2 - v₀) - w₀ + D*max(v₁ - v₀, 0)

    dw₀ = e₀ * (v₀ - g * w₀ + b₀)

    dv₁ = min(max(-2 - v₁, v₁), 2 - v₁) - w₁ + D*max(v₀ - v₁, 0) + 0.1

    dw₁ = e₁ * (v₁ - g * w₁ + b₁)

end D g e₀ b₀ e₁ b₁]


x0 = -2
y0 =-0.6667
x1 = -2.0
y1 = -0.7

u0 = [x0, y0, x1, y1]

D  = -0.047
g  = 0.8
e0 = 0.04
b0 = 0
e1 = 0.04
b1 = 0.2 # was 0.25

p0 = [D, g, e0, b0, e1, b1]

ts = LinRange(1000, 1250, 10000)
sol = solve(ODEProblem(FHND, u0, (0f0, 2000f0), p0), save_idxs = [1, 3])

# plot the solution


# using Plots
#
# plot(
#     sol,
#     legend=false,
#     tspan = (1000.0, 1500.0)
#     )


## try to save
tr = LinRange(1000, 1250, 10000)


sol_a = transpose(hcat(sol.(tr)...))

# plot(sol_a[:, 1], sol_a[:, 3]);

v1 = sol_a[:, 1]

find_period(v1)

plotlyjs()

p1 = plot(arg.(hilbert(v1[:, :])))

p2 = plot!(v1)

plot!(abs.(hilbert(v1[:, :])))

plot(p1)


H = Hilbert.hilbert(v1[:, :]) # the function expects a 2d array but we're giving it a 1d timeseries

ϕ₀ = arg(H[1]) # the initial angle



for i in eachindex(H)[2:end]
    # if the tangent vector has completed a full rotation
    isapprox(arg(H[i]), ϕ₀, atol = 0.005) && return i;
end

plot(arg.(hilbert((x -> sin(x) + 2*cos(2*x)).(0:0.0001:2π)[:, :])))

plot!(abs.(H))

using NumericalIntegration;
integrate(tr, abs.(H))



## tests

using Test

@testset "period estimation" begin

        @test find_period(sin.(0:0.01:8π)) == 629

        @test find_period((x -> sin(x) + cos(2x)).(0:0.01:8π)) == 629

end

using Plots

plot((x -> sin(x) + cos(2x)).(0:0.01:8π))

data = (x -> sin(x) + cos(2x)).(0:0.01:8π)

find_period(data; atol = 0.5)

plot(arg.(hilbert(data[:, :])))

scatter!(arg.(hilbert(data[:, :])))
