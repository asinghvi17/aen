include("intro.jl")

using DifferentialEquations, ParameterizedFunctions

# create ode
FHND = @ode_def begin

    dv₀ = min(max(-2 - v₀, v₀), 2 - v₀) - w₀ + D*max(v₁ - v₀, 0)

    dw₀ = e₀ * (v₀ - g * w₀ + b₀)

    dv₁ = min(max(-2 - v₁, v₁), 2 - v₁) - w₁ + D*max(v₀ - v₁, 0) + 0.1

    dw₁ = e₁ * (v₁ - g * w₁ + b₁)

end D g e₀ b₀ e₁ b₁


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

using DynamicalSystems
T0 = 1000.0
T = 1350.0
ds = ContinuousDynamicalSystem(FHND, u0, p0; t0 = 1000.0)
dsol = trajectory(ds, T)
tr = LinRange(T0, T, size(dsol, 1))

using LombScargle

lsp = LombScargle.plan(tr, dsol[:, 1])
ls = lombscargle(lsp)
findmaxperiod(ls)

using Plots

timeseries = plot(
                tr,
                dsol[:, 1],
                xlabel = L"t",
                ylabel = L"V",
                title = "Timeseries",
                legend = false
            )

periodogram = plot(
    periodpower(ls)...,
    xlims = (0, 50),
    legend=false,
    xlabel = L"T",
    ylabel = L"P",
    title = "Periodogram"
    )

savefig(plot(timeseries, periodogram, layout = @layout [a; b]), "ololol.pdf")

findmaxperiod(ls, [1/(T - T0), Inf])
