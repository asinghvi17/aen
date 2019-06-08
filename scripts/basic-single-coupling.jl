include("intro.jl")

# load packages
using DifferentialEquations,
      ParameterizedFunctions

# create ode
FHNS = @ode_def begin

    dv₀ = min(max(-2 - v₀, v₀), 2 - v₀) - w₀

    dw₀ = e₀ * (v₀ - g * w₀ + b₀)

    dv₁ = min(max(-2 - v₁, v₁), 2 - v₁) - w₁ + D*max(v₀ - v₁, 0) + 0.1

    dw₁ = e₁ * (v₁ - g * w₁ + b₁)

end D g e₀ b₀ e₁ b₁


x0 = -2
y0 =-0.6667
x1 = -2.0
y1 = -0.7

u0 = [x0, y0, x1, y1]

D  = 0.003
g  = 0.8
e0 = 0.04
b0 = 0
e1 = 0.04
b1 = 0.2 # was 0.25

p0 = [D, g, e0, b0, e1, b1]

ts = LinRange(1000, 1250, 10000)
sol = solve(ODEProblem(FHNS, u0, (0f0, 2000f0), p0), save_idxs = [1, 3])

# plot the solution


# using Plots
#
# plot(
#     sol,
#     legend=false,
#     tspan = (1000.0, 1500.0)
#     )


## try to save
tr = LinRange(1000, 1250, 1000)


sol_a = transpose(hcat(sol.(tr)...))

# plot(sol_a[:, 1], sol_a[:, 3]);
# plot!(sol_a[:, 1], sol_a[:, 2])
#
# sol_arr = sol_a |> x -> vcat(["t" "v0" "v1"], x)

# writedlm("$(datadir())/double_plot.dat", sol_arr)

using PGFPlotsX

plt = @pgf Axis(
        {
        legend_pos = "outer north east",
        xlabel = raw"{$t$ (arbitrary)}",
        ylabel = raw"$V$",
        xtick  = raw"\empty",
        title = "Head oscillator (positive coupling)"
        # grid=major,
        },

        PlotInc(
            {
            no_marks,
            very_thick,
            blue
            },
            Coordinates(tr, sol_a[:, 1])
        ),

        LegendEntry(L"V_0"),

        PlotInc(
            {
            no_marks,
            very_thick,
            red
            },
            Coordinates(tr, sol_a[:, 2])
        ),

        LegendEntry(L"V_1"),
    )

# save the plot
pgfsave(plotsdir() * "sim-single/sim-single.tex", plt)
pgfsave(plotsdir() * "sim-single/sim-single.pdf", plt)
pgfsave(plotsdir() * "sim-single/sim-single.png", plt)

plt1 = @pgf Axis(
        {
        legend_pos = "outer north east",
        xlabel = raw"{$t$ (arbitrary)}",
        ylabel = raw"$V$",
        xtick  = raw"\empty",
        title = "Non-sinusoidal oscillator"
        # grid=major,
        },

        PlotInc(
            {
            no_marks,
            very_thick,
            blue
            },
            Coordinates(tr, sol_a[:, 1])
        ),

        # LegendEntry(L"V_0")

        # PlotInc(
        #     {
        #     no_marks,
        #     very_thick,
        #     red
        #     },
        #     Coordinates(tr, sol_a[:, 2])
        # ),
        #
        # LegendEntry(L"V_1"),
    )

pgfsave("giordano.png", plt1)
