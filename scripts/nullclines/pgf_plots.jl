using PGFPlotsX

################################################################################
#                           Fitzhugh-Nagumo dynamics                           #
################################################################################

fh_exc = (x...) -> fitzhugh_nagumo(x...; i_ext = 0.6)
fh_osc = (x...) -> fitzhugh_nagumo(x...)

iso_vosc = isoclines(-2..2, -2..2, fh_osc; ind=1)
iso_wosc = isoclines(-2..2, -2..2, fh_osc; ind=2)

iso_vexc = isoclines(-2..2, -2..2, fh_exc; ind=1)
iso_wexc = isoclines(-2..2, -2..2, fh_exc; ind=2)

lines(iso_vosc)
lines!(iso_wosc)

lines(iso_vexc)
lines!(iso_wexc)

PGFPlotsX.Coordinates(vec::Vector{Point2f0}; kwargs...) = Coordinates(first.(vec), last.(vec); kwargs...)

fig_fh = @pgf GroupPlot(
        {
            group_style = {
                group_size = "2 by 1",
            },
            no_markers,
        },
    {
        title = "Oscillatory regime",
        legend_entries = raw"{{$v$-nullcline}, {$w$-nullcline}}",
        legend_style = {
            at = "{(1.1, -0.1)}",
            anchor = "north"
        },
    },
    PlotInc({very_thick}, Coordinates(iso_vosc)),
    PlotInc({very_thick}, Coordinates(iso_wosc)),

    {title = {"Excitable regime"}},
    PlotInc({very_thick}, Coordinates(iso_vexc)),
    PlotInc({very_thick}, Coordinates(iso_wexc)),

)

mkpath(joinpath(plotsdir(), "Dynamics", "Nullclines", "FHN nullclines"))
pgfsave(joinpath(plotsdir(), "Dynamics", "Nullclines", "FHN nullclines", "fhn_dynamics.tex"), fig_fh)
pgfsave(joinpath(plotsdir(), "Dynamics", "Nullclines", "FHN nullclines", "fhn_dynamics.pdf"), fig_fh)


################################################################################
#                           Neuron model nullclines                            #
################################################################################

xs = Base.LinRange(-60..0, 500)

xu_v_null = Point2f0.(xs, xu_vnull.(xs))
xu_w_null = Point2f0.(xs, xu_wnull.(xs))

ml_v_null = isoclines(-100..100, -1..2, morris_lecar; ind = 1)
ml_w_null = isoclines(-100..100, -1..2, morris_lecar; ind = 2)

fh_v_null = iso_vosc
fh_w_null = iso_wosc

fig_nc = @pgf GroupPlot(
        {
            group_style = {
                group_size = "3 by 1",
            },
            no_markers,
        },
    {
        title = "Real \\textit{C. elegans}",
    },
    PlotInc({very_thick}, Coordinates(xu_v_null)),
    PlotInc({very_thick}, Coordinates(xu_w_null)),

    {
        title = "Morris-Lecar",
        legend_entries = raw"{{$v$-nullcline}, {$w$-nullcline}}",
        legend_style = {
            at = "{(0.5, -0.1)}",
            anchor = "north"
        },
        xmin = -75,
        xmax = 75,
        ymin = -0.05,
        ymax = 0.6,
    },
    PlotInc({very_thick}, Coordinates(ml_v_null)),
    PlotInc({very_thick}, Coordinates(ml_w_null)),

    {
        title = "FitzHugh-Nagumo",
    },
    PlotInc({very_thick}, Coordinates(fh_v_null)),
    PlotInc({very_thick}, Coordinates(fh_w_null)),
)

mkpath(joinpath(plotsdir(), "Dynamics", "Nullclines", "Neuron model nullclines"))
pgfsave(joinpath(plotsdir(), "Dynamics", "Nullclines", "Neuron model nullclines", "nm_dynamics.tex"), fig_nc)
pgfsave(joinpath(plotsdir(), "Dynamics", "Nullclines", "Neuron model nullclines", "nm_dynamics.pdf"), fig_nc)
