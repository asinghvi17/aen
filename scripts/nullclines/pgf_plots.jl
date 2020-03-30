using PGFPlotsX

push!(
    PGFPlotsX.CUSTOM_PREAMBLE,
    raw"""
    % This cycle list encodes the Wong colors, which are visually distinguishable.
    % They also account for colorblind support, and as such are optimal for use in a paper.
    \pgfplotscreateplotcyclelist{wong}{%
        color={rgb, 255 : red, 86  ; green, 180 ; blue, 233}, mark = *\\         % sky blue
        color={rgb, 255 : red, 230 ; green, 159 ; blue, 0},   mark = square*\\   % orange
        color={rgb, 255 : red, 0   ; green, 158 ; blue, 115}, mark = otimes*\\   % blueish green
        color={rgb, 255 : red, 240 ; green, 228 ; blue, 66},  mark = star\\      % yellow
        color={rgb, 255 : red, 0   ; green, 114 ; blue, 178}, mark = diamond*\\  % blue
        color={rgb, 255 : red, 213 ; green, 94  ; blue, 0},   mark = triangle*\\ % vermillion
        color={rgb, 255 : red, 204 ; green, 121 ; blue, 167}, mark = pentagon*\\ % reddish purple
    }
    """
)


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

# Allow Point2f0 to function as a coordinate.  This simply
# decomposes into x and y.
PGFPlotsX.Coordinates(vec::Vector{Point2f0}; kwargs...) = Coordinates(first.(vec), last.(vec); kwargs...)

fig_fh = @pgf GroupPlot(
        {
            group_style = {
                group_size = "2 by 1",
            },
            no_markers,
            xtick  = "\\empty",
            ytick  = "\\empty",
            xlabel = "{\$v\$}",
            ylabel = "{\$w\$}",
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
pgfsave(joinpath(plotsdir(), "Dynamics", "Nullclines", "FHN_nullclines", "fhn_dynamics.tex"), fig_fh; include_preamble = false)
pgfsave(joinpath(plotsdir(), "Dynamics", "Nullclines", "FHN nullclines", "fhn_dynamics.pdf"), fig_fh)|


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
            xtick  = "\\empty",
            ytick  = "\\empty",
            xlabel = "{\$v\$}",
            ylabel = "{\$w\$}",
        },
    {
        title = "Experimental values (Xu et al)",
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

################################################################################
#                         Analog-simulation comparison                         #
################################################################################

########################################
#            Simulated data            #
########################################

sim_data_single = readdlm(datadir("Oscilloscopy", "Comparison data", "Simulation", "single.dat"); skipstart = 1)

sim_single_v0 = Point2f0.(sim_data_single[:, 1], sim_data_single[:, 2])
sim_single_v1 = Point2f0.(sim_data_single[:, 1], sim_data_single[:, 3])

sim_data_double = readdlm(datadir("Oscilloscopy", "Comparison data", "Simulation", "double.dat"); skipstart = 1)

sim_double_v0 = Point2f0.(sim_data_double[:, 1], sim_data_double[:, 2])
sim_double_v1 = Point2f0.(sim_data_double[:, 1], sim_data_double[:, 3])

########################################
#             Analog data              #
########################################

anal_single_v0 = Point2f0.(eachcol(readdlm(datadir("Oscilloscopy", "Comparison data", "Analog", "ch1_s.dat"); skipstart = 1))...)
anal_single_v1 = Point2f0.(eachcol(readdlm(datadir("Oscilloscopy", "Comparison data", "Analog", "ch2_s.dat"); skipstart = 1))...)

anal_double_v0 = Point2f0.(eachcol(readdlm(datadir("Oscilloscopy", "Comparison data", "Analog", "ch1_d.dat"); skipstart = 1))...)
anal_double_v1 = Point2f0.(eachcol(readdlm(datadir("Oscilloscopy", "Comparison data", "Analog", "ch2_d.dat"); skipstart = 1))...)

########################################
#                Figure                #
########################################

fig_cas = @pgf GroupPlot(
        {
            group_style = {
                group_size = "3 by 2",
                horizontal_sep = "2cm",
                vertical_sep   = "1cm",
            },
            no_markers,
            xtick  = "\\empty",
            ytick  = "\\empty",
            xlabel = "{Time (arbitrary)}",
            ylabel = "{\$v\$ (theoretical)}",
        },
    {
        title = "FitzHugh-Nagumo Isoclines",
        domain = "-2:2",
        xlabel = "{\$v\$}",
        ylabel = "{\$w\$}",
    },
    PlotInc({very_thick}, Expression("x - x^3/3")), # FHN v-nullcline
    PlotInc({very_thick}, Expression("x")), # common w-nullcline

    {
        title = "Head oscillator (negative coupling)",
        ytick  = "",
    },
    PlotInc({very_thick}, Coordinates(sim_double_v0)),
    PlotInc({very_thick}, Coordinates(sim_double_v1)),

    {
        title = "Descending pathway (positive coupling)",
        ytick  = "",
    },
    PlotInc({very_thick}, Coordinates(sim_single_v0)),
    PlotInc({very_thick}, Coordinates(sim_single_v1)),

    {
        title = "Keener Isoclines",
        domain = "-2:2",
        xlabel = "{\$v\$}",
        ylabel = "{\$w\$}",
    },
    PlotInc({very_thick}, Coordinates(Point2f0[ # manually expand linear interp
        (-2, -2 - (-2)^3/3),
        (-1, -1 - (-1)^3/3),
        (1, 1 - 1^3/3),
        (2, 2 - (2)^3/3),
    ])), # Keener v-nullcline
    PlotInc({very_thick}, Expression("x")), # common w-nullcline

    {
        xlabel = "Time (s)",
        ylabel = "{\$V\$ (real)}",
        xtick  = "",
        ytick  = "",
    },
    PlotInc({very_thick}, Coordinates(anal_double_v0)),
    PlotInc({very_thick}, Coordinates(anal_double_v1)),

    {
        xlabel = "Time (s)",
        ylabel = "{\$V\$ (real)}",
        xtick  = "",
        ytick  = "",
    },
    PlotInc({very_thick}, Coordinates(anal_single_v0)),
    PlotInc({very_thick}, Coordinates(anal_single_v1)),
)

mkpath(joinpath(plotsdir(), "Dynamics", "Nullclines", "anal_sim_comp"))
pgfsave(joinpath(plotsdir(), "Dynamics", "Nullclines", "anal_sim_comp", "anal_sim_comp.tex"), fig_nc)
pgfsave(joinpath(plotsdir(), "Dynamics", "Nullclines", "anal_sim_comp", "anal_sim_comp.pdf"), fig_nc)
