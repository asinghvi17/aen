AbstractPlotting.set_theme!(
    font = "CMU Serif Roman", # same font as the document.  Can be changed arbitrarily.
    fontsize = 12,            # 12-point font.
    LAxis = (
        xticklabelsize = 10,
        yticklabelsize = 10,
        xtickalign = 1,
        ytickalign = 1,
        xticksize = 3,
        xtickwidth = .5,
        yticksize = 3,
        ytickwidth = .5,
        xlabelpadding = 1,
        ylabelpadding = 1,
        titlesize = 10,
        titlegap = 2,
        xgridwidth = 0.5,
        ygridwidth = 0.5,
        spinewidth = 0.5
    )
)

const a4_resolution = (595, 842)

datadir(args...) = joinpath("data", args...)

################################################################################
#                           Fitzhugh-Nagumo dynamics                           #
################################################################################

fh_exc = (x...) -> fitzhugh_nagumo(x...; i_ext = 0.9)
fh_osc = (x...) -> fitzhugh_nagumo(x...)

iso_vosc = isoclines(-2..2, -2..2, fh_osc; ind=1)
iso_wosc = isoclines(-2..2, -2..2, fh_osc; ind=2)

iso_vexc = isoclines(-2..2, -2..2, fh_exc; ind=1)
iso_wexc = isoclines(-2..2, -2..2, fh_exc; ind=2)

scene, layout = Makie.AbstractPlotting.layoutscene(
    0;                      # padding
    resolution = (400, 300) # resolution in points
)

axs = layout[1, 1:2] = [
            LAxis(scene;
                aspect = AxisAspect(1),
                xlabel = "𝑣",
                ylabel = "𝑤",
                xlabelpadding = 1,
                ylabelpadding = 1,
            )
            for i in 1:2
        ]

hidexdecorations!.(axs); hideydecorations!.(axs)
setproperty!.(axs, :xlabelvisible, true)
setproperty!.(axs, :ylabelvisible, true)
# setproperty!.(axs, :spinewidth, true)
# setproperty!.(axs, :yspinewidth, true)
axs[2].yaxisposition = :right

# setproperty!.(axs, :titlesize, 12)
# `axs[1]` is the excitatory axis.
axs[1].title = "Excitatory mode"

vexc_lines = lines!(axs[1], iso_vexc; linewidth = 1, color = wong[1])
wexc_lines = lines!(axs[1], iso_wexc; linewidth = 1, color = wong[2])

# `axs[2]` is the oscillatory axis.
axs[2].title = "Oscillatory mode"

vosc_lines = lines!(axs[2], iso_vosc; linewidth = 1, color = wong[1])
wosc_lines = lines!(axs[2], iso_wosc; linewidth = 1, color = wong[2])

# Now, we create the legend.

leg = layout[2, 1:2] = LLegend(scene,
                [vexc_lines, wexc_lines],
                [" ", " "];
)

# The legend doesn't need much space...80 points should be enough.
rowsize!(layout, 2, 80)

# Now comes the shady hacking section.  We add layouted TeX to the legend's internal gridlayout.

leg.grid[1, 2] = LTeX(leg.child,
    raw"dv = \frac{v^3}{3} - v + w + I_{ext} = 0"
)

leg.grid[2, 2] = LTeX(leg.child,
    raw"dw = e \cdot \left(v- g \cdot w + b\right) = 0"
)

save("neur_dyn_tmp.pdf", scene)
save("neur_dyn_tmp.png", scene; px_per_unit = 10)

# save(plotsdir("Dynamics/Nullclines/FHN_nullclines/fhn_dynamics_makie.pdf"), scene)

################################################################################
#                           Neuron model nullclines                            #
################################################################################


xs = Base.LinRange(-60..0, 500)

xu_v_null = Point2f0.(xs, xu_vnull.(xs))
xu_w_null = Point2f0.(xs, xu_wnull.(xs))

fh_v_null = iso_vosc
fh_w_null = iso_wosc

scene, layout = layoutscene(
    0;                      # padding
    resolution = (250, 150) # resolution in points
)

axs = layout[1, 1:2] = [
            LAxis(scene;
                # aspect = AxisAspect(1),
                xlabel = "𝑣",
                ylabel = "𝑤",
                xlabelpadding = 1,
                ylabelpadding = 1,
                titlesize = 11,
                titlegap = 2,
                xgridwidth = 0.5,
                ygridwidth = 0.5,
                spinewidth = 0.5
                # autolimitaspect = 1
            )
            for i in 1:2
        ]

hidexdecorations!.(axs); hideydecorations!.(axs)
setproperty!.(axs, :xlabelvisible, true)
setproperty!.(axs, :ylabelvisible, true)
axs[2].yaxisposition = :right

axs[1].title = "Experimental 𝐶. 𝑒𝑙𝑒𝑔𝑎𝑛𝑠"

xu_v_lines = lines!(axs[1], xu_v_null; linewidth = 1, color = wong[1])
xu_w_lines = lines!(axs[1], xu_w_null; linewidth = 1, color = wong[2])

axs[2].title = "FitzHugh-Nagumo"

lines!(axs[2], fh_v_null; linewidth = 1, color = wong[1])
lines!(axs[2], fh_w_null; linewidth = 1, color = wong[2])


# Now, we create the legend.

leg = layout[2, 1:2] = LLegend(scene,
                [xu_v_lines, xu_w_lines],
                ["𝑑𝑣 = 0", "𝑑𝑤 = 0"];
                framevisible = false,
                linewidth = 1,
                linepoints = Point2f0[(.5, .5), (1, .5)],
                rowgap = -6,
)

# The legend doesn't need much space...80 points should be enough.
rowsize!(layout, 2, 30)
rowgap!(layout, 0)

save("exp_null_tmp.pdf", scene)

# save(plotsdir("Dynamics/Nullclines/Neuron model nullclines/nm_nullclines_makie.pdf"), scene)

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
#              Nullclines              #
########################################

keener_vnull = Point2f0[ # manually expand linear interp
    (-2, -2 - (-2)^3/3),
    (-1, -1 - (-1)^3/3),
    (1, 1 - 1^3/3),
    (2, 2 - (2)^3/3),
]

########################################
#                Figure                #
########################################

scene, layout = layoutscene(7;
    # 424 is the text width in the tex document,
    # which can be found using `\the\textwidth`.
    resolution = (424, 220),
    rowgap = 0,
    colgap = 0
)

axs = layout[1:2, 1:3] = [LAxis(scene) for row in 1:2, col in 1:3]
rowgap!(layout, 5)

setproperty!.(axs[1, 2:3], :xlabel, "𝑡 (arb.)")
setproperty!.(axs[1, 2:3], :ylabel, "𝑣")
setproperty!.(axs[1, 2:3], :xticklabelsvisible, false)
setproperty!.(axs[1, 2:3], :yticklabelsvisible, false)
setproperty!.(axs[2, 2:3], :xlabel, "𝑡")
setproperty!.(axs[2, 2:3], :ylabel, "𝑉")

save("tmp.pdf", scene)

##############################
#         Nullclines         #
##############################

fhn_ax = axs[1, 1]
fhn_ax.xlabel = "𝑣"
fhn_ax.ylabel = "𝑤"
fhn_ax.title = "FHN nullclines\n "

lines!(fhn_ax, fh_v_null; linewidth = 1, color = wong[1])
lines!(fhn_ax, fh_w_null; linewidth = 1, color = wong[2])
save("tmp.pdf", scene)


keener_ax = axs[2, 1]
keener_ax.xlabel = "𝑣"
keener_ax.ylabel = "𝑤"
keener_ax.title = "Keener nullclines"

lines!(keener_ax, keener_vnull; linewidth = 1, color = wong[1])
lines!(keener_ax, fh_w_null; linewidth = 1, color = wong[2])

setproperty!.((fhn_ax, keener_ax), :xticklabelsvisible, false)
setproperty!.((fhn_ax, keener_ax), :yticklabelsvisible, false)
save("tmp.pdf", scene)

##############################
#      Head oscillator       #
##############################

sim_double_ax  = axs[1, 2]
anal_double_ax = axs[2, 2]

sim_double_ax.title = "Head oscillator\n(negative coupling)"


lines!(sim_double_ax, sim_double_v0; linewidth = 1, color = wong[1])
lines!(sim_double_ax, sim_double_v1; linewidth = 1, color = wong[2])

lines!(anal_double_ax, anal_double_v0; linewidth = 1, color = wong[1])
lines!(anal_double_ax, anal_double_v1; linewidth = 1, color = wong[2])

anal_double_ax.xticks[] = MakieLayout.LinearTicks(3)

save("tmp.pdf", scene)

##############################
#     Descending pathway     #
##############################

sim_single_ax  = axs[1, 3]
anal_single_ax = axs[2, 3]

sim_single_ax.title = "Descending pathway\n(positive coupling)"


lines!(sim_single_ax, sim_single_v0; linewidth = 1, color = wong[1])
lines!(sim_single_ax, sim_single_v1; linewidth = 1, color = wong[2])

lines!(anal_single_ax, anal_single_v0; linewidth = 1, color = wong[1])
lines!(anal_single_ax, anal_single_v1; linewidth = 1, color = wong[2])

anal_single_ax.xticks[] = MakieLayout.LinearTicks(3)

save("tmp.pdf", scene)

colgap!(layout, 2, 10)
colgap!(layout, 1, 10)
rowgap!(layout, 1, 0)
save("tmp.pdf", scene)

resize!(scene, 424, 250)
anal_single_ax.title = "Real descending pathway"
anal_double_ax.title = "Real head oscillator"
save("tmp.pdf", scene)

labels = [
    LText(
        scene,
        string('(', 'a' + j + (i-1)*3, ')');
        textsize = 10,
        font = "Fira Sans Bold",
        halign = :left,
        padding = (0,0f0,0,0)
    )
    for i in 1:2, j in 0:2
]

for ind in CartesianIndices(labels)
    # ind cannot be bc'ed so convert to Tuple
    layout[Tuple(ind)..., TopLeft()] = labels[ind]
end


save("tmp.pdf", scene)

setproperty!.(labels, :padding, Ref((0f0, 0f0, 0f0, 0f0)))
setproperty!.(labels, :halign, :left)

save("tmp.pdf", scene)


save("papers/Paper/figures/anal_sim_comp/anal_sim_comp.pdf", scene)
