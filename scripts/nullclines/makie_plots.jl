AbstractPlotting.set_theme!(
    font = "CMU Serif Roman", # same font as the document.  Can be changed arbitrarily.
    fontsize = 12,            # 12-point font.
    titlegap = 5,
    xgridwidth = 0.5,
    ygridwidth = 0.5,
)

const a4_resolution = (595, 842)

################################################################################
#                           Fitzhugh-Nagumo dynamics                           #
################################################################################

fh_exc = (x...) -> fitzhugh_nagumo(x...; i_ext = 0.9)
fh_osc = (x...) -> fitzhugh_nagumo(x...)

iso_vosc = isoclines(-2..2, -2..2, fh_osc; ind=1)
iso_wosc = isoclines(-2..2, -2..2, fh_osc; ind=2)

iso_vexc = isoclines(-2..2, -2..2, fh_exc; ind=1)
iso_wexc = isoclines(-2..2, -2..2, fh_exc; ind=2)

scene, layout = layoutscene(
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

leg.grid.content[1].content[1, 2] = LTeX(leg.child,
    raw"dv = \frac{v^3}{3} - v + w + I_{ext} = 0"
)

leg.grid.content[1].content[2, 2] = LTeX(leg.child,
    raw"dw = e \cdot \left(v- g \cdot w + b\right) = 0"
)

save(plotsdir("Dynamics/Nullclines/FHN_nullclines/fhn_dynamics_makie.pdf"), scene)

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
)

# The legend doesn't need much space...80 points should be enough.
rowsize!(layout, 2, 45)

save(plotsdir("Dynamics/Nullclines/Neuron model nullclines/nm_nullclines_makie.pdf"), scene)

################################################################################
#                     Comparison of analog and simulation                      #
################################################################################
