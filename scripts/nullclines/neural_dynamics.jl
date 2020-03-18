include("nullclines.jl")

scene, layout = layoutscene(
    resolution = (1440, 1000),
    font = "Calibri",
    rowsizes = [Auto(), Fixed(100)],
);

axs = [
    LAxis(
        scene;
        aspect = AxisAspect(1),
        # titlesize = 30,
        # titlegap = 20 * 6,
        xticks = WilkinsonTicks(6; k_min = 5),
        yticks = WilkinsonTicks(5;  k_min = 5),
        xticksvisible = false,
        yticksvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,
        xticksize = 0,
        yticksize = 0,
    )
    for _ in 1:3
]

layout[1, 1] = axs[1]
layout[1, 2] = axs[2]
layout[1, 3] = axs[3]
axs[1].title = "Real C. elegans"
axs[2].title = "Morris-Lecar"
axs[3].title = "FitzHugh-Nagumo"

dlw = 3

xs = Base.LinRange(-60..0, 500)
lines!(axs[1], xs, xu_vnull.(xs); color = wong[5], linewidth = dlw)
lines!(axs[1], xs, xu_wnull.(xs); color = wong[6], linewidth = dlw)

ylims!(axs[1], (0.5, 0.9))

lines!(axs[2], isoclines(-100..100, -1..2, morris_lecar; ind = 1), color = wong[5], linewidth = dlw)
lines!(axs[2], isoclines(-100..100, -1..2, morris_lecar; ind = 2), color = wong[6], linewidth = dlw)

lines!(axs[3], isoclines(-2..2, -1..1, fitzhugh_nagumo; ind = 1), color = wong[5], linewidth = dlw)
lines!(axs[3], isoclines(-2..2, -1..1, fitzhugh_nagumo; ind = 2), color = wong[6], linewidth = dlw)

xlims!(axs[2], (-75, 75))
ylims!(axs[2], (-0.05, 0.6))

layout[2, 1:3] = leg = LLegend(
    scene,
    plots(axs[2]),
    ["v-nullcline", "w-nullcline"];
    # linewidth = 12 * 3,
    titlesize = 0,
    titlevisible = false,
    # labelsize = 30 * 6,
    # patchsize = (40f0, 40f0) .* 4,
    valign = :top,
    # height = Fixed(100 * 8),
    strokecolor = RGBAf0(1, 0, 0, 0),
    strokewidth = 0,
    framevisible = false,
    height = Auto(false)
    # rowgap = 100
)

save("neural_dynamics.pdf", scene)


# resize!(scene, (1440, 1000))
