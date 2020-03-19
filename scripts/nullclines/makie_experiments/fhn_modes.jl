include("nullclines.jl")

scene, layout = layoutscene(
    resolution = (1000, 800),
    font = "Calibri",
    rowsizes = [Auto(), Fixed(100)],
);

axs = [
    LAxis(
        scene;
        aspect = AxisAspect(1),
        # titlesize = 40 * 6,
        # titlegap = 20 * 6,
        xticks = WilkinsonTicks(6; k_min = 5),
        yticks = WilkinsonTicks(5;  k_min = 5),
        xlabel = "v",
        ylabel = "w",
        # xlabelsize = 40 * 5,
        # ylabelsize = 40 * 5,
        xticksvisible = false,
        yticksvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,
        xticksize = 0,
        yticksize = 0,
        width = Auto(false),
        height = Auto(false)
    )
    for _ in 1:2
]

layout[1, 1] = axs[1]
layout[1, 2] = axs[2]
# layout[1, 1:2] .= axs[1:2]

axs[1].title = "Excitable regime"
axs[2].title = "Oscillatory regime"
# axs[3].title = "Unstable regime"

#axs[2].ylabelvisible = false
axs[2].yaxisposition = :right

fo = (x...) -> fitzhugh_nagumo(x...; i_ext = 0.5)

iso_vo = isoclines(-2..2, -2..2, fo; ind=1)
iso_wo = isoclines(-2..2, -2..2, fo; ind=2)
lines!(axs[1], iso_vo; color = wong[5], linewidth = 3)
lines!(axs[1], iso_wo; color = wong[6], linewidth = 3)

iso_v = isoclines(-2..2, -2..2, fitzhugh_nagumo; ind=1)
iso_w = isoclines(-2..2, -2..2, fitzhugh_nagumo; ind=2)

# steady-state
iso_vo = isoclines(-2..2, -2..2, fitzhugh_nagumo; ind=1)
iso_wo = isoclines(-2..2, -2..2, fitzhugh_nagumo; ind=2)
lines!(axs[2], iso_v; color = wong[5], linewidth = 3)
lines!(axs[2], iso_w; color = wong[6], linewidth = 3)


# phase dynamics
arrows!(axs[1], LinRange(-2..2, 10), LinRange(-2..2, 10), (x...) -> fo(x...) ./ 7; arrowsize = 0.07, arrowhead = Point2f0[(0, 0), (1, 0), (0.5, 1)])
arrows!(axs[2], LinRange(-2..2, 10), LinRange(-2..2, 10), (x...) -> fitzhugh_nagumo(x...) ./ 7; arrowsize = 0.07, arrowhead = GeometryTypes.Polygon(Point2f0[(0, 0), (1, 0), (0.5, 1)]))

# draw the legend
layout[2, 1:2] = leg = LLegend(
    scene,
    plots(axs[1]),
    ["v-nullcline", "w-nullcline", "gradient"];
    titlesize = 0,
    titlevisible = false,
    valign = :top,
    strokecolor = RGBAf0(0, 0, 0, 0),
    strokewidth = 0,
    framevisible = false,
    height = Auto(true)
)

save("fhn_modes.pdf", scene)
save("fhn_modes.png", scene)

save("poly.png", Makie.poly(Point2f0[(0, 0), (1, 0), (0.5, 1)]))
CairoMakie.activate!()
# # Ideal interface
# scene = Scene()
#
# iso_v = isoclines(-2..2, -2..2, fitzhugh_nagumo; ind=1)
# iso_w = isoclines(-2..2, -2..2, fitzhugh_nagumo; ind=2)
#
# # steady-state
# ss = lines(iso_v; color = wong[5], linewidth = 3, title = "Steady state")
# lines!(iso_w .+ Ref(Point2f0(-0.5, 0)); color = wong[6], linewidth = 3)
#
# # oscillatory
# os = lines(iso_v; color = wong[5], linewidth = 3, title = "Steady state")
# lines!(iso_w; color = wong[6], linewidth = 3)
#
# hbox(vbox(os, ss), legend(os, ["v", "w"]))
arrows(LinRange(-2..2, 10), LinRange(-2..2, 10), (x...) -> fo(x...) ./ 7; arrowsize = 0.07, arrowhead = Point2f0[(0, 0), (1, 0), (0.5, 1)])

@eval CairoMakie begin
    function draw_image(scene, screen, attributes)
        ctx = screen.context
        image = attributes[3][]
        x, y = attributes[1][], attributes[2][]
        model = attributes[:model][]
        imsize = (extrema_nan(x), extrema_nan(y))
        xy_ = project_position(scene, Point2f0(first.(imsize)), model)
        xymax_ = project_position(scene, Point2f0(last.(imsize)), model)
        xy = min.(xy_, xymax_)
        xymax = max.(xy_, xymax_)
        w, h = xymax .- xy
        interp = to_value(get(attributes, :interpolate, true))
        interp = interp ? Cairo.FILTER_BEST : Cairo.FILTER_NEAREST
        s = to_cairo_image(image, attributes)
        Cairo.rectangle(ctx, xy..., w, h)
        Cairo.save(ctx)
        Cairo.translate(ctx, xy[1], xy[2])
        Cairo.scale(ctx, w / s.width, h / s.height)
        Cairo.set_source_surface(ctx, s, 0, 0)
        p = Cairo.get_source(ctx)
        # Set filter doesn't work!?
        Cairo.pattern_set_filter(p, interp)
        Cairo.fill(ctx)
        Cairo.restore(ctx)
    end
end
