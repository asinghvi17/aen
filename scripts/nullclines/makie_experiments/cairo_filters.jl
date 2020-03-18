
using CairoMakie, Cairo

data = rand(10, 10)

for filter in [Cairo.FILTER_BEST, Cairo.FILTER_BILINEAR, Cairo.FILTER_FAST, Cairo.FILTER_GAUSSIAN, Cairo.FILTER_GOOD, Cairo.FILTER_NEAREST]

    @eval CairoMakie begin function draw_image(scene, screen, attributes)
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
            interp = interp ? $filter : Cairo.FILTER_NEAREST
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

    title(heatmap(data; show_axis = false, scale_plot=false, interpolate = true), string(filter)) |> AbstractPlotting.FileIO.save("$(string(filter)).svg")
    title(heatmap(data; show_axis = false, scale_plot=false, interpolate = true), string(filter)) |> AbstractPlotting.FileIO.save("$(string(filter)).png")
    title(heatmap(data; show_axis = false, scale_plot=false, interpolate = true), string(filter)) |> AbstractPlotting.FileIO.save("$(string(filter)).pdf")

    run(`pdftocairo -png $(string(filter)).pdf $(string(filter))_pdf.png`)
    run(`pdftocairo -svg $(string(filter)).pdf $(string(filter))_svg.png`)
end
