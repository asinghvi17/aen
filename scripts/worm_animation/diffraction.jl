trange = LinRange(tspan..., 30000)

timeseries = sol(trange)

signal = eff_signal(timeseries)

animation_range = eachindex(trange)#1000:.1:3000

timeobs = Node(1)

# linobs = @lift signal_to_worm(signal[round(Int, $timeobs * length(trange)/(extrema(trange)[2] - extrema(trange)[1]), RoundUp), :])
linobs = @lift signal_to_worm(signal[$timeobs, :])

sc = lines(
    linobs;
    linewidth = 40,
    color = :gray,
    scale_plot = false,
    limits = Rect2D{Float32}(-4, -10, 14, 20),
    show_axis = false,
    resolution = (800, 800),
);
save("tmp1.png", sc)

@time record(sc, "worm_12.mp4"; framerate = 60) do io
    @time Juno.@progress name="Recording" for i in animation_range
        timeobs[] = i
        recordframe!(io)
    end
end

cairo_screen = CairoMakie.CairoScreen(sc)
img = lift(timeobs) do _
    AbstractPlotting.colorbuffer(cairo_screen)
end

gimg = @lift Gray.($img)

transformed_image = @lift fft(Float64.($gimg))

fftvals = @lift fftshift(abs2.($transformed_image))

ascene, alayout = layoutscene()

worm_ax = alayout[1, 1] = LAxis(ascene; aspect = DataAspect())

hidedecorations!(worm_ax)
hidespines!(worm_ax)
lines!(worm_ax, linobs; linewidth = 40, color = (:black, 0.5))

worm_ax.targetlimits[] = Rect2D{Float32}(-2, -6.5, 14, 13)

fft_ax =
    alayout[1, 2] =
        LAxis(ascene; backgroundcolor = :black, aspect = DataAspect())

heatmap!(fft_ax, @lift log.($fftvals); colormap = [:black, :green])

save("worm_diff.pdf", ascene)
save("worm_diff.png", ascene)

record(ascene, "worm_diff.mp4") do io
    Juno.@progress for t in animation_range
        timeobs[] = t
        recordframe!(io)
    end
end
