
trange = LinRange(tspan..., 10000)

timeseries = sol(trange)

signal = eff_signal(timeseries)

animation_range = 500:0.8:1600

timeobs = Node(0.1)

linobs = @lift signal_to_worm(signal[round(Int, $timeobs * length(trange)/(extrema(trange)[2] - extrema(trange)[1]), RoundUp), :])


scene, layout = layoutscene(10; resolution = (round(Int, 750 * 2.5), 900))
rowgap!(layout, 10)
colgap!(layout, 2)

worm_ax = layout[1:2, 1] = LAxis(
    scene;
    # targetlimits = Rect2D{Float32}(-4, -6, 10, 12),
    aspect = DataAspect(),

)
hidedecorations!(worm_ax)
hidespines!(worm_ax)
# worm_ax.xgridvisible = true
# worm_ax.ygridvisible = true
worm_ax.panbutton    = Mouse.middle
worm_ax.xzoomlock    = true
worm_ax.yzoomlock    = true

ventral_ax = layout[1, 2] = LAxis(scene)
hidedecorations!(ventral_ax)
ventral_ax.xgridvisible = true
ventral_ax.ygridvisible = true

dorsal_ax = layout[2, 2] = LAxis(scene)
hidedecorations!(dorsal_ax)
dorsal_ax.xgridvisible = true
dorsal_ax.ygridvisible = true
dorsal_ax.xticksvisible = true
dorsal_ax.xticklabelsvisible = true
dorsal_ax.xtickalign = 1.0

linkxaxes!(ventral_ax, dorsal_ax)
linkyaxes!(ventral_ax, dorsal_ax)

neuron_colors = cgrad(:isolum)[LinRange(0, 1, NUM_NEURONS÷2)];

ventral_neuron_list = 1:4:(NUM_NEURONS*2)
dorsal_neuron_list  = 3:4:(NUM_NEURONS*2)

# timeseries and signal defined in utils

for neuron in ventral_neuron_list
    lines!(
        ventral_ax,
        trange,
        timeseries[neuron, :];
        color = neuron_colors[round(Int, neuron/4, RoundUp)],
        linewidth = 1
    )
end

for neuron in dorsal_neuron_list
    lines!(
        dorsal_ax,
        trange,
        timeseries[neuron, :];
        color = neuron_colors[round(Int, neuron/4)],
        linewidth = 1
    )
end

# display(scene)

ventral_rect = ventral_ax.targetlimits
dorsal_rect  = dorsal_ax.targetlimits

ventral_timeline = @lift Point2f0[
    ($timeobs, bottom($ventral_rect)),
    ($timeobs,    top($ventral_rect))
]

dorsal_timeline = @lift Point2f0[
    ($timeobs, bottom($dorsal_rect)),
    ($timeobs,    top($dorsal_rect))
]

lines!(ventral_ax, ventral_timeline; linewidth = .7)
lines!(dorsal_ax,  dorsal_timeline;  linewidth = .7)



xlims!(ventral_ax, extrema(animation_range))

timeobs[] = 1400

lines!(worm_ax, linobs; linewidth = 30)
display(scene)
worm_ax.targetlimits[] = Rect2D{Float32}(-7, -10, 14, 20)

timeobs[] = 1410

##
record(scene, "worm_dash_12pairs_cb_e.mp4"; framerate = 30) do io
    Juno.@progress for t in animation_range
        timeobs[] = t
        recordframe!(io)
    end
end
