set_theme!(font = "Fira Sans")

trange = LinRange(1000, 1600, 10000)

timeseries = sol(trange)

signal = eff_signal(timeseries)

animation_range = 1000:0.8:1600

time_to_worm(timeobs) = Point2f0.(signal_to_worm(signal[round(Int, timeobs * length(trange)/(extrema(trange)[2] - extrema(trange)[1]), RoundUp), :]))

scene, layout = layoutscene(10; resolution = (1024, 800))

dv_subgl = layout[1, 1] = GridLayout(scene; alignmode = Outside())

ventral_ax = dv_subgl[1, 1] = LAxis(scene; title = "Ventral neurons")
# dorsal_ax  = dv_subgl[2, 1] = LAxis(scene; title = "Dorsal neurons", ylabel = "𝑣", ylabelfont = "CMU Serif Math")
# linkaxes!(ventral_ax, dorsal_ax)
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

# for neuron in dorsal_neuron_list
#     lines!(
#         dorsal_ax,
#         trange,
#         timeseries[neuron, :];
#         color = neuron_colors[round(Int, neuron/4)],
#         linewidth = 1
#     )
# end

leg = dv_subgl[1, 2] = LLegend(scene, ventral_ax.scene.plots, [string(i) for i in eachindex(neuron_colors)], "Neuron pair\nnumber"; nbanks = 2)

save("temp.png", scene)

worm_subgl = layout[2, 1] = GridLayout(scene)

worm_time_axs = [LAxis(scene; aspect = DataAspect()) for i in 1:6]

worm_subgl[1, 1:6] = worm_time_axs

hidedecorations!.(worm_time_axs)

save("temp2.png", scene)


for (ind, time) in zip(1:6, LinRange(1, 50, 6))
    ax = worm_time_axs[ind]
    for plot in ax.scene.plots
        delete!(ax.scene, plot)
    end
    lines!(ax, time_to_worm(time); linewidth = 10)
    ax.title = "t = $(round(Int, time))"
end

linkaxes!(worm_time_axs...)

picture_subgl = layout[3, 1] = GridLayout(scene)
picture_axs = picture_subgl[1, 1:3] = [LAxis(scene; aspect = DataAspect()) for i in 1:3]
hidedecorations!.(picture_axs)
picture_dir = "figures/Worms/worm_pictures"
pictures = load.(joinpath.(picture_dir, readdir(picture_dir)))
for (ax, pic) in zip(picture_axs, pictures)
    image!(ax, 0..size(pic, 1), 0..size(pic, 2), pic)
end
tightlimits!.(picture_axs)


save("temp3.png", scene)

label_a = layout[1, 1, TopLeft()] = LText(
    scene, "(a)";
    textsize = 20,
    font = "Fira Sans Bold",
    halign = :right
)


label_b = layout[2, 1, TopLeft()] = LText(
    scene, "(b)";
    textsize = 20,
    font = "Fira Sans Bold",
    halign = :right
)

label_c = layout[3, 1, TopLeft()] = LText(
    scene, "(c)";
    textsize = 20,
    font = "Fira Sans Bold",
    halign = :right
)

label_c.padding = (0, 0, 5, 0)

save("temp4.png", scene)

save("figures/Worms/worm_neuron_dash.pdf", scene)
save("papers/IPoLS poster/figures/worm_neuron_dash/worm_neuron_dash.pdf", scene)
