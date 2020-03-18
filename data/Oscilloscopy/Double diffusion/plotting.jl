using DelimitedFiles,
      Plots,
      Loess,
      LaTeXStrings

gr(size = (800, 600))

dir = "../One way positive diffusive coupling/ALL0006"

if !isfile(joinpath(dir, "CH1.csv"))
    run(pipeline(pipeline(`cat $(joinpath(dir, "F000$(dir[end])CH1.csv"))`, `cut -d',' -f4-5`), stdout = joinpath(dir, "CH1.csv")))
end


if !isfile(joinpath(dir, "CH2.csv"))
    run(pipeline(pipeline(`cat $(joinpath(dir, "F000$(dir[end])CH2.csv"))`, `cut -d',' -f4-5`), stdout = joinpath(dir, "CH2.csv")))
end


function readChData(dir::String)

    ch1 = readdlm(joinpath(dir, "CH1.csv"), ',', Float64, '\n')

    ch2 = readdlm(joinpath(dir, "CH2.csv"), ',', Float64, '\n')

    rawset = open(joinpath(dir, "F000$(dir[end])TEK.SET"), "r") do f

                read(f, String) |> x -> split(x, ";")

            end

    return ch1, ch2, rawset

    end

ch1, ch2, set = readChData(dir)

tr = r":HORIZONTAL:MAIN:SCALE *"

tscale = filter(x -> occursin(tr, x), set)... |> x -> x[23:end] |> String |> Meta.parse |> eval

ts = (ch1[:, 1] .- minimum(ch1[:, 1])) .* tscale

## saving

writedlm(joinpath(dir, "ch1.dat"), hcat(ts[1:], ch1[:, 2]))
writedlm(joinpath(dir, "ch2.dat"), hcat(ts, ch2[:, 2]))

## plotting

p = plot(
        ts,
        loess(ch1[:, 1], ch1[:, 2], span=0.08) |> x -> predict(x, ch1[:, 1]),
        # ch1[:, 2],
        linewidth=4,
        xlabel = "Time (s)",
        ylabel = "Voltage (V)",
        label = L"V_0",
        title = "Real circuit response (smoothed)"
        ) # plot channel 1 as the ventral

plot!(
        p,
        ts,
        loess(ch2[:, 1], ch2[:, 2], span=0.08) |> x -> predict(x, ch2[:, 1]),
        # ch2[:, 2],
        linewidth=4,
        label = L"D_0"
        ) # plot channel 2 as the dorsal

savefig(p, joinpath(dir, "smooth_response.pdf"))

smoothData = hcat(ts,loess(ch1[:, 1], ch1[:, 2], span=0.08) |> x -> predict(x, ch1[:, 1]), loess(ch2[:, 1], ch2[:, 2], span=0.08) |> x -> predict(x, ch2[:, 1]))

writedlm(joinpath(dir, "smooth.dat"), vcat(["t" "v0" "v1"], smoothData))
