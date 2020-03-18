using Plots,
      CSV,
      DataFrames,
      PlotThemes

theme(:dark)

gr()

ch1 = CSV.File("CH1.csv") |> DataFrame |> y -> filter(x -> x[1] < 4.0, y)
ch2 = CSV.File("CH2.csv") |> DataFrame |> y -> filter(x -> x[1] < 4.0, y)

ch1[2]

p = plot(ch1[1], ch1[2], xlabel = "t (s)", ylabel = "V", label = "Capacitor")
plot!(ch2[1], ch2[2], label = "Raw")

savefig(p, "output-dark.pdf")
μ
