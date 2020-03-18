using LaTeXStrings
using DynamicalSystems

sys = Systems.roessler()

T = 100.0

sol = trajectory(sys, T)
using Plots
plot(sol)

tr = LinRange(0, T, size(sol, 1))

using BenchmarkTools

@benchmark periodogram(sol[:, 1];  fs = size(sol, 1)/T)
@benchmark lombscargle(tr, sol[:, 1])

function plot_periodogram(tr, tseries, periods, powers; pxlim = ())

    timeseries = plot(
                    tr,
                    tseries,
                    xlabel = L"t",
                    ylabel = L"V",
                    title = "Timeseries",
                    legend = false
                )

    periodogram = plot(
        periods,
        powers,
        # xlims = (0, 50),
        legend=false,
        xlabel = L"T",
        ylabel = L"P",
        title = "Periodogram",
        xlim = pxlim
        )

    plot(timeseries, periodogram, layout = @layout [a;b])

end

using Plots

plot(tr, sol[:, 1])

using LombScargle

lsp = LombScargle.plan(tr, sol[:, 1])
ls = lombscargle(lsp)
findmaxperiod(ls, [T/size(sol, 1)*2, Inf])

plot_periodogram(tr, sol[:, 1], period(ls), LombScargle.power(ls), pxlim = (0, 100))


using DSP



dsp = periodogram(sol[:, 1];  fs = size(sol, 1)/T)

fr = DSP.Periodograms.freq(dsp)

pds = 1 ./ fr

pwrs = DSP.Periodograms.power(dsp)

pds[findmax(pwrs)[2]]

plot_periodogram(tr, sol[:, 1], pds, pwrs)

# """
# """
# function find_period(ds::DynamicalSystem, T::Real, )

using ParameterizedFunctions

FHN = @ode_def begin
    dv = min(max(-2 - v, v), 2 - v) - w
    dw =  e * (v - g * w + b)
end e b g


v0 = -2
w0 =-0.6667
u0 = [v0, w0]

g  = 0.8
e = 0.04
b = 0
p0 = [e, b, g]

T = 750.0

sys = ContinuousDynamicalSystem(FHN, u0, p0)
sol = trajectory(sys, T)
ts = LinRange(0, 750, size(sol, 1))

dsp = periodogram(sol[:, 1];  fs = size(sol, 1)/T)

fr = DSP.Periodograms.freq(dsp)

pds = 1 ./ fr

pwrs = DSP.Periodograms.power(dsp)

pds[findmax(pwrs)[2]]

plot_periodogram(ts, sol[:, 1], pds, pwrs)

plot(ts, [sol[:, 1], sol[:, 2]])
sol[:, 1:2] |> Matrix
