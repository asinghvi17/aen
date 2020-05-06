
using DifferentialEquations,
      ParameterizedFunctions

using Makie, Observables

# create ode

function f(v)
    # min(max(-2 - v, v), 2 - v)
    v - v^3/3
end
FHND = @ode_def begin

    dv₀ = f(v₀) - w₀ + D * Dk * v₁

    dw₀ = e₀ * (v₀ - g * w₀ + b₀)

    dv₁ = f(v₁) - w₁ + D * v₀

    dw₁ = e₁ * (v₁ - g * w₁ + b₁)

end D Dk g e₀ b₀ e₁ b₁

# sl = layout[2, 1] = LSlider(oscene; range = Dkr, startvalue = .5)
Dk = Observable(1.0)

oscene



x0 = -1.0
y0 = -0.51
x1 =  1.0
y1 = -0.49

u0 = [x0, y0, x1, y1]

D  = -0.1
g  = 0.8
e0 = 0.04
b0 = 0.46
e1 = 0.04
b1 = 0.47

p = @lift([D, $Dk, g, e0, b0, e1, b1])

ts = LinRange(1000, 1750, 1000)

eprob = ODEProblem(FHND, u0, (000.0, 3000.0), to_value(p), AutoTsit5(Vern9()))

osol = @lift solve(remake(eprob; p = $p), save_idxs = [1, 3])

sol = @lift $osol(ts)

plot_v0 = @lift $sol[1, :]
plot_v1 = @lift $sol[2, :]




using MakieLayout

oscene, layout = layoutscene()

ax = layout[1, 1] = LAxis(oscene)

txt = layout[0, 1] = LText(oscene, @lift("Dₙₑₜ = " * string(round(D * $Dk; sigdigits=3))); tellwidth = false, textsize = 30)

Dkr = LinRange(0, 1, 240)

lines!(ax, ts, plot_v0; color=:blue)
lines!(ax, ts, plot_v1; color=:red)

record(oscene, "coupling.mp4", LinRange(0, 1, 240)) do i
    Dk[] = i
end
