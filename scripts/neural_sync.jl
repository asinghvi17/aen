
using DifferentialEquations,
      ParameterizedFunctions

using Makie, Observables

# create ode
FHND = @ode_def begin

    dv₀ = min(max(-2 - v₀, v₀), 2 - v₀) - w₀ + Dk*D*max(v₁ - v₀, 0)

    dw₀ = e₀ * (v₀ - g * w₀ + b₀)

    dv₁ = min(max(-2 - v₁, v₁), 2 - v₁) - w₁ + D*max(v₀ - v₁, 0)

    dw₁ = e₁ * (v₁ - g * w₁ + b₁)

end D Dk g e₀ b₀ e₁ b₁


using MakieLayout

oscene, layout = layoutscene()

ax = layout[1, 1] = LAxis(oscene)

Dkr = LinRange(0, 1, 100)

# sl = layout[2, 1] = LSlider(oscene; range = Dkr, startvalue = .5)
Dk = Observable(1.0)

oscene



x0 = -2
y0 =-0.6667
x1 = -2.0
y1 = -0.7

u0 = [x0, y0, x1, y1]

D  = -0.047
g  = 0.8
e0 = 0.04
b0 = 0
e1 = 0.04
b1 = 0.2 # was 0.25

txt = layout[0, 1] = LText(oscene, @lift("Dₙₑₜ = " * string(round(D * $Dk; sigdigits=3))); width = Auto(false), textsize = 30)

p = @lift [D, $Dk, g, e0, b0, e1, b1]

ts = LinRange(1000, 1750, 10000)

prob = ODEProblem(FHND, u0, (000.0, 2000.0), to_value(p), AutoTsit5(Vern9()))

sol = @lift solve(remake(prob; p = $p))(ts)

plot_v0 = @lift $sol[1, :]
plot_v1 = @lift $sol[3, :]

lines!(ax, ts, plot_v0; color=:red)
lines!(ax, ts, plot_v1; color=:blue)

record(oscene, "coupling.mp4", LinRange(0, 1, 240)) do i
    Dk[] = i
end
