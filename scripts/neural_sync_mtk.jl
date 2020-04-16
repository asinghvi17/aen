using ModelingToolkit, DifferentialEquations

using Makie, MakieLayout
@parameters t g e b
@variables v(t) w(t) F(t)
@derivatives D'~t
single_neuron_eqs = [
    D(v) ~ min(max(-2 - v, v), 2 - v) - w - F, # add the flux term
    D(w) ~ e * (v - g * w + b)
]
n1 = ODESystem(single_neuron_eqs, t, [v,w,F],[g,e,b], name = :n1)
n2 = ODESystem(single_neuron_eqs, t, [v,w,F],[g,e,b], name = :n2)
@parameters D Dk
connections = [
    0 ~ n1.F - D * Dk * max(n1.v - n2.v, 0),
    0 ~ n2.F - D * max(n2.v - n1.v, 0)
]
connected = ODESystem(connections,t,[],[D,Dk],systems=[n1,n2],name=:connected)
equations(connected)
u0 = [
    n1.v => -2,
    n1.w => -2/3,
    n1.F => 0,
    n2.v => -2,
    n2.w => -.7,
    n2.F => 0,
]
tspan = (0.0, 1750.0)
p0 = [
    D => 0.047,
    Dk => 0.5,
    n1.g => 0.8,
    n1.e => 0.04,
    n1.b => 0,
    n2.g => 0.8,
    n2.e => 0.04,
    n2.b => 0.2,
]

func = ODEFunction(connected, first.(u0), first.(p0))
prob   = ODEProblem(func, last.(u0), tspan, last.(p0))

scene, layout = layoutscene()

ax = layout[1, 1] = LAxis(scene)

Dkr = LinRange(0, 1, 100)

sl = layout[2, 1] = LSlider(scene; range = Dkr, startvalue = .5)

p = lift(sl.value) do Dkv
    return [
        D => 0.047,
        Dk => Dkv,
        n1.g => 0.8,
        n1.e => 0.04,
        n1.b => 0,
        n2.g => 0.8,
        n2.e => 0.04,
        n2.b => 0.2,
    ]
end

dnet = lift(sl.value) do dk
    "Dₙₑₜ = " * string(round(last(p[][1]) * dk; sigdigits=3))
end

txt = layout[0, 1] = LText(scene, dnet; textsize = 30, width = Auto(false))

ts = LinRange(1000, 1750, 10000)

sol = @lift solve(remake(prob; p = last.($p)) ,Rodas4(),abstol=1e-12,reltol=1e-12)(ts)

plot_v0 = @lift $sol[1, :]
plot_v1 = @lift $sol[4, :]

lines!(ax, ts, plot_v0; linewidth = 1, color = :blue)
lines!(ax, ts, plot_v1; linewidth = 1, color = :red)

scene

record(scene, "mtk_coupling.mp4", LinRange(0, 1, 240)) do i
    sl.value[] = i
end

# TODO raise an issue about the basic init alg on DAEs
