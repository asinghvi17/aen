# Nullclines Draft 1

```julia id=13373429-787d-4bb8-9172-4b3d576ffe61
;rm /root/.julia/environments/v1.3/Project.toml /root/.julia/environments/v1.3/Manifest.toml 
```

```julia id=5c61c477-d20c-4c44-95b6-3dadadace86a
]add GLMakie#master AbstractPlotting#master StatsMakie#master MakieLayout#master GeoMakie#master Observables#master Polynomials; build GLMakie
```

```julia id=eec0d8a5-c89a-47e3-9977-44fad72b5505
using AbstractPlotting, GLMakie
AbstractPlotting.inline!(true) # to display PNGs in the journal
AbstractPlotting.set_theme!(Theme(font = "DejaVu Sans")) # better font
scatter(rand(4)) # test, and compile to get ahead of TTFP
```

![result][nextjournal#output#eec0d8a5-c89a-47e3-9977-44fad72b5505#result]

```julia id=8012a0b5-d4dd-4276-94cd-07cc296c67c0
using Polynomials
using MakieLayout, GeoMakie
import AbstractPlotting: xlims!, ylims!, setlims!

# get a colour palette
wong = AbstractPlotting.default_palettes.color[]

# define differential equations
function morris_lecar(v, w)
    # params for Morris-Lecar
    gca = 4.4
    ϕ   = 0.04
    return [
        # voltage
        (1/20) * (90 - gca * 0.5 * (1 + tanh((v + 1.2)/18)) * (v - 120) - 8 * w * (v + 84) - 2 * (v + 60)),
        # potassium current (gating)
        ϕ * (0.5 * (1 + tanh((v - 2)/30)) - w) * cosh((v - 2)/60)
    ]
end

lin_activ(v) = min(max(-2 - v, v), 2 - v)

function fitzhugh_nagumo(v, w)
    e, g, b = (0.04, 0.8, 0) # unpack parameters
    return [
        w - (-v^3/3 + v),    # dv₀
        e * (v - g * w + b) # dw₀
    ]
end

function nullcline!(scene, f, domain; samples = 500)
    # extract domains
    xmin, xmax = domain[Left()], domain[Right()]
    ymin, ymax = domain[Bottom()], domain[Top()]

    # create axis vectors
    xs = Base.LinRange(xmin, xmax, samples)
    ys = Base.LinRange(ymin, ymax, samples)

    # materialize the direction field of the function
    sols = [Float32.(f(x, y)) for x in xs, y in ys]

    # each variable in one array
    us = getindex.(sols, 1)
    vs = getindex.(sols, 2)

    # extract nullclines using Contour.jl
    u_null_curves = AbstractPlotting.ContoursHygiene.Contour.contour(xs, ys, us, zero(xmin)) |> AbstractPlotting.ContoursHygiene.Contour.lines
    v_null_curves = AbstractPlotting.ContoursHygiene.Contour.contour(xs, ys, vs, zero(ymin)) |> AbstractPlotting.ContoursHygiene.Contour.lines

    # make the nullclines single vectors
    u_null_line = map(x -> Point2f0.(x.vertices), u_null_curves) |> GeoMakie.to_nansep_vec
    v_null_line = map(x -> Point2f0.(x.vertices), v_null_curves) |> GeoMakie.to_nansep_vec

    # plot
    lines!(scene, u_null_line; color = wong[5], linewidth = 14)
    lines!(scene, v_null_line; color = wong[6], linewidth = 14)
end

nullcline!(f::Function, domain; kwargs...) = nullcline!(AbstractPlotting.current_scene(), f, domain; kwargs...)

nullcline(f::Function, domain; kwargs...) = nullcline!(Scene(), f, domain; kwargs...)

# some conversion functions for poly
Base.LinRange(interval::AbstractPlotting.IntervalSets.Interval{:closed, :closed, Num}, length::Integer) where Num <: Real = Base.LinRange(interval.left, interval.right, length)

AbstractPlotting.convert_arguments(::AbstractPlotting.PointBased, interval::AbstractPlotting.IntervalSets.Interval{:closed, :closed, Num}, poly::Polynomials.Poly) where Num <: Real = (Point2f0.(Base.LinRange(interval, 500), poly.(Base.LinRange(interval, 500))),)

# allow setting limits for LAxis
function AbstractPlotting.setlims!(scene::LAxis, lims::NTuple{2, Real}, dim=1)
    ol = scene.limits[]                          # get the Scene's limits as values
    o_origin = ol.origin                         # get the original origin
    o_widths = ol.widths                         # get the original widths
    n_widths = convert(Vector, o_widths)         # convert to mutable form
    n_origin = convert(Vector, o_origin)         # convert to mutable form
    n_origin[dim] = lims[1]                      # set the new origin in dim
    n_widths[dim] = lims[2] - lims[1]            # set the new width in dim
    scene.limits[] = AbstractPlotting.HyperRectangle(Vec{2}(n_origin), Vec{2}(n_widths)) # set the limits of the scene
    return nothing
end

AbstractPlotting.xlims!(scene::LAxis, lims::NTuple{2, Real}) = setlims!(scene, lims, 1)

AbstractPlotting.ylims!(scene::LAxis, lims::NTuple{2, Real}) = setlims!(scene, lims, 2)

empty_plots!(la::LAxis) = empty!(plots(la))
AbstractPlotting.plots(la::LAxis) = la.scene.plots
AbstractPlotting.scene_limits(la::LAxis) = la.limits[]

# gather Xu's data

# v_null_raw = readdlm(joinpath(@__DIR__, "vnull.csv"), ',')
# w_null_raw = readdlm(joinpath(@__DIR__, "wnull.csv"), ',')

xu_vnull = Polynomials.Poly([0.4134453340693177, -0.014481137741191953, 0.0007940304464515647, 4.5698911401249415e-5, 5.246507347827663e-7])
xu_wnull = Polynomials.Poly([0.9380140455857346, 0.0030498536445784227, -6.424069411924451e-5, 4.545559542269378e-7])

using MakieLayout
using GeoMakie: WilkinsonTicks # for those sweet sweet wilkinson ticks

# print(collect('𝑢':'\i'))
```

```julia id=7f805704-1b41-4b86-8c4d-b938755ca873

scene, layout = layoutscene(
    resolution = (2000, 1100) .* 4,
    font = "DejaVu Sans",
    # rowsizes = [Auto(), Relative(0.25)]
);

axes = [
    LAxis(
        scene;
        aspect = AxisAspect(1),
        titlesize = 40 * 6,
        titlegap = 20 * 6,
        xticks = WilkinsonTicks(; k_min = 5, k_ideal = 6),
        yticks = WilkinsonTicks(;  k_min = 5, k_ideal = 5),
        xticksvisible = false,
        yticksvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,

    )
    for _ in 1:3
]

layout[1, 1:3] = axes
axes[1].title = "Real 𝐶. 𝐸𝑙𝑒𝑔𝑎𝑛𝑠"
axes[2].title = "Morris-Lecar"
axes[3].title = "FitzHugh-Nagumo"

# nullcline!(axes[1], xu, HyperRectangle(Vec2f0(-100, -1.0), Vec2f0(200, 2)))
nullcline!(axes[2], morris_lecar, AbstractPlotting.HyperRectangle(Vec2f0(-100, -1), Vec2f0(200, 2)))
nullcline!(axes[3], fitzhugh_nagumo, AbstractPlotting.HyperRectangle(Vec2f0(-2, -1), Vec2f0(4, 2)))

xlims!(axes[2], (-75, 75))
ylims!(axes[2], (-0.05, 0.6))

layout[2, 1:3] = leg = LLegend(
    scene,
    plots(axes[2]),
    ["𝑣-nullcline", "𝑤-nullcline"];
    linewidth = 12 * 3,
    titlesize = 0,
    titlevisible = false,
    labelsize = 30 * 6,
    patchsize = (40f0, 40f0) .* 4,
    valign = :top,
    height = Fixed(100 * 8),
    width = Relative(1.2),
    strokecolor = RGBAf0(0, 0, 0, 0),
    rowgap = 100
)

# empty!(plots(axes[1]))

xs = Base.LinRange(-60..0, 500)
lines!(axes[1], xs, xu_vnull.(xs); color = wong[5], linewidth = 14)
lines!(axes[1], xs, xu_wnull.(xs); color = wong[6], linewidth = 14)

ylims!(axes[1], (0.5, 0.9))

scene
```

![result][nextjournal#output#7f805704-1b41-4b86-8c4d-b938755ca873#result]

```julia id=a6088903-193d-48ab-bb52-104a3cd48df4
# Returns a Vector{Point2f0} of points which define the isocline
function isoclines(xs::AbstractVector{<: Real}, ys::AbstractVector{<: Real}, zs::AbstractMatrix{<: Real}, isoval::Real = 0.0)
		return AbstractPlotting.Contours.contour(xs, ys, zs, eltype(xs)(isoval)) |> AbstractPlotting.Contours.lines |> x -> map(y -> Point2f0.(y.vertices), x) |> GeoMakie.to_nansep_vec
end
```

```julia id=83b8e98c-7ddf-4250-b5ac-c1ef80818062
function isoclines(xs::AbstractPlotting.IntervalSets.ClosedInterval{<: Real}, ys::AbstractPlotting.IntervalSets.ClosedInterval{<: Real}, f::Function, isoval::Real = 0.0; samples = 500, ind = 1)
  	xr = LinRange(xs, samples)
    yr = LinRange(ys, samples)
    return isoclines(xr, yr, getindex.([f(x, y) for x in xr, y in yr], ind), isoval)
end
```

```julia id=108fec84-90d8-4d96-bee6-66f257a0fd7c
scene, layout = layoutscene()

la = layout[1, 1] = LAxis(scene; aspect = AxisAspect(1))

lines!(la, isoclines(-2..2, -2..2, fitzhugh_nagumo; ind=1))
lines!(la, isoclines(-2..2, -2..2, fitzhugh_nagumo; ind=2) .+ Ref(Point2f0(-0.5, 0)))

scene
```

![result][nextjournal#output#108fec84-90d8-4d96-bee6-66f257a0fd7c#result]

```julia id=2a1daa0e-a82c-4cc0-be6a-c9989528f9bb
scene, layout = layoutscene(
    resolution = (2000, 1100) .* 4 ./ (3/2, 1) .|> x -> round(Int, x),
    font = "DejaVu Sans",
    # rowsizes = [Auto(), Relative(0.25)]
);

axes = [
    LAxis(
        scene;
        aspect = AxisAspect(1),
        titlesize = 40 * 6,
        titlegap = 20 * 6,
        xticks = WilkinsonTicks(; k_min = 5, k_ideal = 6),
        yticks = WilkinsonTicks(;  k_min = 5, k_ideal = 5),
        xlabel = "𝑣",
        ylabel = "𝑤",
        xlabelsize = 40 * 5,        
        ylabelsize = 40 * 5,
        xticksvisible = false,
        yticksvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,

    )
    for _ in 1:2
]

layout[1, 1:2] = axes

axes[1].title = "Excitable regime"
axes[2].title = "Oscillatory regime"
# axes[3].title = "Unstable regime"

#axes[2].ylabelvisible = false
axes[2].yaxisposition = :right

iso_v = isoclines(-2..2, -2..2, fitzhugh_nagumo; ind=1)
iso_w = isoclines(-2..2, -2..2, fitzhugh_nagumo; ind=2)

# steady-state
lines!(axes[1], iso_v; color = wong[5], linewidth = 14)
lines!(axes[1], iso_w .+ Ref(Point2f0(-0.5, 0)); color = wong[6], linewidth = 14)

# oscillatory
lines!(axes[2], iso_v; color = wong[5], linewidth = 14)
lines!(axes[2], iso_w; color = wong[6], linewidth = 14)

# unstable
#lines!(axes[3], iso_v; color = wong[5], linewidth = 14)
#lines!(axes[3], iso_w .+ Ref(Point2f0(0.5, 0)); color = wong[6], linewidth = 14)



# draw the legend
layout[2, 1:2] = leg = LLegend(
    scene,
    plots(axes[1]),
    ["𝑣-nullcline", "𝑤-nullcline"];
    linewidth = 12 * 3,
    titlesize = 0,
    titlevisible = false,
    labelsize = 30 * 6,
    patchsize = (40f0, 40f0) .* 4,
    valign = :top,
    height = Fixed(100 * 8),
    width = Relative(1.2),
    strokecolor = RGBAf0(0, 0, 0, 0),
    strokewidth = 0,
    framevisible = false,
    rowgap = 100
)

scene
```

![result][nextjournal#output#2a1daa0e-a82c-4cc0-be6a-c9989528f9bb#result]

```julia id=3eb748a7-3d25-4024-8ae0-8b85938481a0
```

```julia id=3d8b480d-644f-4923-ad9a-afbbc2d301ab
```

[nextjournal#output#eec0d8a5-c89a-47e3-9977-44fad72b5505#result]:
<https://nextjournal.com/data/QmSkWk41V5XCPtwm46XMxw5mbaSSPYBowGYx5tJ2fS62V4?content-type=image%2Fpng>

[nextjournal#output#7f805704-1b41-4b86-8c4d-b938755ca873#result]:
<https://nextjournal.com/data/Qmdy1KNzYU3F7XcA5Ve3AtcLNJuTdq4jzjzsmc7g1Msx7W?content-type=image%2Fpng>

[nextjournal#output#108fec84-90d8-4d96-bee6-66f257a0fd7c#result]:
<https://nextjournal.com/data/QmacnrvdhJZB4biWPU5Bd5Py6hyUhAyTwtQpEKvgX3DTho?content-type=image%2Fpng>

[nextjournal#output#2a1daa0e-a82c-4cc0-be6a-c9989528f9bb#result]:
<https://nextjournal.com/data/Qma2VQvr1kZVVjkmMFM48yFhdEgcM6hwuKKdWGQQPQ9xVN?content-type=image%2Fpng>

<details id="com.nextjournal.article">
<summary>This notebook was exported from <a href="https://nextjournal.com/a/MEXVrNSGuMmKnQPvgUYRc?change-id=CeiCEjWw1zaW35TvgZD2Lo">https://nextjournal.com/a/MEXVrNSGuMmKnQPvgUYRc?change-id=CeiCEjWw1zaW35TvgZD2Lo</a></summary>

```edn nextjournal-metadata
{:article
 {:nodes
  {"108fec84-90d8-4d96-bee6-66f257a0fd7c"
   {:compute-ref #uuid "1016fdc7-6b48-4d55-8be1-5e126a925603",
    :exec-duration 2913,
    :id "108fec84-90d8-4d96-bee6-66f257a0fd7c",
    :kind "code",
    :output-log-lines {},
    :runtime [:runtime "289932cd-4e71-41c7-a499-35d8f07f9193"]},
   "13373429-787d-4bb8-9172-4b3d576ffe61"
   {:compute-ref #uuid "5bc955a8-4b3e-4b13-a376-988816f29cb1",
    :exec-duration 513,
    :id "13373429-787d-4bb8-9172-4b3d576ffe61",
    :kind "code",
    :output-log-lines {},
    :runtime [:runtime "289932cd-4e71-41c7-a499-35d8f07f9193"]},
   "289932cd-4e71-41c7-a499-35d8f07f9193"
   {:environment
    [:environment
     {:article/nextjournal.id
      #uuid "5b460d39-8c57-43a6-8b13-e217642b0146",
      :change/nextjournal.id
      #uuid "5e293bbb-3982-45f5-99ed-3aae670e5e6b",
      :node/id "39e3f06d-60bf-4003-ae1a-62e835085aef"}],
    :environment? true,
    :id "289932cd-4e71-41c7-a499-35d8f07f9193",
    :kind "runtime",
    :language "julia",
    :name "Nullclines in Makie",
    :resources
    {:machine-type "n1-standard-4",
     :accelerator-type "nvidia-tesla-k80",
     :accelerator-count 1},
    :type :nextjournal,
    :docker/environment-image
    "docker.nextjournal.com/environment@sha256:4c23168550b2c67aef742cd02700204eab2f85ff2232904d610fbfd8f985bc38"},
   "2a1daa0e-a82c-4cc0-be6a-c9989528f9bb"
   {:compute-ref #uuid "e4fd5080-abd7-4329-914a-bbd53b7c62af",
    :exec-duration 15579,
    :id "2a1daa0e-a82c-4cc0-be6a-c9989528f9bb",
    :kind "code",
    :output-log-lines {},
    :runtime [:runtime "289932cd-4e71-41c7-a499-35d8f07f9193"]},
   "3d8b480d-644f-4923-ad9a-afbbc2d301ab"
   {:id "3d8b480d-644f-4923-ad9a-afbbc2d301ab",
    :kind "code",
    :runtime [:runtime "289932cd-4e71-41c7-a499-35d8f07f9193"]},
   "3eb748a7-3d25-4024-8ae0-8b85938481a0"
   {:compute-ref #uuid "53ece365-4848-47f6-ac19-b31da458cb13",
    :exec-duration 11514,
    :id "3eb748a7-3d25-4024-8ae0-8b85938481a0",
    :kind "code",
    :output-log-lines {},
    :runtime [:runtime "289932cd-4e71-41c7-a499-35d8f07f9193"]},
   "5c61c477-d20c-4c44-95b6-3dadadace86a"
   {:compute-ref #uuid "96a2d295-de7c-47df-b3b4-e48050ee1929",
    :exec-duration 192390,
    :id "5c61c477-d20c-4c44-95b6-3dadadace86a",
    :kind "code",
    :output-log-lines {:stdout 310},
    :runtime [:runtime "289932cd-4e71-41c7-a499-35d8f07f9193"]},
   "7f805704-1b41-4b86-8c4d-b938755ca873"
   {:compute-ref #uuid "7379db0d-1b95-4c26-9b90-7a9de95c6d59",
    :exec-duration 56231,
    :id "7f805704-1b41-4b86-8c4d-b938755ca873",
    :kind "code",
    :output-log-lines {},
    :runtime [:runtime "289932cd-4e71-41c7-a499-35d8f07f9193"]},
   "8012a0b5-d4dd-4276-94cd-07cc296c67c0"
   {:compute-ref #uuid "0c762254-ea6c-477b-8d08-3702e03b574e",
    :exec-duration 56763,
    :id "8012a0b5-d4dd-4276-94cd-07cc296c67c0",
    :kind "code",
    :output-log-lines {},
    :runtime [:runtime "289932cd-4e71-41c7-a499-35d8f07f9193"]},
   "83b8e98c-7ddf-4250-b5ac-c1ef80818062"
   {:compute-ref #uuid "adba9dc4-9520-4774-ba02-1f9dd51c8aef",
    :exec-duration 247,
    :id "83b8e98c-7ddf-4250-b5ac-c1ef80818062",
    :kind "code",
    :output-log-lines {},
    :runtime [:runtime "289932cd-4e71-41c7-a499-35d8f07f9193"]},
   "a6088903-193d-48ab-bb52-104a3cd48df4"
   {:compute-ref #uuid "0d313276-5d40-4681-b7f7-da01a82dbd14",
    :exec-duration 562,
    :id "a6088903-193d-48ab-bb52-104a3cd48df4",
    :kind "code",
    :output-log-lines {},
    :runtime [:runtime "289932cd-4e71-41c7-a499-35d8f07f9193"]},
   "eec0d8a5-c89a-47e3-9977-44fad72b5505"
   {:compute-ref #uuid "673addb6-0575-44d5-980f-b09c39707824",
    :exec-duration 256343,
    :id "eec0d8a5-c89a-47e3-9977-44fad72b5505",
    :kind "code",
    :output-log-lines {:stdout 5},
    :runtime [:runtime "289932cd-4e71-41c7-a499-35d8f07f9193"]}},
  :nextjournal/id #uuid "02d33602-e0a0-414f-9dcc-6291deaca2b7",
  :article/change
  {:nextjournal/id #uuid "5e586e69-3188-46ee-949e-3192620f9920"}}}

```
</details>
