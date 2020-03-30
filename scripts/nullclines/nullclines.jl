using Polynomials
using MakieLayout, Makie
import AbstractPlotting: xlims!, ylims!, setlims!
using GeoMakie

using DelimitedFiles
# backend
# using CairoMakie
# get a colour palette
wong = AbstractPlotting.wong_colors

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

function fitzhugh_nagumo(v, w; i_ext = 0.0)
    e, g, b = (0.04, 0.8, 0) # unpack parameters
    return [
        w - (-v^3/3 + v) + i_ext,    # dv₀
        e * (v - g * w + b)          # dw₀
    ]
end

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

## isocline utilities

# Returns a Vector{Point2f0} of points which define the isocline
function isoclines(xs::AbstractVector{<: Real}, ys::AbstractVector{<: Real}, zs::AbstractMatrix{<: Real}, isoval::Real = 0.0)
		return AbstractPlotting.Contours.contour(xs, ys, zs, eltype(xs)(isoval)) |> AbstractPlotting.Contours.lines |> x -> map(y -> Point2f0.(y.vertices), x) |> GeoMakie.to_nansep_vec
end

# extract isoclines from a multivariate function
function isoclines(xs::AbstractPlotting.IntervalSets.ClosedInterval{<: Real}, ys::AbstractPlotting.IntervalSets.ClosedInterval{<: Real}, f::Function, isoval::Real = 0.0; samples = 500, ind = 1)
  	xr = LinRange(xs, samples)
    yr = LinRange(ys, samples)
    return isoclines(xr, yr, getindex.([f(x, y) for x in xr, y in yr], ind), isoval)
end

##################
# arrows overloads

function AbstractPlotting.convert_arguments(::Type{<: Arrows}, x::AbstractVector, y::AbstractVector, u::Function)
    u_out = u.(x, y')
    return (vec(Point2f0.(x, y')), vec(Vec2f0.(getindex.(u_out, 1), getindex.(u_out, 2))))
end

function MakieLayout.legendelements(plot::Arrows)
    MakieLayout.LegendElement[ArrowElement(
        linecolor = plot.linecolor, arrowhead = plot.arrowhead,
        linestyle = plot.linestyle, arrowcolor = plot.arrowcolor)]
end

struct ArrowElement <: MakieLayout.LegendElement
    attributes::Attributes
end

ArrowElement(;kwargs...) = ArrowElement(Attributes(kwargs...))

# function MakieLayout.legendsymbol!(scene, element::ArrowElement, bbox, defaultattrs::Attributes)
#     merge!(element.attributes, defaultattrs)
#     attrs = element.attributes
#
# 	arrowsize = get(attrs, :arrowsize, lift(x -> x * 3.5, attrs.linewidth))
#
#     fracpoints = attrs.linepoints
#
#     points = @lift(MakieLayout.fractionpoint.(Ref($bbox), $fracpoints))
#
#     origin = @lift([$points[1]])
# 	diff_bbox = @lift(diff($points))
#     widths = @lift([$diff_bbox[1] .- Point2f0($arrowsize / 40 * $diff_bbox[1])]) # get the arrow's scale
#
#
#     return arrows!(
#         scene,
#         origin, widths;
#         linecolor = attrs.linecolor, arrowhead = attrs.arrowhead,
#         linestyle = attrs.linestyle, arrowcolor = attrs.arrowcolor,
#         arrowsize = arrowsize, raw = true
#     )[end]
# end


# gather Xu's data

# v_null_raw = readdlm(joinpath(@__DIR__, "vnull.csv"), ',')
# w_null_raw = readdlm(joinpath(@__DIR__, "wnull.csv"), ',')

xu_vnull = Polynomials.Poly([0.4134453340693177, -0.014481137741191953, 0.0007940304464515647, 4.5698911401249415e-5, 5.246507347827663e-7])
xu_wnull = Polynomials.Poly([0.9380140455857346, 0.0030498536445784227, -6.424069411924451e-5, 4.545559542269378e-7])


# unit stuff
# struct Unit{sym} end
# Base.:*(x, @nospecialize(::Unit{Sym})) where Sym = (x, sym)
#
# const px = Unit{:px}()
#
# 10px # (10, :px)

# wong_palette = RGB[
#     RGB{N0f8}((( 86, 180, 233) ./ 255)...), # sky blue
#     RGB{N0f8}(((230, 159,   0) ./ 255)...), # orange
#     RGB{N0f8}(((  0, 158, 115) ./ 255)...), # blueish green
#     RGB{N0f8}(((240, 228,  66) ./ 255)...), # yellow
#     RGB{N0f8}(((  0, 114, 178) ./ 255)...), # blue
#     RGB{N0f8}(((213,  94,   0) ./ 255)...), # vermillion
#     RGB{N0f8}(((204, 121, 167) ./ 255)...), # reddish purple
#     ]
#
# function print2int(c::RGB)
#     r = Int(red(c) * 255)
#     g = Int(green(c) * 255)
#     b = Int(blue(c) * 255)
#
#     println("{rgb,255:red,$r:green,$g:blue,$b}")
# end
