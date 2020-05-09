
function eff_signal(
        ts;
        σ = 40,
        kernel = reflect(ImageFiltering.Kernel.gaussian((σ,), (σ * 4 * 2 + 1,))),
    )

    num_neurons = size(ts, 1) ÷ 2 # 2 eqs/neuron
    num_pairs   = num_neurons  ÷ 2 # 2 neurons per pair

    effective_signal = zeros(size(ts, 2), num_pairs)

    for nn in 1:num_pairs
        neuron_number = nn - 1 # offset for indexing
        effective_signal[:, nn]  = (ts[4*neuron_number + 1, :]  .+ abs.(ts[4*neuron_number + 1, :])) ./ 2
        effective_signal[:, nn] -=(ts[4*neuron_number+2 + 1, :] .+ abs.(ts[4*neuron_number+2 + 1, :])) ./ 2
        effective_signal[:, nn]  = imfilter(effective_signal[:, nn], kernel)
    end

    return effective_signal
end

function rotmat(θ)
    return AbstractPlotting.Mat2(cos(θ), sin(θ), -sin(θ), cos(θ))
end

function signal_to_worm(
        Δθ::AbstractVector;
        tightness_of_bend = 0.5,
    )

    xpos = 0
    ypos = 0

    xcoords = Float64[xpos]
    ycoords = Float64[ypos]

    θ = 0

    number_of_joints = length(Δθ)

    for j in 1:number_of_joints
        δθ = Δθ[j]
        θ += tightness_of_bend * δθ

        xpos += cos(θ)
        ypos += sin(θ)
        append!(xcoords, xpos)
        append!(ycoords, ypos)
    end
        # print(i,deltatheta, xpos, ypos)
        # @show θ

    # spline fit
    nodes=range(0, length = number_of_joints+1)

    csx = CubicSplineInterpolation(nodes,xcoords)
    csy = CubicSplineInterpolation(nodes,ycoords)

    worm_coords = range(0, stop = number_of_joints, step = 0.01)
    # I am incrementing the position in steps of
    # 0.01*(distance between nodes)

    smoothwormx=csx.(worm_coords)
    smoothwormy=csy.(worm_coords)

    points = Point2f0.(smoothwormx, smoothwormy)


    rotation = rotmat(-tightness_of_bend * Δθ[1])

    return Ref(rotation) .* points
end
