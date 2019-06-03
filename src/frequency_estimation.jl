using Hilbert

"""
    arg(z::Complex)

Computes the "angle" of the given complex number
as a vector away from the real axis, counterclockwise.
"""
arg(z::Complex) = atan(z.im, z.re)

"""
    find_period(data)

Uses the Hilbert transform to determine
the periodicity of the dataset.
This function assumes that the data is known to be
periodic.
"""
function find_period(data)

    # the function expects a 2d array but we're giving it a 1d timeseries
    H = Hilbert.hilbert(data[:, :])

    ϕ₀ = arg(H[1]) # the initial angle

    for i in eachindex(H)[2:end]

        # if the tangent vector has completed a full rotation
        arg(H[i]) ≈ ϕ₀ && return i;

    end

    return -1

end

find_period(sin.(0:0.01:4))


data = sin.(0:π/10:3π)

H = Hilbert.hilbert(data[:, :]) # the function expects a 2d array but we're giving it a 1d timeseries

ϕ₀ = arg(H[1]) # the initial angle

arg(H[21]) ≈ ϕ₀

using Plots

plot(arg.(H[:]), data)

for i in eachindex(H)[2:end]
    # if the tangent vector has completed a full rotation
    arg(H[i]) ≈ ϕ₀ && print(i);
end
