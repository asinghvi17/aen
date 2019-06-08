"""
    find_period(data; atol = 0.005, ptol = 5*atol)

Uses the Hilbert transform to determine
the periodicity of the dataset.
This function assumes that the data is known to be
periodic.

Returns the index of the first entry that is 2π away from the initial entry.

Arguments
***

- atol - The error tolerance for the phase difference between indices.
- ptol - The error tolerance for what constitutes "one rotation" of the tangent vector.

!!!note
Currently implemented only for one dimension.
"""
function find_period(data; atol = 0.005, ptol = 5*atol)

    # the function expects a 2d array but we're giving it a 1d timeseries
    H = Hilbert.hilbert(data[:, :])

    ϕ = arg.(H)

    ϕ₀ = ϕ[1]

    Δϕ = 0

    for i in eachindex(ϕ)[2:end]

        Δϕ += abs(ϕ[i] - ϕ[i-1])

        # if the tangent vector has completed a full rotation
        if isapprox(Δϕ, 4π; atol = ptol) # why 4π?

            println(i) # debuggingg info

            println(Δϕ) # debugging info

            isapprox(ϕ[i], ϕ₀; atol = atol) && return i # second check - phase angle should be accurate...

        end

    end

    return -1

end
