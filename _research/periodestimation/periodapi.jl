using DSP, LombScargle, ChaosTools

function pefunc(alg::Symbol)

end

function pgram(times::AbstractVector, data::AbstractVector{Real}; alg = :auto, verify::Vector{Symbol} = [], kwargs...)

    primary = pefunc(alg)(times, data; kwargs...)

    vfuncs = pefunc.(verify)

    vres = zeros(axes(vfuncs))

    for i ∈ axes(vfuncs)[1]

        vres[i] = vfuncs[i](times, data; kwargs...)

    end

end
