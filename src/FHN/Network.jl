module Network

    using ParameterizedFunctions

    export Neuron

    struct Neuron
        "The parameters of the neuron"
        params::Array{<:Real, 1} # g, e, b

        "Initial values of the neuron"
        initcond::Array{<:Real, 1} # initial conditions

        "Connection specification (entry is index in array)"
        connections::Array{Int, 1}

        ""
        diffusions::Array{<:Real, 1}

    end

end
