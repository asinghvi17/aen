using ModelingToolkit



struct ODESpec

    dimension::Int

    constructor::Function

    names::Tuple{Symbol}

end

struct ODENode

    vars::Tuple{Variable}

    params::Array{Variable, 1}

    eqs::Array{Equation, 1}

end

function buildFitzhughNagumo
