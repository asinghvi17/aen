using ModelingToolkit

"""
    construct_fhn(vars, params)
vars   : :v, :w
params : :g, :e, :b
"""
function constructFHN(vars, params, D)
    return [
        min(max(-2 - vars[1], vars[1]), 2 - vars[1]) - vars[2],

        params[2] * (vars[1] - params[1] * vars[2] + params[3])
    ]
end

function constructSubscript(i::Int; subscripts = ["₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"])

    return join(subscripts[reverse(digits(i)) .+ 1])

end

function constructCPG(numpairs = 8)

vardim = 2

paramdim = 5

numparams = paramdim * numpairs*2

numvars = numpairs * 2 * vardim

varnames = [:v, :w]

paramnames = [:g, :e, :b, :Dₕ, :Dᵥ]

t = Variable(:t; known = true)()

D = Differential(t)

allvars = [
    Variable(Symbol(varnames[i % vardim + 1], constructSubscript(i÷vardim)))(t)
    for i in vardim:numvars+vardim-1
    ]

allparams = [
    Variable(Symbol(paramnames[i % paramdim + 1], constructSubscript(i÷paramdim)); known = true)()
    for i in paramdim : numparams + paramdim - 1
    ]

# TODO add the actual functionality - building a system of diffeqs!

alleqs = constructFHN(@view(allvars[1:2]), @view(allparams[1:5]), D)

for i in 2:numpairs*2

    append!(alleqs, constructFHN(@view(allvars[(1:vardim) .+ i*vardim .- vardim]), @view(allparams[(1:paramdim) .+ i*paramdim .- paramdim]), D))

    alleqs[i*vardim] += allparams[i*paramdim + ] # TODO FIXME make the system first top then bottom.  therefore


end

alleqs



end

sys1 = constructFHN(@view(allvars[1:2]), @view(allparams[1:5]), D)

sys2 = constructFHN(@view(allvars[3:4]), @view(allparams[6:10]), D)

sys2[1] += allparams[4]*max(allvars[1] - allvars[3], 0) + 0.1

tsys = cat(sys1, sys2; dims=1)

fhns = @. ~(D($@view(allvars[1:4])), tsys)

sys = ODESystem(fhns)

generate_function(sys, allvars[1:4], allparams[1:10])
