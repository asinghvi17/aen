module Network

    using LightGraphs, MetaGraphs

    struct FormulaSpec

        template::FormulaTemplate

        coupling::CouplingTemplate

    end

    struct DESpec

        lhs::Vector{Symbol}

        rhs::Vector{FormulaSpec}

        cvars::Vector{Int}

    end

    macro coupling_def(arg) end

end
