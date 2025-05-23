
module FactorUtils

using Requires
using Reexport
@reexport using MultivariateStats
@reexport using FactorRotations
@reexport using NamedArrays
using PrettyTables
using StatsBase
using DataFrames: DataFrame, disallowmissing
using MakieCore

#mostly for method extensions
import Statistics: mean, var, cov
import LinearAlgebra: eigvals, eigvecs, diag
import StatsAPI: predict, r2

include("factorresults.jl")
include("highlighters.jl")
include("pretty_printing.jl")
include("dataanalysis.jl")
include("plots.jl")

loadings = MultivariateStats.loadings
rotate! = FactorRotations.rotate!
setnames! = NamedArrays.setnames!

export 
    FactorResults,
    cos2_ind, cos2_var, unique_variance,
    pretty,
    normalise, prep_data,
    fa, pca, efa,
    biplotarrows, indscatter,
    mean, var, cov, eigvals, eigvecs, predict, r2

function __init__()
    @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" begin
        include("lavaan.jl")
        export efa_lavaan
    end

    #Makie required for vlines and labels for indplot and biplot
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin

        function get_var_explained(fa::FactorResults{<:PCA}, dims=[1,2])
            eigs = eigvals(fa)
            var_exp = eigs ./ sum(eigs)
            return var_exp[dims]
        end

        function biplot(fa::FactorResults; dims = [1,2], kwargs...)
            fig, ax = biplotarrows(fa; dims = dims, kwargs...)
            Makie.vlines!(ax, [0]; color=:darkgrey)
            Makie.hlines!(ax, [0]; color=:darkgrey)
            ax.aspect = Makie.DataAspect()
            if typeof(fa) <: FactorResults{<:PCA}
                var_explained = round.(get_var_explained(fa) * 100; digits = 1)
                ax.xlabel = "PC$(dims[1]) ($(var_explained[dims[1]])%)"
                ax.ylabel = "PC$(dims[2]) ($(var_explained[dims[2]])%)"
            else
                ax.xlabel = "Factor $(dims[1])"
                ax.ylabel = "Factor $(dims[2])"
            end
            fig
        end

        function indplot(fa::FactorResults; dims = [1,2], kwargs...)
            fig, ax = indscatter(fa; dims = dims, kwargs...)
            Makie.vlines!(ax, [0]; color=:darkgrey)
            Makie.hlines!(ax, [0]; color=:darkgrey)
            ax.aspect = Makie.DataAspect()
            if typeof(fa) <: FactorResults{<:PCA}
                var_explained = round.(get_var_explained(fa) * 100; digits = 1)
                ax.xlabel = "PC$(dims[1]) ($(var_explained[dims[1]])%)"
                ax.ylabel = "PC$(dims[2]) ($(var_explained[dims[2]])%)"
            else
                ax.xlabel = "Factor $(dims[1])"
                ax.ylabel = "Factor $(dims[2])"
            end
            fig
        end

        export indplot, biplot
    end
end

end
