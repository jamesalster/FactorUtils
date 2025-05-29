
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
include("datatransformations.jl")
include("dataanalysis.jl")
include("plot_recipes.jl")

loadings = MultivariateStats.loadings
rotate! = FactorRotations.rotate!
setnames! = NamedArrays.setnames!

export FactorResults,
    cos2_ind,
    cos2_var,
    unique_variance,
    variance_explained,
    pretty,
    zscore_transform,
    prep_data,
    fa,
    pca,
    efa,
    biplotarrows,
    indscatter,
    mean,
    var,
    cov,
    eigvals,
    eigvecs,
    predict,
    r2

function __init__()
    @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" begin
        include("lavaan.jl")
        export efa_lavaan
    end

    #Makie required for vlines and labels for indplot and biplot
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
        include("plots.jl")
        export indplot, biplot
    end
end

end
