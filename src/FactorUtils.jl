
module FactorUtils

using Requires
using Reexport
using MultivariateStats
using FactorRotations
@reexport using DimensionalData
using DimensionalData.Dimensions: dims, label
using PrettyTables
using StatsBase
using DataFrames: DataFrame, disallowmissing

#mostly for method extensions
import Statistics: mean, var, cov
import LinearAlgebra: eigvals, eigvecs
import StatsAPI: predict, r2
import MultivariateStats: loadings, projection, reconstruct,
    principalvars, tprincipalvar, tresidualvar
import FactorRotations: rotate, rotate!

include("factorresults.jl")
include("highlighters.jl")
include("pretty_printing.jl")
include("datatransformations.jl")
include("dataanalysis.jl")
include("plots.jl")

export FactorResults,
    cos2_ind,
    cos2_var,
    unique_variance,
    variance_explained,
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
    r2,
    loadings,
    rotate,
    rotate!,
    reconstruct,
    projection,
    principalvars,
    tprincipalvar,
    tresidualvar,
    Geomin,
    Quartimin,
    Quartimax,
    Varimax,
    Promax,
    Oblimin,
    Infomax


function __init__()
    @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" begin
        include("lavaan.jl")
        export efa_lavaan
    end
end

end
