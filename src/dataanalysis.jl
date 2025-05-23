
"""
    normalise(x)

Z-score a numeric vector to give it mean of 0 and standard deviation of 1.
    Uses StatsBase.ZScoreTransform.
"""
function normalise(x)
    x_float = eltype(x) <: AbstractFloat ? x : convert.(Float64, x)
    trans = fit(ZScoreTransform, x_float)
    return StatsBase.transform(trans, x_float)
end

"""
    prep_data(df::DataFrame)::DataFrame

Prepare data for Factor Analysis: take a DataFrame, check there are no missing values, 
    and drop (with warning) columns where the sum is 0, or element type is not numeric.
"""
function prep_data(df::DataFrame)::DataFrame
    #no missing
    df_nomissing = disallowmissing(df)
    #drop if not numeric
    not_numeric = findall(col -> !<:(eltype(col), Number), eachcol(df))
    !isempty(not_numeric) && @warn """Dropping columns because type is not numeric: $(join(names(df)[not_numeric], "\n"))"""
    #drop if sum is 0
    sum_is_zero = findall(col -> (eltype(col) <: Number) && (sum(col) ≈ 0.0), eachcol(df))
    !isempty(sum_is_zero) && @warn """Dropping columns because sum is 0: $(join(names(df)[sum_is_zero], "\n"))"""
    cols_to_keep = setdiff(1:size(df, 2), union(sum_is_zero, not_numeric))
    #return
    return df_nomissing[:,cols_to_keep]
end

"""
    pca(df::DataFrame, nfactors::Int; scale=true)

Perform a principal component analysis, with scaling if scale=true. 
    A wrapper around MultivariateStats' `fit(PCA, X)` to return a `FactorResults` object.

The dataframe should be observations * variables (the transpose of the matrix taken by MultivariateStats)
"""
function pca(df::DataFrame, nfactors::Int; scale=true)::FactorResults
    df_fit = prep_data(df)
    nm = names(df_fit) #get names
    mat = Matrix(df_fit)
    X = scale ? mapslices(normalise, mat; dims = 1) : mat
    pca = fit(PCA, X'; maxoutdim=nfactors)
    fa_obj = FactorResults(pca, X', nm)
    return fa_obj
end

"""
    fa(df::DataFrame, nfactors::Int; scale=true, method = :cm)

Perform a factor analysis, with scaling if scale=true. 
    A wrapper around MultivariateStats' `fit(FactorAnalysis, X)` to return a `FactorResults` object.

The dataframe should be observations * variables (the transpose of the matrix taken by MultivariateStats)

By default the method is `:cm``, which seems more reliable than `:em`. See MultivariateStats.FactorAnalysis for details.
"""
function fa(df::DataFrame, nfactors::Int; scale=true, method=:cm)::FactorResults
    df_fit = prep_data(df)
    nm = names(df_fit) #get names
    mat = Matrix(df_fit)
    X = scale ? mapslices(normalise, mat; dims = 1) : mat
    fa = fit(FactorAnalysis, X'; maxoutdim=nfactors, method=method)
    return FactorResults(fa, X', nm)
end

"""
    efa(fa::FactorResults, rotation::RotationMethod; kwargs...)

Perform a factor rotation on the `FactorResults` object. Defaults are set to a high number of
iterations and random starts. A wrapper to `rotate()`.

`...kwargs` are passed to `FactorRotations.rotate!`
"""
function efa(
    fa::FactorResults, 
    rotation::RotationMethod; 
    kwargs...
)::FactorResults
    return rotate(fa, rotation; maxiter1=10_000, maxiter2=1_000, randomstarts=1_000, kwargs...)
end

