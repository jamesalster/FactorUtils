
"""
    pca(df::DataFrame; scale=true)

Perform a principal component analysis, with scaling if scale=true. 
    A wrapper around MultivariateStats' `fit(PCA, X)` to return a `FactorResults` object.

The dataframe should be observations * variables (the transpose of the matrix taken by MultivariateStats)
"""
function pca(df::DataFrame; scale=true)::FactorResults
    df_fit = prep_data(df)
    transform_fun = scale ? zscore_transform(df_fit) : identity
    df_scaled = transform_fun(df_fit)
    nm = names(df_scaled) #get names
    X = Matrix(df_scaled)
    pca = fit(PCA, X')
    fa_obj = FactorResults(pca, X', nm, transform_fun)
    return fa_obj
end

"""
    fa(df::DataFrame, nfactors::Int; scale=true, method = :cm)

Perform a factor analysis, with scaling if scale=true. 
    A wrapper around MultivariateStats' `fit(FactorAnalysis, X)` to return a `FactorResults` object.

The dataframe should be observations * variables (the transpose of the matrix taken by MultivariateStats)

By default the method is `:cm`, which seems more reliable than `:em`. See MultivariateStats.FactorAnalysis for details.
"""
function fa(df::DataFrame, nfactors::Int; scale=true, method=:cm)::FactorResults
    df_fit = prep_data(df)
    transform_fun = scale ? zscore_transform(df_fit) : identity
    df_scaled = transform_fun(df_fit)
    nm = names(df_scaled) #get names
    X = Matrix(df_scaled)
    fa = fit(FactorAnalysis, X'; maxoutdim=nfactors, method=method)
    return FactorResults(fa, X', nm, transform_fun)
end

"""
    efa(fa::FactorResults, rotation::RotationMethod; kwargs...)

Perform a factor rotation on the `FactorResults` object. Defaults are set to a high number of
iterations and random starts. A wrapper to `rotate()`.

`...kwargs` are passed to `FactorRotations.rotate!`
"""
function efa(fa::FactorResults, rotation::RotationMethod; kwargs...)::FactorResults
    return rotate(
        fa, rotation; maxiter1=10_000, maxiter2=1_000, randomstarts=1_000, kwargs...
    )
end
