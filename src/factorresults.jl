
"""
    FactorResults

Type to hold data, results and variable names from a factor analysis, to allow 
    using these to produce named and labelled outputs.

Implements all methods from MultivariateStats for types `FactorAnalysis` and `PCA` as well
    as `rotate` and `rotate!` from FactorRotations. It also has `cos2_ind` and `cos2_var` methods 
    for individul and variable importance scores, and `unique_variance`.

## Fields

- fa::MultivariateStats.AbstractDimensionalityReduction: FactorAnalysis object from MultivariateStats 
- X::AbstractArray{<:Real}: the input data with dimensions **obseravtions x variables** (NB)
- nm::Vector{String}: variable names
- trans::Function: function that transforms a dataframe by the same scaling transformation used on the original data.
"""
mutable struct FactorResults{T<:MultivariateStats.AbstractDimensionalityReduction}
    fa::T
    X::AbstractArray{<:Real}
    nm::Vector{String}
    trans::Function
end

#extend FactorAnalysis and PCA methods
Base.size(fa::FactorResults) = size(fa.fa)
Base.size(fa::FactorResults, i) = size(fa.fa, i)

function Base.show(io::IO, obj::FactorResults)
    n_vars = size(obj.X, 1)
    n_obs = size(obj.X, 2)
    n_factors = size(obj, 2)

    print(io, "FactorResults with $n_factors factors")
    print(io, "\n  Variables: $n_vars")
    print(io, "\n  Observations: $n_obs")
    print(io, "\n  Method: $(typeof(obj.fa))")
end

function FactorRotations.rotate(
    fa::FactorResults{T}, rot::RotationMethod; kwargs...
) where {T}
    fa_new = deepcopy(fa)
    FactorRotations.rotate!(fa_new.fa, rot; kwargs...)
    return fa_new
end

function FactorRotations.rotate!(fa::FactorResults, rot::RotationMethod; kwargs...)
    FactorRotations.rotate!(fa.fa, rot; kwargs...)
    return fa
end

function MultivariateStats.loadings(fa::FactorResults)
    loads = MultivariateStats.loadings(fa.fa)
    return NamedArray(
        loads; dimnames=(:variable, :factor), names=(fa.nm, ["f$i" for i in 1:size(fa, 2)])
    )
end

function MultivariateStats.projection(fa::FactorResults)
    proj = MultivariateStats.projection(fa.fa)
    return NamedArray(
        proj; dimnames=(:variable, :factor), names=(fa.nm, ["f$i" for i in 1:size(fa, 2)])
    )
end

""" 
    predict(fa::FactorResults, X; apply_scaling=true)
    predict(fa::FactorResults)

Provide the position of the rows in X according to the factor representation in `fa`.

This re-applies the original data normalization using `fa.trans()` if was set to `true` in 
    the original `fa()` or `pca()` call.

NB that the dimensionality expected is observations * variables, unlike 
`MultivariateStats.predict(::FactorAnalysis, X)` which expects variables * observations.

`predict()` without the array X takes the fit data for `fa` as the points to represent.
"""
function predict(fa::FactorResults, X; apply_scaling=true)

    #Check dims
    size(X, 2) != size(fa, 1) &&
        throw(DimensionMismatch("X should be of dimension obs * vars for FactorUtils"))

    # Convert to df if necessary for transformation
    X_df = X isa DataFrame ? X : DataFrame(X, fa.nm) # assume names are correct
    #Transform if necessary and back to matrix
    X_trans = apply_scaling ? Matrix(fa.trans(X_df)) : Matrix(X_df)

    preds = predict(fa.fa, X_trans')'
    return NamedArray(
        Matrix(preds);
        dimnames=(:row, :factor),
        names=(1:size(preds, 1), ["f$i" for i in 1:size(fa, 2)]),
    )
end
predict(fa::FactorResults) = predict(fa, fa.X'; apply_scaling=false)

function reconstruct(fa::FactorResults, z)
    size(z, 2) != size(fa, 2) &&
        throw(DimensionMismatch("z should be of dimension obs * vars for FactorUtils"))
    recon = MultivariateStats.reconstruct(fa.fa, z')'
    return NamedArray(
        Matrix(recon); dimnames=(:row, :factor), names=(1:size(recon, 1), fa.nm)
    )
end

function cov(fa::FactorResults)
    covs = MultivariateStats.cov(fa.fa)
    return NamedArray(covs; dimnames=(:x, :y), names=(fa.nm, fa.nm))
end

function var(fa::FactorResults)
    vec = MultivariateStats.var(fa.fa)
    return NamedArray(vec; dimnames=(:variable,), names=(fa.nm,))
end

function mean(fa::FactorResults)
    vec = MultivariateStats.mean(fa.fa)
    return NamedArray(vec; dimnames=(:variable,), names=(fa.nm,))
end

function eigvals(fa::FactorResults{<:PCA})
    eigs = eigvals(fa.fa)
    return NamedArray(eigs; dimnames=(:factor,), names=(["f$i" for i in 1:size(fa, 2)],))
end
function eigvecs(fa::FactorResults{<:PCA})
    eigs = eigvecs(fa.fa)
    return NamedArray(
        eigs; dimnames=(:variable, :factor), names=(fa.nm, ["f$i" for i in 1:size(fa, 2)])
    )
end
MultivariateStats.principalvars(fa::FactorResults{<:PCA}) = eigvals(fa)
MultivariateStats.tprincipalvar(fa::FactorResults{<:PCA}) = tprincipalvar(fa.fa)
MultivariateStats.tresidualvar(fa::FactorResults{<:PCA}) = tresidualvar(fa.fa)
r2(fa::FactorResults) = r2(fa.fa)

# Add cos2 methods
"""
    cos2_ind(fa::FactorResults)

Get cos2 scores for individuals for a factor analysis.
"""
function cos2_ind(fa::FactorResults)
    factor_scores = predict(fa)
    squared_scores = factor_scores .^ 2
    squared_distances = sum(squared_scores; dims=1)  # 1×n_obs vector
    cos2 = squared_scores ./ squared_distances
    return NamedArray(
        Matrix(cos2);
        dimnames=(:row, :factor),
        names=(1:size(cos2, 1), ["f$i" for i in 1:size(fa, 2)]),
    )
end

"""
    cos2_var(fa::FactorResults)

Get cos2 scores for variables for a factor analysis.
"""
function cos2_var(fa::FactorResults)
    return loadings(fa) .^ 2
end

"""
    unique_variance(fa::FactorResults; nfactors=size(fa, 2))

Get unique variance for variables for a factor analysis across the first nfactors.

This is 1 - the sum of the squared loadings across factors.
"""
function unique_variance(fa::FactorResults; nfactors=size(fa, 2))
    unique_var = vec(1 .- sum(loadings(fa)[:, 1:nfactors] .^ 2; dims=2))
    return NamedArray(unique_var; dimnames=(:variable,), names=(fa.nm,))
end

"""
    variance_explained(fa::FactorResults{<:PCA})

Get variance explained by PCA dimension (= eigvals ./ sum(eigvals))
"""
function variance_explained(fa::FactorResults{<:PCA})
    eigs = eigvals(fa)
    return eigs ./ sum(eigs)
end
