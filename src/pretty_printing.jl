
default_options = (
    formatters=ft_printf("%5.3f"),
    header_crayon=crayon"black bold",
    row_label_header_crayon=crayon"black bold",
    crop=:horizontal,
)

"""
    pretty(io::IO, arr::DimArray; topcorner=nothing, highlighters=factor_highlighters, kwargs...)

Provide a pretty table representation of a DimArray. The default highlighters are those defined by 
    `FactorUtils` for factor analysis, suited for quick understanding of loadings on **normalized** data. 
    Pass `higlighters=()` to see un-highlighted data.

Will throw an error if the array `arr` does not have two dimensions.

## Default highlighters:

- 0.7: strong positive 
- 0.4 to 0.7: positive 
- 0.3 to 0.4: weak positive 
- -0.3 to 0.3: semi-hidden 
- -0.3 to -0.4 weak negative 
- -0.4 to -0.7 negative
- <-0.7 strong negative

## Arguments

- topcorner: The label to print in the top corner. Defaults to the dimnames of the DimArray.
- highlighters: the default highlighters defined by FactorUtils (see above).
- kwargs: passed to `PrettyTables.pretty_table()` Note especially "crop" which is `:horizontal` by default
    but could be `:vertical` or `:none`.
"""
function pretty(
    io::IO, arr::DimArray; topcorner=nothing, highlighters=factor_highlighters, kwargs...
)
    ndims(arr) != 2 && throw(ArgumentError("pretty() only works on 2D Arrays"))
    top_corner = isnothing(topcorner) ? join(label(dims(arr)), " / ") : topcorner
    return pretty_table(
        io,
        parent(arr);
        header=Array(dims(arr, 2)),
        row_labels=Array(dims(arr, 1)),
        row_label_column_title=top_corner,
        highlighters=highlighters,
        default_options...,
        kwargs...,
    )
end

"""
    pretty(io::IO, fa::FactorResults{<:PCA}; nfactors=5, kwargs...)

Provide a pretty representation of a `FactorResults` holding a `PCA`. 
    Shows a loadings table, with unique variance, and the variance explained on each factor.
    Shows only up to `nfactors` factors, default is 5.

Kwargs are passed to `PrettyTables.pretty_table()`.
"""
function pretty(io::IO, fa::FactorResults{<:PCA}; nfactors=5, kwargs...)
    println(io, crayon"bold", "PCA results, showing the first $nfactors dimensions:\n")
    ## Loadings
    loads = loadings(fa)[:, 1:nfactors]
    println(io, crayon"bold", "Factor Loadings")
    pretty(io, loads; topcorner="Variable", kwargs...)
    ## Variance explained
    eigs = eigvals(fa)
    eig_pct = eigs ./ sum(eigs)
    arr = hcat(eigs[1:nfactors], eig_pct[1:nfactors], cumsum(eig_pct[1:nfactors]; dims=1))
    pretty(
        io,
        arr;
        topcorner="Factor",
        highlighters=(),
        header=["Eigenvalue", "% Variance Explained", "Cumulative % of Variance"],
        kwargs...,
    )
end

"""
    pretty(io::IO, fa::FactorResults{<:FactorAnalysis}; kwargs...)

Provide a pretty representation of a `FactorResults` holding a `FactorAnalysis`. 
    Shows a loadings table, with unique variance, and the (empirical) latent
    variable correlations, i.e. the correlation of the predicted individual positions.

Kwargs are passed to `PrettyTables.pretty_table()`.
"""
function pretty(io::IO, fa::FactorResults{<:FactorAnalysis}; kwargs...)
    println(io, crayon"bold", "Factor Analysis results:\n")
    ## Loadings
    loads = set(loadings(fa), Dim{:factor} => DimensionalData.Dimensions.Unordered)
    # broken, messy
    unique_var = unique_variance(fa)
    unique_var = DimArray(
        reshape(unique_var, :, 1), (dims(unique_var, 1), Dim{:factor}(["Unique Var"]))
    )
    arr = cat(loads, unique_var; dims=2)
    println(io, crayon"bold", "Factor Loadings")
    pretty(io, arr; topcorner="Variable", kwargs...)
    ## Latent variable correlations
    println(io, crayon"bold", "Empirical Latent Variable Correlations")
    pretty(io, cor(predict(fa)); topcorner="Correlations", kwargs...)
end

# Catch-all method for non-IO calls
function pretty(fa::FactorResults, args...; kwargs...)
    pretty(stdout, fa, args...; kwargs...)
end
function pretty(da::DimArray, args...; kwargs...)
    pretty(stdout, da, args...; kwargs...)
end
