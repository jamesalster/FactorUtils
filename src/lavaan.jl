
"""
    efa_lavaan(df::DataFrame, nfactors::Int, rotation::String; scale=true)

Use the R package `lavaan` to return a FactorResults object, using the different method in
    `lavaan::efa()`. Rotation string should match `lavaan` docs.

N.B. that this places the loadings and means into a julia FactorAnalysis object to allow using MultivariateStats methods.
    However the underlying implementation is different so the outputs of these methods do not match `lavaan`s own outputs.

Note you will need an R installation with lavaan installed, and to import `RCall` into your session.
"""
function efa_lavaan(df::DataFrame, nfactors::Int, rotation::String; scale=true)
    @warn "Factor analysis from lavaan reproduces the loadings but other methods may not match lavaan output."
    df_fit = prep_data(df)
    nm = names(df_fit) #get names
    mat = Matrix(df_fit)
    X = scale ? mapslices(normalise, mat; dims = 1) : mat
    #fit in lavaan
    X_df = DataFrame(X, nm) #lavaan needs df
    RCall.R"""
        require(lavaan)
        efa_fit = efa($X_df, nfactors = $nfactors, rotation = $rotation, output = "lavaan")
        fac = lavInspect(efa_fit, "est")$lambda
        resid = lavInspect(efa_fit, "est")$theta
    """
    fac = RCall.@rget fac; resid = diag(RCall.@rget resid)
    mn = vec(mean(X; dims = 1))
    fa_obj = FactorAnalysis{eltype(fac)}(mn, fac, resid)
    return FactorResults(fa_obj, X', nm)
end