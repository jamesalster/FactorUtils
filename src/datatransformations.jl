
"""
    zscore_transform(df::DataFrame)

Returns a function that z-scores a new DataFrame in the same way as `df`.
    Uses StatsBase.ZScoreTransform. Non-float columns are converted to Float64. 
    Non-numeric columns will error.
"""
function zscore_transform(df::DataFrame)::Function
    transform_dict = Dict()
    for (nm, col) in pairs(eachcol(df))
        col_float = eltype(col) <: AbstractFloat ? col : convert.(Float64, col)
        transform_dict[nm] = fit(ZScoreTransform, col_float)
    end
    function transform_fun(DF::DataFrame)::DataFrame
        new_df = DataFrame()
        names_mismatch = setdiff(string.(names(DF)), string.(keys(transform_dict)))
        length(names_mismatch) > 0 &&
            throw(ArgumentError("Column names do not match: $names_mismatch"))
        for (nm, col) in pairs(eachcol(DF))
            col_float = eltype(col) <: AbstractFloat ? col : convert.(Float64, col)
            new_df[:, nm] = StatsBase.transform(transform_dict[nm], col_float)
        end
        return new_df
    end
    return transform_fun
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
    not_numeric = findall(col -> !<:(eltype(col),Number), eachcol(df))
    !isempty(not_numeric) &&
        @warn """Dropping columns because type is not numeric: $(join(names(df)[not_numeric], "\n"))"""
    #drop if sum is 0
    sum_is_zero = findall(col -> (eltype(col) <: Number) && (sum(col) ≈ 0.0), eachcol(df))
    !isempty(sum_is_zero) &&
        @warn """Dropping columns because sum is 0: $(join(names(df)[sum_is_zero], "\n"))"""
    cols_to_keep = setdiff(1:size(df, 2), union(sum_is_zero, not_numeric))
    #return
    return df_nomissing[:, cols_to_keep]
end
