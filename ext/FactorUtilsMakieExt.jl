module FactorUtilsMakieExt

using FactorUtils
using Makie
using MultivariateStats: PCA
import FactorUtils: biplotarrows, biplotarrows!, indscatter, indscatter!, biplot, indplot
import FactorUtils: FactorResults, loadings, predict, variance_explained, PCA

Makie.@recipe(BiPlotArrows, loadings, var_names) do scene
    Makie.Theme(; dims=[1, 2], max_labels=10, circle=true, color=nothing)
end

"""
    biplotarrows(loadings, var_names; dims=[1,2], max_labels=10, circle=true, color=nothing, kwargs...)

Create a biplot with arrows and text showing how variables relate to principal components.

Shows arrows pointing from the origin to each variable's position.
    Longer arrows mean the variable has more influence on those components.
    The unit circle helps you see which variables are strongly related.

# Arguments
- `loadings`: Matrix where each row is a variable and columns are principal components
- `var_names`: Names of your variables
- `dims`: Which two components to show (default shows PC1 vs PC2)
- `max_labels`: Max number of labels to show
- `circle`: Show a circle?
- `color`: Vector of colors, or numbers for a scale, defaults to heatmap on arrow length
- `kwargs...`: Extra plotting options like figure size
"""
function Makie.plot!(bp::BiPlotArrows)
    loadings_ = bp[1][]
    var_names = bp[2][]
    dims = bp.dims[]
    circle = bp.circle[]
    color = bp.color[]
    max_labels = first(bp.max_labels[])

    if circle
        θ = 0:0.01:2π
        Makie.lines!(bp, cos.(θ), sin.(θ); color=:gray, alpha=0.3)
    end

    x, y = loadings_[:, dims[1]], loadings_[:, dims[2]]
    dist = sqrt.(x .^ 2 .+ y .^ 2)
    if isnothing(color)
        Makie.arrows2d!(bp, zeros(length(x)), zeros(length(y)), x, y; color=dist, colormap=:heat)
    elseif eltype(color) <: Number
        Makie.arrows2d!(bp, zeros(length(x)), zeros(length(y)), x, y; color=color, colormap=:heat)
    else
        Makie.arrows2d!(bp, zeros(length(x)), zeros(length(y)), x, y; color=color)
    end

    idx = if length(dist) > max_labels
        partialsortperm(dist, 1:max_labels; rev=true)
    else
        1:length(dist)
    end
    for i in idx
        Makie.text!(bp, x[i], y[i]; text=var_names[i], offset=(0, 5), alpha=dist[i])
    end
    bp
end

"""
    biplotarrows(fa::FactorResults; kwargs...)

Draw biplot arrows and text directly from a fitted factor analysis model.

`kwargs...` are forwarded to the `biplotarrows` recipe — e.g. `dims`, `max_labels`,
`circle`, `color` — plus any Figure-level options (e.g. `figure=(; size=(600,400))`).
"""
biplotarrows(fa::FactorResults; kwargs...) = biplotarrows(loadings(fa), fa.nm; kwargs...)

Makie.@recipe(IndScatter, coords) do scene
    Makie.Theme(; dims=[1, 2], color=nothing)
end

"""
    indscatter(coords; dims=[1,2], color=nothing, kwargs...)

Plot individual data points in principal component space.

# Arguments
- `coords`: Matrix where each row is an observation and columns are principal components
- `dims`: Which two components to plot (default is PC1 vs PC2)
- `color`: How to color the points - can be numbers (creates a color scale) or categories
- `kwargs...`: Extra plotting options
"""
function Makie.plot!(ip::IndScatter)
    coords = ip.coords[]
    dims = ip.dims[]
    color = ip.color[]

    x = collect(coords[:, dims[1]])
    y = collect(coords[:, dims[2]])

    if isnothing(color)
        Makie.scatter!(ip, x, y; alpha=0.3)
    elseif eltype(color) <: Number
        Makie.scatter!(ip, x, y; alpha=0.3, color=color, colormap=:Spectral)
    else
        Makie.scatter!(ip, x, y; alpha=0.3, color=color)
    end
    ip
end

"""
    indscatter(fa::FactorResults; kwargs...)

Plot individuals directly from a factor fit and your original data.

`kwargs...` are forwarded to the `indscatter` recipe — e.g. `dims`, `color` —
plus any Figure-level options.
"""
indscatter(fa::FactorResults; kwargs...) = indscatter(predict(fa); kwargs...)

# Shared axis dressing for biplot/indplot: crosshair lines, equal aspect,
# and axis labels (% variance for PCA, factor number otherwise).
function _decorate_axis!(ax, fa::FactorResults, dims)
    Makie.vlines!(ax, [0]; color=:darkgrey)
    Makie.hlines!(ax, [0]; color=:darkgrey)
    ax.aspect = Makie.DataAspect()
    if typeof(fa) <: FactorResults{<:PCA}
        var_explained = round.(variance_explained(fa) * 100; digits=1)
        ax.xlabel = "PC$(dims[1]) ($(var_explained[dims[1]])%)"
        ax.ylabel = "PC$(dims[2]) ($(var_explained[dims[2]])%)"
    else
        ax.xlabel = "Factor $(dims[1])"
        ax.ylabel = "Factor $(dims[2])"
    end
    ax
end

"""
    biplot(fa::FactorResults; dims=[1,2], kwargs...)

Draw a biplot (arrows + variable labels) directly from a fitted factor model,
with crosshair axis lines, equal aspect ratio, and labeled axes (variance
explained for PCA, factor number otherwise).

`kwargs...` are passed to `Figure()`.
"""
function biplot(fa::FactorResults; dims=[1, 2], kwargs...)
    fig, ax = biplotarrows(fa; dims=dims, kwargs...)
    _decorate_axis!(ax, fa, dims)
    fig
end

"""
    indplot(fa::FactorResults; dims=[1,2], kwargs...)

Plot individual data points from a fitted factor model, with crosshair axis
lines, equal aspect ratio, and labeled axes (variance explained for PCA,
factor number otherwise).

`kwargs...` are passed to `Figure()`.
"""
function indplot(fa::FactorResults; dims=[1, 2], kwargs...)
    fig, ax = indscatter(fa; dims=dims, kwargs...)
    _decorate_axis!(ax, fa, dims)
    fig
end

end