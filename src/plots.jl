
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
MakieCore.@recipe(BiPlotArrows, loadings, var_names) do scene
    MakieCore.Theme(
        dims = [1, 2],
        max_labels = 10,
        circle = true,
        color = nothing
    )
end

function MakieCore.plot!(bp::BiPlotArrows)
    loadings = bp[1][]  # Extract value from observable
    var_names = bp[2][]  # Extract value from observable
    dims = bp.dims[]     # Extract value from observable
    circle = bp.circle[]
    color = bp.color[]
    max_labels = first(bp.max_labels[])
    
    # Unit circle
    if circle
        θ = 0:0.01:2π
        MakieCore.lines!(bp, cos.(θ), sin.(θ), color=:gray, alpha=0.3)
    end
    
    # Arrows and labels
    x, y = loadings[:, dims[1]], loadings[:, dims[2]]
    dist = sqrt.(x.^2 .+ y.^2)
    if isnothing(color)
        MakieCore.arrows!(bp, zeros(length(x)), zeros(length(y)), x, y, color=dist, colormap=:heat)
    elseif eltype(color) <: Number
        MakieCore.arrows!(bp, zeros(length(x)), zeros(length(y)), x, y, color=color, colormap=:heat)
    else
        MakieCore.arrows!(bp, zeros(length(x)), zeros(length(y)), x, y, color=color)
    end
    
    #crop text to only most important vars if too many
    idx = length(dist) > max_labels ? partialsortperm(dist, 1:max_labels, rev=true) : 1:length(dist)
    for i in idx
        MakieCore.text!(bp, x[i], y[i], text=var_names[i], offset=(0, 5), alpha = dist[i])
    end
    bp
end

"""
    biplotarrows(fa::FactorResults; kwargs...)

Draw biplot arrows and text directly from a fitted factor analysis model.

Takes your already-fitted model and variable names, then makes the biplot for you.
"""
biplotarrows(fa::FactorResults; kwargs...) = biplotarrows(loadings(fa), fa.nm; kwargs...)

"""
    indplot(coords; dims=[1,2], color=nothing, kwargs...)

Plot individual data points in principal component space.

# Arguments
- `coords`: Matrix where each row is an observation and columns are principal components
- `dims`: Which two components to plot (default is PC1 vs PC2)  
- `color`: How to color the points - can be numbers (creates a color scale) or categories
- `kwargs...`: Extra plotting options
"""
MakieCore.@recipe(IndScatter, coords) do scene
    MakieCore.Theme(
        dims = [1, 2],
        color = nothing
    )
end

function MakieCore.plot!(ip::IndScatter)
    coords = ip.coords[]  # Extract value from observable
    dims = ip.dims[]   # Extract value from observable
    color = ip.color[] # Extract value from observable

    # Add dark grey lines
    # Not possible to have hlines or vlines without importing Makie
    
    if isnothing(color)
        MakieCore.scatter!(ip, coords[:, dims[1]], coords[:, dims[2]], alpha = 0.3)
    elseif eltype(color) <: Number
        MakieCore.scatter!(ip, coords[:, dims[1]], coords[:, dims[2]], alpha = 0.3, color=color, colormap = :Spectral)
        # Note: Colorbar would need to be handled at the figure level
    else
        MakieCore.scatter!(ip, coords[:, dims[1]], coords[:, dims[2]], alpha = 0.3, color=color)
    end
    ip
end

"""
    indscatter(fa::FactorResults; kwargs...)

Plot individuals directly from a factor fit and your original data.

Takes your fitted factors and original data, then plots it.
"""
indscatter(fa::FactorResults; kwargs...) = indscatter(predict(fa); kwargs...)