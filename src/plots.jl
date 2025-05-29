
"""
    biplot(fa::FactorResults; dims=[1,2], max_labels=10, circle=true, color=nothing, kwargs...)

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
- `kwargs...`: Passed to Figure()
"""
function biplot(fa::FactorResults; dims=[1, 2], kwargs...)
    fig, ax = biplotarrows(fa; dims=dims, kwargs...)
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
    fig
end

"""
    indplot(fa::FactorResults; dims=[1,2], color=nothing, kwargs...)

Plot individual data points in principal component space.

# Arguments
- `coords`: Matrix where each row is an observation and columns are principal components
- `dims`: Which two components to plot (default is PC1 vs PC2)  
- `color`: How to color the points - can be numbers (creates a color scale) or categories
- `kwargs...`: Extra plotting options passed to Figure()
"""
function indplot(fa::FactorResults; dims=[1, 2], kwargs...)
    fig, ax = indscatter(fa; dims=dims, kwargs...)
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
    fig
end
