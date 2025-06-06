
# FactorUtils

Utilities for Factor Analysis in Julia. 

This package is essentially a high-level wrapper around [MultivariateStats.jl](https://juliastats.org/MultivariateStats.jl/dev/) and [FactorRotations.jl](https://p-gw.github.io/FactorRotations.jl/stable/). 

It centres on methods `fa` and `pca` that take `DataFrame`s, and return `FactorResults` objects which store not only the `FactorAnalysis` output from *MultivariateStats* but also the variable names.

These objects then produce output as prettier `DimArrays` from [DimensionalData.jl](https://rafaqz.github.io/DimensionalData.jl/stable/) and can be called with `pretty()` to provide informative output about the model. `biplot()` and `indplot()` methods are also provided for Makie.

NB while *MultivariateStats* uses **variable * observation** matrices, this package uses **observation * variable** matrices.

## Methods exported

See the documentation for arguments.

* **fa**: Factor Analysis on a DataFrame
* **efa**: EFA on the output of fa
* **pca**: Principal Component Analysis on a DataFrame

All methods from MultivariateStats (e.g. loadings, mean, eigvals etc.) work
on the **FactorResults** outputs from those functions.

* **prep_data**: Checks missing-ness and drops non-numeric or sum-to-0 columns
* **normalise**: A wrapper around StatsBase.ZScoreTransform

* **pretty**: Highlgihted pretty output for a DimArray

* **cos2_ind**: get cos2 (contribution) score by individual
* **cos2_var**: get cos2 (contribution) score by variable (= squared loadings)
* **unique_variance**: get unique variance not captured for each variable (= 1 - sum of squared loadings).
* **variance_explained**: get unique variance not captured for each variable (= 1 - sum of squared loadings).

* **biplot**: Quick variable biplot for `FactorResults` object. Loaded with `Makie` if a Makie backend is imported.
* **indplot**: Quick individuals plot for `FactorResults` object. Loaded with `Makie` if a Makie backend is imported.

* **biplotarrows**: Internal plot method for variable biplot. Requires importing a Makie backend
* **indscatter**: quick individual plot on factor axes. Requires importing a Makie backend

* **efa_lavaan**: uses R to take the output from the `lavaan` implementation of EFA. Do `using RCall` first to make this method available.

## Example

```julia
using FactorUtils
using RDatasets

bfi_all = dataset("psych", "bfi")[:,2:26]
# Handle missing data properly, this is just an e.g
bfi = dropmissing(bfi_all)

# PCA
pca1 = pca(bfi, 5)

# Show loadings and variance explained
pretty(pca1)

# Get some statistics
eigvals(pca1)
loadings(pca1) |> pretty
cos2_ind(pca1)
variance_explained(pca1)

# Predict calls on data used for prediction
predict(pca1)
# or can be passed new data as observations * vars, scaling is applied
predict(pca1, bfi[1:10,:])
predict(pca1, bfi[1:10,:]; apply_scaling=false)

# show biplot
using GLMakie
biplot(pca1; max_labels=20)
indplot(pca1)

# show individuals plot and customise
fig, ax = indplot(pca1)
ax.xlabel = "My own label..."

# FA
fa1 = fa(bfi, 5)
# rotate
# preferred to FactorUtils.rotate! since it uses many repeats by default
efa1 = efa(fa1, Geomin()) 
pretty(efa1)
loadings(efa1)
predict(efa1)

# Lavaan efa - note the warning
using RCall
efa2 = efa_lavaan(bfi, 5, "geomin")
pretty(loadings(efa2))
```

## TODO

Add proper tests for Makie plots.

Allow biplot!() and indplot!() methods to work.