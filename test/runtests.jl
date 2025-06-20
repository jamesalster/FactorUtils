
using FactorUtils
using Test

using RCall: @R_str, @rget, rcopy
using RDatasets
using DataFrames: dropmissing
using Statistics: std, mean
using StatsBase
using Makie: Makie

@info "Setting up test environment..."

# Set up data
bfi_all = dataset("psych", "bfi")[:, 2:26]
bfi = dropmissing(bfi_all)

# Make objects from this package
pca_fu = pca(bfi)
fa_fu = fa(bfi, 5)
fa_fu2 = fa(bfi, 5; scale=false)

@testset "Namespace" begin
    @test loadings === MultivariateStats.loadings
    @test FactorUtils.rotate! === FactorRotations.rotate!
end

@testset "DataProcessing" begin

    #normalisation
    trans_fun = zscore_transform(bfi)
    @test trans_fun isa Function
    bfi_trans = @test_nowarn trans_fun(bfi)
    bfi_trans2 = @test_nowarn trans_fun(bfi[1:5, :])
    @test bfi_trans2 == bfi_trans[1:5, :]
    @test eltype(bfi_trans[:, :A1]) === Float64
    @test isapprox(mean(bfi_trans[:, :A1]), 0.0; atol=1e-10)
    @test isapprox(std(bfi_trans[:, :A1]), 1.0; atol=1e-10)

    new_df = DataFrame(:a => 1:5, :b => 2:6)
    @test_throws ArgumentError trans_fun(new_df)
    new_df2 = DataFrame(:A1 => Float32.(1:5), :C2 => Float32.(2:6))
    new_df2_trans = @test_nowarn trans_fun(new_df2)
    @test eltype(new_df2_trans[:, :A1]) === Float32

    #data processing
    @test_throws ArgumentError prep_data(bfi_all)
    bfi2 = deepcopy(bfi)
    bfi2[:, :newcol] .= "a"
    prepped2 = @test_warn r"newcol" prep_data(bfi2)
    @test size(prepped2, 2) == 25
    bfi2[:, 2] .= 0
    prepped3 = @test_warn r"A2" prep_data(bfi2)
    @test size(prepped3, 2) == 24
end

@testset "Methods" begin
    bfi_mat = mapslices(
        x -> StatsBase.transform(fit(ZScoreTransform, x), x),
        Matrix(convert.(Float64, bfi));
        dims=1,
    )
    # multivariate stats
    fa_ms = fit(FactorAnalysis, bfi_mat', maxoutdim=5, method=:cm)

    # object properties
    @test fa_fu isa FactorResults
    @test fa_fu.fa isa FactorAnalysis
    @test fa_fu.nm == names(bfi)
    @test fa_fu.trans isa Function
    @test all(isapprox.(mean(fa_fu.X; dims=2), 0.0, atol=1e-10))
    @test all(isapprox.(std(fa_fu.X; dims=2), 1.0, atol=1e-10))
    @test size(fa_fu.X) == size(bfi_mat')
    @test fa_fu2.trans == identity

    # methods
    @test size(fa_fu) == size(fa_ms)
    @test loadings(fa_fu) ≈ MultivariateStats.loadings(fa_ms)
    @test loadings(fa_fu) isa DimArray
    @test var(fa_fu) ≈ MultivariateStats.var(fa_ms)
    @test cov(fa_fu) ≈ MultivariateStats.cov(fa_ms)
    @test size(predict(fa_fu)) == (size(bfi, 1), 5)
    @test predict(fa_fu) ≈ MultivariateStats.predict(fa_ms, bfi_mat')' #nb shape
    @test projection(fa_fu) ≈ MultivariateStats.projection(fa_ms)
    @test mean(fa_fu) ≈ MultivariateStats.mean(fa_ms)

    # test prediciton on new data, applying scaling. 
    @test predict(fa_fu, bfi[1:9, :]) ≈ predict(fa_fu)[1:9, :]
    # won't work if not scaled
    @test !isapprox(
        predict(fa_fu, bfi[1:9, :]; apply_scaling=false), predict(fa_fu)[1:9, :]
    )

    # new methods 
    @test cos2_var(fa_fu) ≈ MultivariateStats.loadings(fa_ms) .^ 2
    ## TODO add test for cos2 ind
    @test unique_variance(fa_fu) ≈
        1 .- dropdims(sum(MultivariateStats.loadings(fa_ms) .^ 2; dims=2); dims=2)

    # now for pca
    pca_ms = fit(PCA, bfi_mat')

    # object properties
    @test pca_fu isa FactorResults
    @test pca_fu.fa isa PCA
    @test pca_fu.nm == names(bfi)
    @test all(isapprox.(mean(pca_fu.X; dims=2), 0.0, atol=1e-10))
    @test all(isapprox.(std(pca_fu.X; dims=2), 1.0, atol=1e-10))
    @test size(pca_fu.X) == size(bfi_mat')

    # methods
    @test size(pca_fu) == size(pca_ms)
    @test loadings(pca_fu) ≈ MultivariateStats.loadings(pca_ms)
    @test loadings(pca_fu) isa DimArray
    @test size(predict(pca_fu)) == (size(bfi, 1), 25)
    @test predict(pca_fu) ≈ MultivariateStats.predict(pca_ms, bfi_mat')' #nb shape
    @test projection(pca_fu) ≈ MultivariateStats.projection(pca_ms)
    @test mean(pca_fu) ≈ MultivariateStats.mean(pca_ms)

    # new methods 
    @test cos2_var(pca_fu) ≈ MultivariateStats.loadings(pca_ms) .^ 2
    ## TODO add test for cos2 ind
    @test unique_variance(pca_fu) ≈
        1 .- dropdims(sum(MultivariateStats.loadings(pca_ms) .^ 2; dims=2); dims=2)

    ## Factor rotations methods
    fa_new = deepcopy(fa_fu)

    old_loadings = loadings(fa_fu)
    FactorUtils.rotate!(fa_new, Geomin())
    @test fa_new.fa isa FactorAnalysis
    @test !isapprox(loadings(fa_new) .- old_loadings, zeros(size(fa_fu)), atol=1e-5)
    FactorRotations.rotate!(fa_ms, Geomin())
    @test all(loadings(fa_new) .≈ loadings(fa_ms))

    fa_new = efa(fa_new, Varimax())
    @test fa_new isa FactorResults
    @test fa_new.fa isa FactorAnalysis
end

@testset "Pretty Printing" begin

    # Create test data
    data = [0.8 0.3 -0.1; 0.2 0.9 -0.4; -0.6 0.1 0.7]
    arr = loadings(fa_fu)

    # Test basic functionality
    @test_nowarn pretty(arr)

    # Test with no highlighters
    @test_nowarn pretty(arr; highlighters=())

    # Test error for non-2D arrays
    arr_1d = DimArray([1, 2, 3], Dim{:dim1}(["a", "b", "c"]))
    @test_throws ArgumentError pretty(arr_1d)

    arr_3d = DimArray(
        rand(2, 2, 2),
        (Dim{:dim1}(["a", "b"]), Dim{:Dim2}(["c", "d"]), Dim{:dim3}(["e", "f"])),
    )
    @test_throws ArgumentError pretty(arr_3d)

    # Capture output using redirect
    output = sprint(pretty, arr)

    @test contains(output, "f1")
    @test contains(output, "0.307")
    @test contains(output, "variable / factor")
    @test contains(output, "A1")

    output2 = @test_nowarn sprint((io, x) -> pretty(io, x; topcorner="Custom Label"), arr)
    @test contains(output2, "Custom Label")

    output_fa = @test_nowarn sprint(pretty, fa_fu)
    @test contains(output_fa, "Factor Loadings")
    @test contains(output_fa, "Variable")
    @test contains(output_fa, "C1")
    @test contains(output_fa, "f1")
    @test contains(output_fa, "0.285")
    @test contains(output_fa, "Empirical Latent Variable Correlations")
    @test contains(output_fa, "1.000")

    output_pca = @test_nowarn sprint(pretty, pca_fu)
    @test contains(output_pca, "Factor Loadings")
    @test contains(output_pca, "Variable")
    @test contains(output_pca, "C1")
    @test contains(output_pca, "f1")
    @test contains(output_pca, "0.376")
    @test contains(output_pca, "Variance Explained")
    @test contains(output_pca, "Eigenvalue")
    @test contains(output_pca, "% of Variance")
    @test contains(output_pca, "Cumulative % of Variance")
    @test contains(output_pca, "5.134")
end

@testset "Plots" begin
    @warn "No tests implemented for Plotting functionality. Test indplot() and biplot() manually."
end

@testset "Lavaan Functionality" begin
    # Test R is installed and accessible
    @test_nowarn R"""version"""
    @test rcopy(R"""require(lavaan)""")
    r_version = rcopy(R"R.version.string")
    @test startswith(r_version, "R version 4")
    lav_efa = @test_warn r"lavaan" efa_lavaan(bfi, 5, "geomin")
    r_loadings = let
        R"""
            efa_fit = lavaan::efa($bfi, nfactors = 5, rotation = "geomin")
        """
        efa_fit = @rget efa_fit
        efa_fit[:loadings]
    end
    @test lav_efa isa FactorResults
    @test lav_efa.nm == names(bfi)
    # only approximate
    @test isapprox(abs.(r_loadings), abs.(loadings(lav_efa)); atol=1e-2)
end
