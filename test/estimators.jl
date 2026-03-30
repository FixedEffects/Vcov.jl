@testset "SimpleCovariance" begin
    v = simple()
    @test sprint(show, v) == "Simple covariance estimator"
end

@testset "RobustCovariance" begin
    v = robust()
    @test sprint(show, v) == "Heteroskedasticity-robust (HC1) covariance estimator"
    @test sprint(show, Vcov.hc(1)) == "Heteroskedasticity-robust (HC1) covariance estimator"
    @test sprint(show, Vcov.hc(2)) == "Heteroskedasticity-robust (HC2) covariance estimator"
    @test sprint(show, Vcov.hc(3)) == "Heteroskedasticity-robust (HC3) covariance estimator"
    @test_throws ArgumentError Vcov.hc(4)
end

@testset "ClusterCovariance" begin
    @test_throws MethodError cluster()
    c1 = cluster(:a)
    @test names(c1) == (:a,)
    @test length(c1) == 1
    c2 = cluster(:a, :b)
    @test names(c2) == (:a, :b)
    @test length(c2) == 2

    @test sprint(show, c1) == "Cluster-robust covariance estimator"
    @test sprint(show, MIME("text/plain"), c1) == """
        1-way cluster-robust covariance estimator:
          a"""
    @test sprint(show, MIME("text/plain"), c2) == """
        2-way cluster-robust covariance estimator:
          a
          b"""

    N = 10
    a1 = collect(1:N)
    a2 = convert(Vector{Union{Float64, Missing}}, a1)
    a2[1] = missing
    cols = (a=a1, b=a2)

    g1 = GroupedArray(a1)
    c1 = cluster((:a,), (g1,))
    @test nclusters(c1) == (a=N,)
    c1m = materialize(cols, cluster(:a))
    @test c1m.clusternames == c1.clusternames
    @test c1m.clusters[1] == g1 

    @test completecases(cols, c1) == trues(N)
    c2 = materialize(cols, cluster(:b))
    @test completecases(cols, c2) == (1:N.!=1)
    c3 = materialize(cols, cluster(:a, :b))
    @test completecases(cols, c3) == (1:N.!=1)
end

##############################################################################
## Vcov computation tests
##############################################################################

# Set up a reproducible OLS problem
rng = MersenneTwister(42)
N, K = 100, 3
X = randn(rng, N, K)
e = randn(rng, N)
XtX = Symmetric(X' * X)
invXtX = inv(XtX)

vcovdata = Vcov.VcovData(X, XtX, invXtX, e, N - K)

@testset "SimpleCovariance vcov" begin
    V = vcov(vcovdata, simple())
    expected = Symmetric(invXtX .* sum(abs2, e) / (N - K))
    @test V ≈ expected
end

@testset "SimpleCovariance S_hat" begin
    S = Vcov.S_hat(vcovdata, simple())
    expected = Symmetric(XtX .* sum(abs2, e))
    @test S ≈ expected
end

@testset "RobustCovariance vcov" begin
    V = vcov(vcovdata, robust())
    # HC1: inv(X'X) * X' diag(e^2) X * (N/(N-K)) * inv(X'X)
    meat = X' * Diagonal(e .^ 2) * X .* (N / (N - K))
    expected = Symmetric(invXtX * meat * invXtX)
    @test V ≈ expected
end

@testset "RobustCovariance S_hat" begin
    S = Vcov.S_hat(vcovdata, robust())
    expected = Symmetric(X' * Diagonal(e .^ 2) * X .* (N / (N - K)))
    @test S ≈ expected
end

@testset "HC2 vcov" begin
    V = vcov(vcovdata, Vcov.hc(2))
    h = vec(sum(X .* (X * invXtX), dims=2))
    meat = X' * Diagonal(e .^ 2 ./ (1 .- h)) * X
    expected = Symmetric(invXtX * meat * invXtX)
    @test V ≈ expected
end

@testset "HC3 vcov" begin
    V = vcov(vcovdata, Vcov.hc(3))
    h = vec(sum(X .* (X * invXtX), dims=2))
    meat = X' * Diagonal(e .^ 2 ./ (1 .- h) .^ 2) * X
    expected = Symmetric(invXtX * meat * invXtX)
    @test V ≈ expected
end

@testset "ClusterCovariance vcov" begin
    # 10 clusters of 10 observations each
    groups = repeat(1:10, inner=10)
    g = GroupedArray(groups)
    cv = cluster((:cl,), (g,))
    vcovdata_cl = Vcov.VcovData(X, XtX, invXtX, e, N - K)
    V = vcov(vcovdata_cl, cv)
    @test size(V) == (K, K)
    @test V ≈ V'
    # Singleton clusters: cluster S_hat should match robust S_hat
    g_single = GroupedArray(collect(1:N))
    cv_single = cluster((:cl,), (g_single,))
    S_cl = Vcov.S_hat(vcovdata, cv_single)
    S_rob = Vcov.S_hat(vcovdata, robust())
    @test S_cl ≈ S_rob
end

@testset "ranktest" begin
    # Edge case: empty matrices return NaN
    @test isnan(Vcov.ranktest!(zeros(10, 0), zeros(10, 1), zeros(0, 1), simple(), 0, 0))
    @test isnan(Vcov.ranktest!(zeros(10, 1), zeros(10, 0), zeros(1, 0), simple(), 0, 0))

    # Smoke test with simple covariance: just-identified IV
    rng2 = MersenneTwister(123)
    n = 200
    Z = randn(rng2, n, 1)
    Xendo = Z * [0.8] .+ randn(rng2, n, 1) * 0.5
    Pi = Z \ Xendo
    stat = Vcov.ranktest!(Xendo, Z, Pi, simple(), 0, 0)
    @test isfinite(stat)
    @test stat > 0

    # Smoke test with robust covariance
    stat_r = Vcov.ranktest!(Xendo, Z, Pi, robust(), 0, 0)
    @test isfinite(stat_r)
    @test stat_r > 0
end
