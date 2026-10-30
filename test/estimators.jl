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

@testset "pinvertible" begin
    # a tiny variance is not a negative eigenvalue
    A = Symmetric([1e-20 1e-11; 1e-11 1.0])
    @test Vcov.pinvertible(A) === A
    # a singular matrix is left alone
    A = Symmetric([1.0 1.0; 1.0 1.0])
    @test Vcov.pinvertible(A) === A
    A = Symmetric([1.0 0.0; 0.0 0.0])
    @test Vcov.pinvertible(A) === A
    # negative eigenvalues are set to zero
    A = Symmetric([1.0 2.0; 2.0 1.0])
    @test (@test_logs (:warn,) Vcov.pinvertible(A)) ≈ [1.5 1.5; 1.5 1.5]
    A = Symmetric([1.0 0.0; 0.0 -1e-20])
    @test (@test_logs (:warn,) Vcov.pinvertible(A)) ≈ [1.0 0.0; 0.0 0.0]
    # standard errors are equivariant to the scale of a regressor
    rng = MersenneTwister(1)
    n = 1000
    X = hcat(ones(n), randn(rng, n), randn(rng, n))
    res = randn(rng, n) .* (1 .+ abs.(X[:, 2]))
    g = rand(rng, 1:50, n)
    function se(s)
        Xs = X .* [1 1 s]
        XX = Symmetric(Xs' * Xs)
        data = Vcov.VcovData(Xs, XX, Symmetric(inv(XX)), res, n - 3)
        return [sqrt.(diag(vcov(data, v))) for v in (robust(), materialize((g = g,), cluster(:g)))]
    end
    for s in (1e6, 1e9, 1e11)
        for (a, b) in zip(se(1.0), se(s))
            @test b ≈ a ./ [1, 1, s] rtol = 1e-6
        end
    end
end

@testset "pinvertible with a non positive semi-definite matrix" begin
    # the adjustment does not depend on the units of the regressors
    A = Symmetric([1.0 2.0; 2.0 1.0])
    d = [1.0, 1e8]
    B = Symmetric(A ./ (d * d'))
    a = @test_logs (:warn,) Vcov.pinvertible(A)
    b = @test_logs (:warn,) Vcov.pinvertible(B, d)
    @test b ≈ a ./ (d * d')
    # its scale is the norm of each regressor, which is 1 for a zero column
    X = [3.0 0.0; 4.0 0.0]
    XX = Symmetric(X' * X)
    @test Vcov.regressor_scale(Vcov.VcovData(X, XX, XX, zeros(2), 0)) == [5.0, 1.0]
end

@testset "single cluster" begin
    n = 20
    X = hcat(ones(n), collect(1.0:n))
    XX = Symmetric(X' * X)
    data = Vcov.VcovData(X, XX, Symmetric(inv(XX)), randn(n), n - 2)
    @test all(isnan, diag(vcov(data, materialize((g = ones(Int, n),), cluster(:g)))))
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

    # Several endogenous regressors. With simple covariance the statistic equals n * s_k^2, where
    # s_k is the smallest singular value of theta (the Cragg-Donald form), and it must not depend
    # on the order of the endogenous variables
    rng3 = MersenneTwister(42)
    n = 500
    for (k, l) in ((3, 3), (2, 4))
        Z = randn(rng3, n, l)
        Xendo = Z * randn(rng3, l, k) .+ randn(rng3, n, k) * randn(rng3, k, k)
        Pi = Z \ Xendo
        X = Xendo - Z * Pi    # first-stage residuals, as passed by FixedEffectModels
        Fm = cholesky(Symmetric(Z' * Z)).U
        Gm = cholesky(Symmetric(X' * X)).U
        s_k = minimum(svdvals(Fm * (Gm' \ Pi')'))
        @test Vcov.ranktest!(copy(X), copy(Z), copy(Pi), simple(), 0, 0) ≈ n * s_k^2
        # the degrees of freedom enter as in the Cragg-Donald F statistic
        @test Vcov.ranktest!(copy(X), copy(Z), copy(Pi), simple(), l, 3) ≈ (n - l - 3) * s_k^2
        perm = [k; 1:(k - 1)]
        for v in (simple(), robust())
            @test Vcov.ranktest!(copy(X), copy(Z), copy(Pi), v, 0, 0) ≈
                  Vcov.ranktest!(X[:, perm], copy(Z), Pi[:, perm], v, 0, 0)
        end
    end
end
