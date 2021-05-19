@testset "SimpleCovariance" begin
    v = simple()
    @test sprint(show, v) == "Simple covariance estimator"
end

@testset "RobustCovariance" begin
    v = robust()
    @test sprint(show, v) == "Heteroskedasticity-robust covariance estimator"
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

    g1 = group(a1)
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
