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
    g1 = group(a1)
    c1 = cluster((:a,), (g1,))
    @test nclusters(c1) == (a=N,)
end
