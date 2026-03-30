@compile_workload begin
    N, K = 10, 2
    X = randn(N, K)
    Z = randn(N, K)
    e = randn(N)
    XtX = Symmetric(X' * X)
    invXtX = inv(XtX)
    Pi = Z \ X
    for v in (SimpleCovariance(), RobustCovariance(:hc1), RobustCovariance(:hc2), RobustCovariance(:hc3))
        vcovdata = VcovData(X, XtX, invXtX, e, N - K)
        S_hat(vcovdata, v)
        StatsAPI.vcov(vcovdata, v)
    end
    ranktest!(X, Z, Pi, SimpleCovariance(), K, 0)
    ranktest!(X, Z, Pi, RobustCovariance(), K, 0)
end
