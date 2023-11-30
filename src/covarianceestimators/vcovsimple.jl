struct SimpleCovariance <: CovarianceEstimator end

"""
    simple()

Estimate variance-covariance matrix with a simple estimator.
"""
simple() = SimpleCovariance()

Base.show(io::IO, ::SimpleCovariance) =
    print(io, "Simple covariance estimator")

function S_hat(x::RegressionModel, ::SimpleCovariance)
	Symmetric(crossmodelmatrix(x) .* sum(abs2, residuals(x)))
end

function StatsAPI.vcov(x::RegressionModel, ::SimpleCovariance)
    xtx = invcrossmodelmatrix(x)
    Symmetric(xtx .* (sum(abs2, residuals(x)) /  dof_residual(x)))
end

