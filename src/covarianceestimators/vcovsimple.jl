struct SimpleCovariance <: CovarianceEstimator end

"""
    simple()

Estimate variance-covariance matrix with a simple estimator.
"""
simple() = SimpleCovariance()

Base.show(io::IO, ::SimpleCovariance) =
    print(io, "Simple covariance estimator")

function S_hat(x::RegressionModel, ::SimpleCovariance)
	rmul!(crossmodelmatrix(x), sum(abs2, residuals(x)))
end

function StatsAPI.vcov(x::RegressionModel, ::SimpleCovariance)
    invcrossmodelmatrix = Matrix(inv(crossmodelmatrix(x)))
    rmul!(invcrossmodelmatrix, sum(abs2, residuals(x)) /  dof_residual(x))
    Symmetric(invcrossmodelmatrix)
end

