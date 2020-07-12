
The package define types that inherit from `CovarianceEstimator`. Each type defines the following methods: 
```julia
# take the data needed to compute the standard errors
materialize(df::AbstractDataFrame, v::CovarianceEstimator) = v
# return a vector with false for observations that are missing to compute the standard error
completecases(df::AbstractDataFrame, ::CovarianceEstimator) = trues(size(df, 1))
# return the S_hat vector corresponding to the standard error
S_hat(x::RegressionModel, ::CovarianceEstimator) = error("S_hat not defined for this type")
vcov(x::RegressionModel, ::CovarianceEstimator) = error("vcov not defined for this type")
df_FStat(x::RegressionModel, ::CovarianceEstimator, hasintercept::Bool) = dof_residual(x) - hasintercept
```

The type `RegressionModel` must define four methods: `modelmatrix`, `crossmodelmatrix`, `residuals`, and `dof_residuals`.


