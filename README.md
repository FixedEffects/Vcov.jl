
The package define types that inherit from `CovarianceEstimator`. Each type defines the following methods: 
```julia
materialize(df::AbstractDataFrame, v::CovarianceEstimator) = v
completecases(df::AbstractDataFrame, ::CovarianceEstimator) = trues(size(df, 1))
S_hat(x::RegressionModel, ::CovarianceEstimator) = error("S_hat not defined for this type")
df_FStat(x::RegressionModel, ::CovarianceEstimator, hasintercept::Bool) = dof_residual(x) - hasintercept
```

The type `RegressionModel` must define four methods: `modelmatrix`, `crossmodelmatrix`, `residuals`, and `dof_residuals`.




