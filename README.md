

This package acts on RegressionModels. It uses four methods: `modelmatrix`, `crossmodelmatrix`, `residuals`, and `dof_residuals`.


All the standard-errors types defined in this package create the following methods: 
```
materialize(df::AbstractDataFrame, v::CovarianceEstimator) = v
completecases(df::AbstractDataFrame, ::CovarianceEstimator) = trues(size(df, 1))
S_hat(x::RegressionModel, ::CovarianceEstimator) = error("S_hat not defined for this type")
df_FStat(x::RegressionModel, ::CovarianceEstimator, hasintercept::Bool) = dof_residual(x) - hasintercept
```


