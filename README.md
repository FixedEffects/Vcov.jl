This package should be used by package developers to compute standard errors. It allows users to specify a `::CovarianceEstimator` argument in the `fit` function of your package. See `FixedEffectModels` for an example.


Each type defined in this package defines the following methods: 
```julia
# return a vector indicating non-missing observations for standard errors
completecases(table, ::CovarianceEstimator) = trues(size(df, 1))
# materialize a CovarianceEstimator by using the data needed to compute the standard errors
materialize(table, v::CovarianceEstimator) = v
# return variance-covariance matrix
vcov(x::RegressionModel, ::CovarianceEstimator) = error("vcov not defined for this type")
# returns the degree of freedom for the F-statistic
df_FStat(x::RegressionModel, ::CovarianceEstimator, hasintercept::Bool) = dof_residual(x) - hasintercept
```

For now, it includes `Vcov.simple()`, `Vcov.robust()`, and `Vcov.cluster(...)`.

