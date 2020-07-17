This package should be used by package developers to compute standard errors. The goal is to allow users to specify a `::CovarianceEstimator` argument in the `fit` function of your package. See `FixedEffectModels` for an example.


Each type defined in this package defines the following methods: 
```julia
# take the data needed to compute the standard errors
materialize(table, v::CovarianceEstimator) = v
# return a vector with false for observations that are missing to compute the standard error
completecases(table, ::CovarianceEstimator) = trues(size(df, 1))
# return variance-covariance matrix
vcov(x::RegressionModel, ::CovarianceEstimator) = error("vcov not defined for this type")
# returns the degree of freedom for the F-statistic
df_FStat(x::RegressionModel, ::CovarianceEstimator, hasintercept::Bool) = dof_residual(x) - hasintercept
```

