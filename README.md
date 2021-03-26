[![Build Status](https://travis-ci.com/matthieugomez/Vcov.jl.svg?branch=master)](https://travis-ci.com/matthieugomez/Vcov.jl)

This package should be used as a backend by package developers. 
It allows developers to add a `::CovarianceEstimator` argument in the `fit` method defined by their package. See `FixedEffectModels` for an example.


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

Authors: Matthieu Gomez, Valentin Haddad, Erik Loualiche