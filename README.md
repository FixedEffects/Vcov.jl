[![Build status](https://github.com/FixedEffects/Vcov.jl/workflows/CI/badge.svg)](https://github.com/FixedEffects/Vcov.jl/actions)

This package should be used as a backend by package developers. 
It allows developers to add a `::CovarianceEstimator` argument in the `fit` method defined by their package. See `FixedEffectModels` for an example.


Each type defined in this package defines the following methods:
```julia
# return a vector indicating non-missing observations for standard errors
completecases(table, ::CovarianceEstimator) = trues(length(Tables.rows(table)))
# materialize a CovarianceEstimator by using the data needed to compute the standard errors
materialize(table, v::CovarianceEstimator) = v
# return the "meat" of the sandwich estimator
S_hat(x::RegressionModel, ::CovarianceEstimator)
# return variance-covariance matrix
vcov(x::RegressionModel, ::CovarianceEstimator)
# return degrees of freedom for t/F statistics
dof_residual(x::RegressionModel, ::CovarianceEstimator) = dof_residual(x)
```

## Covariance Estimators

- `Vcov.simple()` — iid errors (classical OLS variance)
- `Vcov.robust()` — heteroskedasticity-robust (HC1, default)
- `Vcov.hc(1)`, `Vcov.hc(2)`, `Vcov.hc(3)` — explicit HC variant:
  - **HC1**: White's estimator with small-sample correction *N/(N-K)*
  - **HC2**: residuals scaled by *1/√(1 - hᵢᵢ)*, where *hᵢᵢ* is the leverage
  - **HC3**: residuals scaled by *1/(1 - hᵢᵢ)*
- `Vcov.cluster(:x)` — one-way cluster-robust
- `Vcov.cluster(:x, :y)` — two-way cluster-robust

## References
Kleibergen, F, and Paap, R. (2006) *Generalized reduced rank tests using the singular value decomposition.* Journal of econometrics

Kleibergen, F. and Schaffer, M.  (2007) *RANKTEST: Stata module to test the rank of a matrix using the Kleibergen-Paap rk statistic*. Statistical Software Components, Boston College Department of Economics.
