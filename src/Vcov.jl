module Vcov

using Combinatorics: combinations
using GroupedArrays: GroupedArray
using LinearAlgebra: cholesky!, Symmetric, Hermitian, svd, rmul!, eigen, Diagonal
using StatsAPI: StatsAPI, RegressionModel, modelmatrix, crossmodelmatrix, residuals, dof_residual
using StatsBase: CovarianceEstimator
using Tables: Tables
using Base: @propagate_inbounds

##############################################################################
##
## Mimimum RegressionModel used in Vcov
##
##############################################################################
struct VcovData{T, Tu, N} <: RegressionModel
    modelmatrix::Matrix{Float64} # X
    crossmodelmatrix::T          # X'X 
    invcrossmodelmatrix::Tu      # inv(X'X)
    residuals::Array{Float64, N} # vector or matrix of residuals (matrix in the case of IV, residuals of Xendo on (Z, Xexo))
    dof_residual::Int
end
StatsAPI.modelmatrix(x::VcovData) = x.modelmatrix
StatsAPI.crossmodelmatrix(x::VcovData) = x.crossmodelmatrix
invcrossmodelmatrix(x::VcovData) = x.invcrossmodelmatrix
StatsAPI.residuals(x::VcovData) = x.residuals
StatsAPI.dof_residual(x::VcovData) = x.dof_residual

##############################################################################
##
## Any type used for standard errors must define the following methods:
##
##############################################################################
function completecases(table, ::CovarianceEstimator)
	Tables.istable(table) || throw(ArgumentError("completecases requires a table input"))
	return trues(length(Tables.rows(table)))
end
materialize(table, v::CovarianceEstimator) = v
S_hat(x::RegressionModel, ::CovarianceEstimator) = error("S_hat not defined for this type")
StatsAPI.dof_residual(x::RegressionModel, ::CovarianceEstimator) = dof_residual(x)


include("utils.jl")
include("covarianceestimators/vcovsimple.jl")
include("covarianceestimators/vcovrobust.jl")
include("covarianceestimators/vcovcluster.jl")
include("ranktest.jl")

include("precompile.jl")
_precompile_()

end
