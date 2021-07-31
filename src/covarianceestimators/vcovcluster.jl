##############################################################################
##
## Cluster Covariance
## 
##############################################################################

struct ClusterCovariance <: CovarianceEstimator
    clusternames::Tuple{Symbol, Vararg{Symbol}} # names are not allowed to be empty
    clusters::Union{Tuple{Vararg{GroupedArray}}, Nothing}
end

"""
    cluster(names::Symbol...)

Estimate variance-covariance matrix with a cluster-robust estimator.

# Arguments
- `names::Symbol...`: column names of variables that represent the clusters.
"""
cluster(ns::Symbol...) = ClusterCovariance((ns...,), nothing)
cluster(ns::NTuple{N, Symbol}, gs::NTuple{N, GroupedArray}) where N =
    ClusterCovariance(ns, gs)

"""
    names(vce::ClusterCovariance)

Return column names of variables used to form clusters for `vce`.
"""
names(v::ClusterCovariance) = v.clusternames

length(v::ClusterCovariance) = length(names(v))

show(io::IO, ::ClusterCovariance) = print(io, "Cluster-robust covariance estimator")

function show(io::IO, ::MIME"text/plain", v::ClusterCovariance)
    print(io, length(v), "-way cluster-robust covariance estimator:")
    for n in names(v)
        print(io, "\n  ", n)
    end
end

function Vcov.completecases(table, v::ClusterCovariance)
    Tables.istable(table) || throw(ArgumentError("completecases requires a table input"))
    N = length(Tables.rows(table))
    out = trues(N)
    aux = BitVector(undef, N)
    columns = Tables.columns(table)
    for name in names(v)
        col = Tables.getcolumn(columns, name)
        if Missing <: eltype(col)
            # The use of aux follows DataFrames.completecases for performance reasons
            # See https://github.com/JuliaData/DataFrames.jl/pull/2726
            aux .= .!ismissing.(col)
            out .&= aux
        end
    end
    return out
end

function materialize(table, v::ClusterCovariance)
    Tables.istable(table) || throw(ArgumentError("materialize requires a table input"))
    ns = names(v)
    return cluster(ns, map(n->GroupedArray(Tables.getcolumn(table, n)), ns))
end

"""
    nclusters(vce::ClusterCovariance)

Return the number of clusters for each dimension/way of clustering.
"""
function nclusters(v::ClusterCovariance)
    NamedTuple{names(v)}(map(x -> x.ngroups, v.clusters))
end

function df_FStat(x::RegressionModel, v::ClusterCovariance, ::Bool)
    minimum(nclusters(v)) - 1
end

function S_hat(x::RegressionModel, v::ClusterCovariance) 
    # Cameron, Gelbach, & Miller (2011): section 2.3
    dim = size(modelmatrix(x), 2) * size(residuals(x), 2)
    S = zeros(dim, dim)
    for c in combinations(1:length(v))
        g = GroupedArray((v.clusters[i] for i in c)...)
        S += (-1)^(length(c) - 1) * helper_cluster(modelmatrix(x), residuals(x), g)
    end
    # scale total vcov estimate by ((N-1)/(N-K)) * (G/(G-1))
    # another option would be to adjust each matrix given by helper_cluster by number of its categories
    # both methods are presented in Cameron, Gelbach and Miller (2011), section 2.3
    # I choose the first option following reghdfe
    G = minimum(nclusters(v))
    rmul!(S, (size(modelmatrix(x), 1) - 1) / dof_residual(x) * G / (G - 1))
end

# res is a Vector in OLS, Matrix in IV
function helper_cluster(X::Matrix, res::Union{Vector, Matrix}, g::GroupedArray)
    X2 = zeros(eltype(X), g.ngroups, size(X, 2) * size(res, 2))
    idx = 0
    for k in 1:size(res, 2)
        for j in 1:size(X, 2)
            idx += 1
            @inbounds @simd for i in 1:size(X, 1)
                X2[g.refs[i], idx] += X[i, j] * res[i, k]
            end
        end
    end
    return Symmetric(X2' * X2)
end

function StatsBase.vcov(x::RegressionModel, v::ClusterCovariance)
    xtx = inv(crossmodelmatrix(x))
    pinvertible(Symmetric(xtx * S_hat(x, v) * xtx))
end

