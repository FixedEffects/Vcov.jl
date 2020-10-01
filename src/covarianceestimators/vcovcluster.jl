
##############################################################################
##
## group combines multiple refs
## Missings have a ref of 0
## 
##############################################################################

mutable struct GroupedArray{N} <: AbstractArray{UInt32, N}
    refs::Array{UInt32, N}   # refs must be between 0 and n. 0 means missing
    n::Int                   # Number of potential values (= maximum(refs))
end
group(xs::GroupedArray) = xs

function group(xs::AbstractArray)
    refs = Array{UInt32}(undef, size(xs))
    invpool = Dict{eltype(xs), UInt32}()
    n = 0
    has_missing = false
    @inbounds for i in eachindex(xs)
        x = xs[i]
        if x === missing
            refs[i] = 0
            has_missing = true
        else
            lbl = get(invpool, x, 0)
            if lbl !== 0
                refs[i] = lbl
            else
                n += 1
                refs[i] = n
                invpool[x] = n
            end
        end
    end
    return GroupedArray{ndims(xs)}(refs, n)
end

function group(args...)
    g1 = deepcopy(group(args[1]))
    for j in 2:length(args)
        gj = group(args[j])
        length(g1.refs) == length(gj.refs) || throw(DimensionError())
        combine!(g1, gj)
    end
    factorize!(g1)
end

function combine!(g1::GroupedArray, g2::GroupedArray)
    @inbounds for i in eachindex(g1.refs, g2.refs)
        # if previous one is missing or this one is missing, set to missing
        g1.refs[i] = (g1.refs[i] == 0 || g2.refs[i] == 0) ? 0 : g1.refs[i] + (g2.refs[i] - 1) * g1.n
    end
    g1.n = g1.n * g2.n
    return g1
end

function factorize!(x::GroupedArray{N}) where {N}
    refs = x.refs
    uu = sort!(unique(refs))
    has_missing = uu[1] == 0
    ngroups = length(uu) - has_missing
    dict = Dict{UInt32, UInt32}(zip(uu, UInt32(1-has_missing):UInt32(ngroups)))
    @inbounds for i in eachindex(refs)
        refs[i] = dict[refs[i]]
    end
    GroupedArray{N}(refs, ngroups)
end

##############################################################################
##
## Cluster Covariance
## 
##############################################################################

struct ClusterCovariance <: CovarianceEstimator
    clusters
end

cluster(x::Symbol) = ClusterCovariance((x,))
cluster(args...) = ClusterCovariance(args)


function Vcov.completecases(table, v::ClusterCovariance)
    Tables.istable(table) || throw(ArgumentError("completecases requires a table input"))
    out = trues(length(Tables.rows(table)))
    columns = Tables.columns(table)
    for name in v.clusters
        out .&= .!ismissing.(Tables.getcolumn(columns, name))
    end
    return out
end

function materialize(table, v::ClusterCovariance)
    Tables.istable(table) || throw(ArgumentError("completecases requires a table input"))
    columns = Tables.columns(table)
    ClusterCovariance(
        NamedTuple{v.clusters}(
        ntuple(i -> group(Tables.getcolumn(columns, v.clusters[i])), length(v.clusters))
        ))
end

function nclusters(v::ClusterCovariance)
    map(x -> x.n, v.clusters)
end

function df_FStat(x::RegressionModel, v::ClusterCovariance, ::Bool)
    minimum(nclusters(v)) - 1
end

function S_hat(x::RegressionModel, v::ClusterCovariance) 
    # Cameron, Gelbach, & Miller (2011): section 2.3
    dim = size(modelmatrix(x), 2) * size(residuals(x), 2)
    S = zeros(dim, dim)
    for c in combinations(keys(v.clusters))
        g = group((v.clusters[var] for var in c)...)
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
    X2 = zeros(eltype(X), g.n, size(X, 2) * size(res, 2))
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

