
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
Base.size(g::GroupedArray) = size(g.refs)
@propagate_inbounds Base.getindex(g::GroupedArray, i::Int) = getindex(g.refs, i)
@propagate_inbounds Base.getindex(g::GroupedArray, I) = getindex(g.refs, I)

group(xs::GroupedArray) = xs

function group(xs::AbstractArray)
    _group(DataAPI.refarray(xs), DataAPI.refpool(xs))
end


function _group(xs, ::Nothing)
    refs = Array{UInt32}(undef, size(xs))
    invpool = Dict{eltype(xs), UInt32}()
    n = UInt32(0)
    z = UInt32(0)
    @inbounds for i in eachindex(xs)
        x = xs[i]
        if x === missing
            refs[i] = 0
        else
            lbl = get(invpool, x, z)
            if !iszero(lbl)
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

function _group(ra, rp)
    refs = Array{UInt32}(undef, size(ra))
    fira = firstindex(ra)
    firp = firstindex(rp)
    hashes = Array{UInt32}(undef, length(rp))
    n = 0
    for i in eachindex(hashes)
        if rp[i+firp-1] === missing
            hashes[i] = UInt32(0)
        else
            n += 1
            hashes[i] = n
        end
    end
    @inbounds for i in eachindex(refs)
        refs[i] = hashes[ra[i+fira-1]-firp+1]
    end
    return GroupedArray{ndims(ra)}(refs, n)
end

function group(args...)
    g1 = deepcopy(group(args[1]))
    for j in 2:length(args)
        gj = group(args[j])
        size(g1) == size(gj) || throw(DimensionMismatch(
            "cannot match array of size $(size(g1)) with array of size $(size(gj))"))
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

# An in-place version of group() that relabels the refs
function factorize!(g::GroupedArray{N}) where {N}
    refs = g.refs
    invpool = Dict{UInt32, UInt32}()
    n = 0
    z = UInt32(0)
    @inbounds for i in eachindex(refs)
        x = refs[i]
        if !iszero(x)
            lbl = get(invpool, x, z)
            if !iszero(lbl)
                refs[i] = lbl
            else
                n += 1
                refs[i] = n
                invpool[x] = n
            end
        end
    end
    return GroupedArray{N}(refs, n)
end

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
    out = trues(length(Tables.rows(table)))
    columns = Tables.columns(table)
    for name in names(v)
        out .&= .!ismissing.(Tables.getcolumn(columns, name))
    end
    return out
end

function materialize(table, v::ClusterCovariance)
    Tables.istable(table) || throw(ArgumentError("materialize requires a table input"))
    ns = names(v)
    return cluster(ns, map(n->group(Tables.getcolumn(table, n)), ns))
end

"""
    nclusters(vce::ClusterCovariance)

Return the number of clusters for each dimension/way of clustering.
"""
function nclusters(v::ClusterCovariance)
    NamedTuple{names(v)}(map(x -> x.n, v.clusters))
end

function df_FStat(x::RegressionModel, v::ClusterCovariance, ::Bool)
    minimum(nclusters(v)) - 1
end

function S_hat(x::RegressionModel, v::ClusterCovariance) 
    # Cameron, Gelbach, & Miller (2011): section 2.3
    dim = size(modelmatrix(x), 2) * size(residuals(x), 2)
    S = zeros(dim, dim)
    for c in combinations(1:length(v))
        g = group((v.clusters[i] for i in c)...)
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

