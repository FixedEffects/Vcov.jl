function _precompile_()
    Base.precompile(Tuple{typeof(Vcov.ranktest!),Array{Float64,2},Array{Float64,2},Array{Float64,2},Vcov.SimpleCovariance,Int64,Int64})
    Base.precompile(Tuple{typeof(Vcov.ranktest!),Array{Float64,2},Array{Float64,2},Array{Float64,2},Vcov.ClusterCovariance,Int64,Int64})
    Base.precompile(Tuple{typeof(Vcov.ranktest!),Array{Float64,2},Array{Float64,2},Array{Float64,2},Vcov.RobustCovariance,Int64,Int64})
end
