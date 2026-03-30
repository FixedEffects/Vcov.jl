struct RobustCovariance <: CovarianceEstimator
    correction::Symbol # :hc1, :hc2, :hc3
    RobustCovariance(correction::Symbol = :hc1) = new(correction)
end

"""
    robust()

Estimate variance-covariance matrix with a heteroskedasticity-robust (HC1) estimator.
Use `hc(n)` to select a specific variant (HC1, HC2, or HC3).
"""
robust() = RobustCovariance(:hc1)

"""
    hc(n::Integer)

Estimate variance-covariance matrix with an HC`n` heteroskedasticity-robust estimator.
- `hc(1)`: White's estimator with small-sample correction N/(N-K).
- `hc(2)`: residuals divided by sqrt(1 - hᵢᵢ), where hᵢᵢ is the leverage.
- `hc(3)`: residuals divided by (1 - hᵢᵢ), where hᵢᵢ is the leverage.
"""
function hc(n::Integer)
    n in 1:3 || throw(ArgumentError("HC$n is not supported. Use hc(1), hc(2), or hc(3)."))
    RobustCovariance(Symbol("hc$n"))
end

function Base.show(io::IO, v::RobustCovariance)
    print(io, "Heteroskedasticity-robust ($(uppercase(string(v.correction)))) covariance estimator")
end

# Compute leverage values hᵢᵢ = xᵢ'(X'X)⁻¹xᵢ
function leverage(m::Matrix, invXtX)
    P = m * invXtX
    vec(sum(m .* P, dims=2))
end

function S_hat(x::RegressionModel, v::RobustCovariance)
    m = modelmatrix(x)
    r = residuals(x)
    if v.correction == :hc1
        if r isa AbstractVector
            X2 = m .* r
        else
            # IV case: columns are m[:,j] .* r[:,k] for each (j,k)
            X2 = zeros(size(m, 1), size(m, 2) * size(r, 2))
            idx = 0
            for k in 1:size(r, 2)
                for j in 1:size(m, 2)
                    idx += 1
                    @views X2[:, idx] .= m[:, j] .* r[:, k]
                end
            end
        end
        S2 = X2' * X2
        return Symmetric(rmul!(S2, size(m, 1) / dof_residual(x)))
    elseif v.correction == :hc2
        r isa AbstractVector || throw(ArgumentError("HC2 is only supported for OLS (vector residuals)"))
        h = leverage(m, invcrossmodelmatrix(x))
        X2 = m .* (r ./ sqrt.(1 .- h))
        return Symmetric(X2' * X2)
    elseif v.correction == :hc3
        r isa AbstractVector || throw(ArgumentError("HC3 is only supported for OLS (vector residuals)"))
        h = leverage(m, invcrossmodelmatrix(x))
        X2 = m .* (r ./ (1 .- h))
        return Symmetric(X2' * X2)
    else
        throw(ArgumentError("Unknown correction type: $(v.correction). Use :hc1, :hc2, or :hc3."))
    end
end

function StatsAPI.vcov(x::RegressionModel, v::RobustCovariance)
    xtx = invcrossmodelmatrix(x)
    pinvertible(Symmetric(xtx * S_hat(x, v) * xtx))
end
