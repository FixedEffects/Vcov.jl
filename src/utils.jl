# The eigenvalues of A are not comparable across coefficients with different units, so test
# for positive semi-definiteness on the correlation scale. A negative variance has a -1 on the
# diagonal there, so tol, which only absorbs roundoff, never hides one.
function pinvertible(A::Symmetric, scale = ones(size(A, 1)), tol = sqrt(eps(real(float(one(eltype(A)))))))
    d = sqrt.(abs.(diag(A)))
    d[iszero.(d)] .= 1
    correlation = Symmetric(Diagonal(d) \ Matrix(A) / Diagonal(d))
    if any(<(-tol), eigvals(correlation))
        @warn "estimated covariance matrix of moment conditions not positive-semidefinite.
                adjustment from Cameron, Gelbach & Miller (2011, JBES) applied.
                 model tests should be interpreted with caution."
        # on the scale of the regressors, so that the adjustment does not depend on their units
        eigval, eigvect = eigen(Symmetric(Diagonal(scale) * Matrix(A) * Diagonal(scale)))
        eigval[eigval .< 0] .= 0
        return Symmetric(Diagonal(scale) \ (eigvect * Diagonal(eigval) * eigvect') / Diagonal(scale))
    else
        return A
    end
end

function regressor_scale(x::RegressionModel)
    scale = vec(sqrt.(sum(abs2, modelmatrix(x); dims = 1)))
    scale[iszero.(scale)] .= 1
    return scale
end
