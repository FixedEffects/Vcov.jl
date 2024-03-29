function pinvertible(A::Symmetric, tol = eps(real(float(one(eltype(A))))))
    eigval, eigvect = eigen(A)
    small = eigval .<= tol
    if any(small)
        @warn "estimated covariance matrix of moment conditions not positive-semidefinite.
                adjustment from Cameron, Gelbach & Miller (2011, JBES) applied.
                 model tests should be interpreted with caution."
        eigval[small] .= 0
        return Symmetric(eigvect * Diagonal(eigval) * eigvect')
    else
        return A
    end
end