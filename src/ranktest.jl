##############################################################################
##
## Generalized reduced rank tests using the singularvalue decomposition
##
## See also the Stata command ranktest
## RANKTEST: Stata module to test the rank of a matrix using the Kleibergen-Paap rk statistic
## Authors: Frank Kleibergen, Mark E Schaffer
## More precisely, it corresponds to the Stata command:  ranktest  (X) (Z), wald full

## For Kleibergen-Paap first-stage F-statistics
## See "Generalized reduced rank tests using the singular value decomposition"
## Journal of Econometrics, 2006, 133:1
## Frank Kleibergen and Richard Paap
## doi:10.1016/j.jeconom.2005.02.011
##############################################################################

function ranktest!(X::Matrix{Float64}, 
                    Z::Matrix{Float64}, 
                    Pi::Matrix{Float64}, 
                    vcov_method::CovarianceEstimator, 
                    dof_small::Int, 
                    dof_fes::Int)
    
    if (size(X, 2) == 0) | (size(Z, 2) == 0)
        return NaN
    end

    # Compute theta
    Fmatrix = cholesky!(Symmetric(Z' * Z), check = false).U
    Gmatrix = cholesky!(Symmetric(X' * X), check = false).U
    theta = Fmatrix * (Gmatrix' \ Pi')'

    # compute lambda
    svddecomposition = svd(theta, full = true) 
    u = svddecomposition.U
    vt = svddecomposition.Vt

    k = size(X, 2) 
    l = size(Z, 2) 

    u_sub = u[k:l, k:l]
    vt_sub = vt[k,k]


    # see p.102 of KP
    # the operation u_sup \ sqrt(u_sup * u_sub') may fail. For now, we handle the case k = l but not the more general one
    a_qq = iszero(u_sub) ? u[1:l, k:l] : u[1:l, k:l]  * (u_sub \ sqrt(u_sub * u_sub'))
    b_qq = iszero(vt_sub) ? vt[1:k, k]' : sqrt(vt_sub * vt_sub') * (vt_sub' \ vt[1:k, k]')

    kronv = kron(b_qq, a_qq')
    lambda = kronv * vec(theta)

    # compute vhat
    if vcov_method isa Vcov.SimpleCovariance
        vlab = cholesky!(Hermitian((kronv * kronv') ./ size(X, 1)), check = false)
    else
        K = kron(Gmatrix, Fmatrix)'
        vcovmodel = Vcov.VcovData(Z, K, nothing, X, size(Z, 1) - dof_small - dof_fes) 
        matrix_vcov2 = Vcov.S_hat(vcovmodel, vcov_method)
        vhat = K \ (K \ matrix_vcov2)'
        vlab = cholesky!(Hermitian(kronv * vhat * kronv'), check = false)
    end
    r_kp = lambda' * (vlab \ lambda)
    return r_kp[1]
end
