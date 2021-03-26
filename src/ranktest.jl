##############################################################################
##
## Generalized reduced rank tests using the singularvalue decomposition
## Kleibergen Paap, Journal of Econometrics (2006)
##
## See also the Stata command ranktest
## RANKTEST: Stata module to test the rank of a matrix using the Kleibergen-Paap rk statistic
## Authors: Frank Kleibergen, Mark E Schaffer
##
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
                    df_small::Int, 
                    df_absorb::Int)
    
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

    if k==l
        # see p.102 of KP
        a_qq = u_sub[1] >= 0 ? u[1:l, k:l] : -u[1:l, k:l] 
        b_qq = vt_sub[1] >= 0 ? vt[1:k, k]' : -vt[1:k, k]'
    else
        # there might be something do to about the sign here as well
        a_qq = u[1:l, k:l]  * (u_sub \ sqrt(u_sub * u_sub'))
        b_qq = sqrt(vt_sub * vt_sub') * (vt_sub' \ vt[1:k, k]')
    end 

    kronv = kron(b_qq, a_qq')
    lambda = kronv * vec(theta)

    # compute vhat
    if vcov_method isa Vcov.SimpleCovariance
        vlab = cholesky!(Hermitian((kronv * kronv') ./ size(X, 1)), check = false)
    else
        K = kron(Gmatrix, Fmatrix)'
        vcovmodel = Vcov.VcovData(Z, K, X, size(Z, 1) - df_small - df_absorb) 
        matrix_vcov2 = Vcov.S_hat(vcovmodel, vcov_method)
        vhat = K \ (K \ matrix_vcov2)'
        vlab = cholesky!(Hermitian(kronv * vhat * kronv'), check = false)
    end
    r_kp = lambda' * (vlab \ lambda)
    return r_kp[1]
end
