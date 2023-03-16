update_L <- function( ThetaSR, Q, R, gamma, I ) {
    # update L (square root of Psi)
    L <- tcrossprod( MASS::ginv( crossprod( Q ) + gamma * I ) %*% crossprod( Q, ThetaSR ), R ) # regularized L
    #L <- tcrossprod( MASS::ginv( Q ) %*% ThetaSR, R ) # unregularized L
    
    # project L to non-negative Cholesky space
    L[ lower.tri( L ) ] <- 0
    L[ L < 0 ] <- 0
    L[ L > 1 ] <- 1
    
    return( L )
}

