update_Q <- function( ThetaSR, L, R, delta, I ) {
    # update Q (admixture)
    # Q <- tcrossprod( ThetaSR, R ) %*% MASS::ginv( L ) # no Q regularization
    Q <- tcrossprod( ThetaSR, L %*% R ) %*% MASS::ginv( tcrossprod( L ) + delta * I ) # regularized
    
    # project Q to simplex
    Q <- t( apply( Q, 1L, projsplx ) )
    ## # old projection hack, did not perform as well
    ## for (j in 1:n) {
    ##     Q[j,] <- ifelse( Q[j,] < 0.0, Q[j,] - min(Q[j,]) + 0.0001, Q[j,] )
    ## }
    ## Q <- Q / rowSums(Q)

    return( Q )
}
