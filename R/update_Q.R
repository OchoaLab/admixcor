update_Q <- function( ThetaSR, L, R, delta, I, algorithm = c('original', 'nnls', 'bvls', 'glmnet') ) {
    # process options
    algorithm <- match.arg( algorithm )

    if ( algorithm == 'original' ) {
        # Q <- tcrossprod( ThetaSR, R ) %*% MASS::ginv( L ) # no Q regularization
        Q <- tcrossprod( ThetaSR, L %*% R ) %*% MASS::ginv( tcrossprod( L ) + delta * I ) # regularized
    
        # project Q to simplex
        Q <- t( apply( Q, 1L, projsplx ) )
        ## # old projection hack, did not perform as well
        ## for (j in 1:n) {
        ##     Q[j,] <- ifelse( Q[j,] < 0.0, Q[j,] - min(Q[j,]) + 0.0001, Q[j,] )
        ## }
        ## Q <- Q / rowSums(Q)
    } else if ( algorithm == 'nnls' || algorithm == 'bvls' || algorithm == 'glmnet' ) {
        # uses non-negative (nnls) or bounded variable least squares (bvls), which are close to our problem without regularization; glmnet adds regularization; all are improvement over original which has no constraints at all (at first), but none of them directly restrict to simplex so projection is still necessary afterwards
        
        # first set up the problem to resemble a linear regression for Q, of the form Ax=b
        # Q %*% ( L %*% R ) = ThetaSR
        # transpose to best match desired eq setup
        # t( L %*% R ) %*% t( Q ) = t( ThetaSR )
        #B <- t( ThetaSR ) # no need to explicitly transpose it
        # in this case A is the same for all rows
        A <- t( L %*% R )
        # start constructing output Q in parts
        n <- nrow( ThetaSR )
        K <- ncol( ThetaSR )
        Q <- matrix( 0, n, K )
        # solve each row of Q (i.e. column of t(Q)) separately
        for ( i in 1L : n ) {
            # get corresponding row of ThetaSR (column of its transpose)
            b <- ThetaSR[ i, ]
            # solve system of equations with non-negative constraint, or both [0,1] bounds
            if ( algorithm == 'nnls' ) {
                x <- nnls::nnls( A, b )$x
            } else if ( algorithm == 'bvls' ) {
                x <- bvls::bvls( A, b, rep.int( 0, K ), rep.int( 1, K ) )$x
            } else if ( algorithm == 'glmnet' ) {
                # `alpha = 0`: use "ridge" penalty
                # reduced default tolerance of 1e-7 because exact solutions were not good enough
                x <- glmnet::glmnet( A, b, alpha = 0, lambda = delta, intercept = FALSE, lower.limits = 0, upper.limits = 1, thresh = 1e-20 )$beta[,1]
            }
            # now project to simplex (fixes upper bound cap for nnls, only case that doesn't enforce it directly)
            x <- projsplx( x )
            # save column in desired place
            Q[ i, ] <- x
        }

    }

    return( Q )
}
