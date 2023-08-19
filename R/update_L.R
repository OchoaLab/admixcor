# update L (square root of Psi)
update_L <- function( ThetaSR, Q, R, gamma = 0, I = NULL, algorithm = c('original', 'nnls', 'bvls', 'glmnet') ) {
    # process options
    algorithm <- match.arg( algorithm )

    if ( algorithm == 'original' ) {
        # solves using unconstrained ridge regression, then projects to desired space by truncation
        # does not force lower triangle to be zero during optimization, only afterwards by force
        
        L <- tcrossprod( MASS::ginv( crossprod( Q ) + gamma * I ) %*% crossprod( Q, ThetaSR ), R ) # regularized L
        #L <- tcrossprod( MASS::ginv( Q ) %*% ThetaSR, R ) # unregularized L
        
        # project L to non-negative Cholesky space
        L[ lower.tri( L ) ] <- 0
        L[ L < 0 ] <- 0
        L[ L > 1 ] <- 1

        # alpha = 1 / sqrt( max( L %*% t( L ) ) )
    } else if ( algorithm == 'nnls' || algorithm == 'bvls' || algorithm == 'glmnet' ) {
        # uses non-negative (nnls) or bounded variable least squares (bvls), which are equivalent to our problem without regularization; glmnet in theory solves our problem with regularization!
        # nnls solves exactly the requirement that the lower bound be zero, and exactly forces lower triangle to be zero during optimization, but does not set upper bound directly, we set strict cap of 1 as before
        
        # first set up the problem to resemble a linear regression for L, of the form Ax=b
        B <- tcrossprod( ThetaSR, R )
        # start constructing the output L
        K <- ncol( Q )
        L <- matrix( 0, K, K )
        # solve each column of L separately
        for ( i in 1L : K ) {
            # subset the matrices of interest
            # make sure A is a matrix even for i==1L; b is always a vector
            A <- Q[ , 1L : i, drop = FALSE ]
            b <- B[ , i ]
            # solve system of equations with non-negative constraint, or both [0,1] bounds
            if ( algorithm == 'nnls' ) {
                x <- nnls::nnls( A, b )$x
                # only this case requires additional crude enforcement of upper bound (not good enough in some cases!)
                x[ x > 1 ] <- 1
            } else if ( algorithm == 'bvls' ) {
                x <- bvls::bvls( A, b, rep.int( 0, i ), rep.int( 1, i ) )$x
            } else if ( algorithm == 'glmnet' ) {

                if ( i == 1L ) {
                    # glmnet doesn't work to solve a single variable, super annoying!
                    # I wrote my own solver for that nearly trivial case
                    x <- glmnet_one_bounded( A, b, gamma )
                } else {
                    # `alpha = 0`: use "ridge" penalty
                    # reduced default tolerance of 1e-7 because exact solutions were not good enough
                    x <- glmnet::glmnet( A, b, alpha = 0, lambda = gamma, intercept = FALSE, lower.limits = 0, upper.limits = 1, thresh = 1e-20 )$beta[,1]
                }
            }
            # save column in desired place
            L[ 1L : i , i ] <- x
        }
    }
    
    return( L )
}

glmnet_one_bounded <- function( a, b, gamma ) {
    # this is very silly, but glmnet won't work if there's a single variable to solve for
    # so let's do it explicitly
    # this is objective, for reference:
    # O = ( a * x - b ) %*% t( a * x - b ) + gamma * x^2
    # this is unbounded solution
    #x <- ( b %*% a ) / ( crossprod( a ) + gamma )
    x <- crossprod( a, b ) / ( crossprod( a ) + gamma )

    # if solution is in bounds, that's the desired solution!
    if ( 0 < x && x < 1 )
        return( x )
    
    # else, if solution is out of bounds, for this single variable case only other solutions are edge cases, i.e. 0 and 1
    obj0 <- crossprod( b )
    obj1 <- crossprod( a - b ) + gamma
    return( if ( obj0 > obj1 ) 0 else 1 )
}
