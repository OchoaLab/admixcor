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
        # uses non-negative or bounded variable least squares, which are equivalent to our problem without regularization
        # solves exactly the requirement that the lower bound be zero, and exactly forces lower triangle to be zero during optimization, but does not set upper bound directly, we set strict cap of 1 as before

        # first set up the problem to resemble a linear regression for L
        # L has to be flattened into a vector
        # this is the value to solve for expressed as a vector, which is easy from matrix:
        b <- as.numeric( tcrossprod( ThetaSR, R ) )
        # the coefficients matrix are from the admixture proportions, but they have to be arranged in a very weird way to match flat L
        # NOTE: for glmnet we could use sparse matrix!
        n <- nrow( Q )
        K <- ncol( Q )
        A <- matrix( 0, nrow = n * K, ncol = K * ( K + 1L ) / 2L )
        for ( u in 1L : K ) {
            # place a copy of a subset of Q in desired location
            # it's hard to explain here why these formulas work, definitely need a picture for that!
            indexes_rows <- ( 1L : n ) + ( u - 1L ) * n
            indexes_cols <- ( 1L : u ) + u * ( u - 1L ) / 2L
            A[ indexes_rows, indexes_cols ] <- Q[ , 1L : u ]
        }
        # solve system of equations with non-negative constraint, or both [0,1] bounds
        # this L is flat one
        if ( algorithm == 'nnls' ) {
            x <- nnls::nnls( A, b )$x
            # only this case requires additional crude enforcement of upper bound (not good enough in some cases!)
            x[ x > 1 ] <- 1
        } else if ( algorithm == 'bvls' ) {
            x <- bvls::bvls( A, b, rep.int( 0, ncol(A) ), rep.int( 1, ncol(A) ) )$x
        } else if ( algorithm == 'glmnet' ) {
            # `alpha = 0`: use "ridge" penalty
            # reduced default tolerance of 1e-7 because exact solutions were not good enough
            x <- glmnet::glmnet( A, b, alpha = 0, lambda = gamma, intercept = FALSE, lower.limits = 0, upper.limits = 1, thresh = 1e-20 )$beta[,1]
        }
        # unflatten L!
        L <- matrix( 0, K, K )
        # this setup was confirmed to work in toy examples (fills by columns too, as desired)
        L[ upper.tri( L, diag = TRUE ) ] <- x
    }
    
    return( L )
}

