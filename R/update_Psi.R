# NOTE: alpha must be zero for bvls version
update_Psi <- function( Theta, Q, alpha = 0, algorithm = c('glmnet', 'bvls') ) {
    # process options
    algorithm <- match.arg( algorithm )
    
    # get dimensions
    K <- ncol( Q )

    # start by formulating the ordinary Kronecker trick
    b <- as.vector( Theta )
    # this is the non-trivial part
    A <- Q %x% Q

    # the rest is to not double-state/count the off-diagonal, which I suspect could result in singularity, though regardless shrinks the problem by almost a factor of 2 so its bound to result in faster and/or better solutions

    # make a fake Psi to pass to lower.tri only
    # actually this set that results in complement indexes is more useful
    Psi <- t( matrix( 1 : K^2, K, K ) )
    # for the coefficients of A that get folded to match subsetted x, we have to to extra work because A itself isn't symmetric!  We have to sum those coefficients instead of just doubling them
    indexes_rows <- as.vector( lower.tri( Theta, diag = TRUE ) )
    indexes_cols <- as.vector( lower.tri( Psi, diag = TRUE ) )
    indexes_cols_nodiag <- which( as.vector( lower.tri( Psi ) ) )
    # hacky way to find out where the paired columns are, for before subsetting
    indexes_cols_nodiag_comp <- Psi[ lower.tri( Psi ) ]
    # subset all appropriately!
    b <- b[ indexes_rows ]
    # for A, do it in two parts
    # first throw away rows, very easy and speeds up second step
    A <- A[ indexes_rows, ]
    # sum the two different coefficients of x (paired columns), then subset columns
    A[ , indexes_cols_nodiag ] <- A[ , indexes_cols_nodiag ] + A[ , indexes_cols_nodiag_comp ]
    A <- A[ , indexes_cols ]

    # now get solutions!
    # solve system of equations with non-negative constraint, or both [0,1] bounds
    if ( algorithm == 'bvls' ) {
        if ( alpha != 0 )
            stop( '`L_algorithm == "bvls"` requires `alpha = 0`' )
        # length of x, needed to set bounds
        lx <- K * ( K + 1 ) / 2
        x <- bvls::bvls( A, b, rep.int( 0, lx ), rep.int( 1, lx ) )$x
    } else if ( algorithm == 'glmnet' ) {
        if ( K == 1L ) {
            # glmnet doesn't work to solve a single variable, super annoying!
            # I wrote my own solver for that nearly trivial case
            x <- glmnet_one_bounded( A, b, alpha )
        } else {
            # `alpha = 0`: use "ridge" penalty
            # reduced default tolerance of 1e-7 because exact solutions were not good enough
            x <- glmnet::glmnet( A, b, alpha = 0, lambda = alpha, intercept = FALSE, lower.limits = 0, upper.limits = 1, thresh = 1e-20 )$beta[,1]
        }
    }
    # lastly, reform x into Psi, return that!
    # first copy lower triangle and diagonal from linearized version
    Psi[ lower.tri( Psi, diag = TRUE ) ] <- x
    # then set upper triangle too
    Psi[ upper.tri( Psi ) ] <- t( Psi )[ upper.tri( Psi ) ]
    
    return( Psi )
}
